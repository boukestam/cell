import { Workbook } from "exceljs";
import GenbankParser from "genbank-parser";
import { ParseConfig, parse } from 'papaparse';
import { Datasheet } from "./Datasheet";
import { XMLParser } from "fast-xml-parser";
import { SBML, parseSBML } from "./SBML";

export async function loadExcel(url: string, sheetName: string, header: boolean, skipRows: number = 0) {
  const request = await fetch(url);
  const buffer = await request.arrayBuffer();

  const workbook = new Workbook();
  await workbook.xlsx.load(buffer);

  const sheet = workbook.getWorksheet(sheetName);

  const rows = sheet.getRows(1 + skipRows, sheet.rowCount - skipRows).map(row =>
    (row.values as any[]).slice(1).map(value => value?.result || value)
  );

  return new Datasheet(rows.slice(1), rows[0]);
}

export async function loadGenbank(url: string) {
  const request = await fetch(url);
  const text = await request.text();

  return GenbankParser(text)[0];
}

export async function loadCSV(url: string, header: boolean, columns?: string[]): Promise<Datasheet> {
  if (header && columns) throw new Error("Cannot specify columns when header is true");

  try {
    const request = await fetch(url);
    let text = await request.text();

    if (columns) {
      text = columns.join(',') + '\n' + text;
    }

    const config: ParseConfig = { header: header, skipEmptyLines: true, delimiter: ',' };

    const result = parse(text, config);

    if (header && !columns) columns = result.meta.fields;

    const rows = columns ? result.data.map(row => columns.map(column => row[column])) : result.data;

    return new Datasheet(rows, columns);
  } catch (e) {
    console.error(`Error loading CSV from ${url}`);
    throw e;
  }
}

export type TSV = Map<string, Datasheet>;

export async function loadTSV(url: string): Promise<TSV> {
  const request = await fetch(url);
  const text = await request.text();

  const sheets = new Map<string, Datasheet>();

  const lines = text.split('\n');

  let sheetName = '';
  let rows: string[][] = [];
  let columns: string[] = [];

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];

    if (line.startsWith("!!SBtab")) {
      // ... TableName='Sheet1' ...
      sheetName = line.match(/TableName='(.*?)'/)?.[1];

      // !Column1 !Column2 !Column3
      columns = lines[i + 1].split('\t').map(column => column.slice(1));

      i += 1;
    } else if (line.trim().length > 0) {
      // Value1 Value2 Value3
      rows.push(line.split('\t'));
    } else {
      if (sheetName) {
        sheets.set(sheetName, new Datasheet(rows, columns));
      }

      sheetName = '';
      rows = [];
      columns = [];
    }
  }

  return sheets;
}

export async function loadSBML(url: string): Promise<SBML> {
  const request = await fetch(url);
  const text = await request.text();

  return parseSBML(text);
}