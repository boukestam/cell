function isNumber(s: string) {
  // Check if string is a number with a regex
  const regex = /^-?\d+\.?\d*$/;
  return regex.test(s);
}

export class Datasheet {

  rows: any[][];
  columns: string[];

  constructor(rows: any[][], columns?: string[]) {
    this.rows = rows;
    this.columns = columns;

    for (const row of rows) {
      for (let i = 0; i < row.length; i++) {
        const value = row[i];

        if (value === undefined) continue;

        if (typeof value === "number") continue;
        if (typeof value === "boolean") continue;

        // Try to parse numbers
        if (typeof value === 'string' && isNumber(value)) {
          row[i] = parseFloat(value);
          continue;
        }

        // Try to parse booleans
        if (value.toLowerCase() === 'true') {
          row[i] = true;
          continue;
        } else if (value.toLowerCase() === 'false') {
          row[i] = false;
          continue;
        }
      }
    }
  }

  findRowIndex(column: string | number, value: any) {
    const columnIndex = typeof column === "number" ? column : this.columns.indexOf(column);
    if (columnIndex === -1) {
      throw new Error(`Column not found: ${column}`);
    }

    return this.rows.findIndex(row => row[columnIndex] == value);
  }

  findRow(column: string | number, value: any) {
    const rowIndex = this.findRowIndex(column, value);
    if (rowIndex === -1) return undefined;

    return this.getRow(rowIndex);
  }

  getRow(index: number): any {
    if (index < 0 || index >= this.rows.length) throw new Error(`Index out of bounds: ${index}`);

    if (this.columns) {
      return this.columns.reduce((row, column, i) => {
        row[column] = this.rows[index][i];
        row[i] = this.rows[index][i];
        return row;
      }, {});
    } else {
      return this.rows[index];
    }
  }

  getRows() {
    return this.rows.map((row, i) => this.getRow(i));
  }

  getColumn(column: string | number) {
    const columnIndex = typeof column === "number" ? column : this.columns.indexOf(column);
    if (columnIndex === -1) {
      throw new Error(`Column not found: ${column}`);
    }

    return this.rows.map(row => row[columnIndex]);
  }

  setColumn(column: string | number, values: any[]) {
    const columnIndex = typeof column === "number" ? column : this.columns.indexOf(column);
    if (columnIndex === -1) {
      throw new Error(`Column not found: ${column}`);
    }

    for (let i = 0; i < values.length; i++) {
      this.rows[i][columnIndex] = values[i];
    }
  }

  filter(callback: (row: any, index: number) => boolean) {
    const rows = this.getRows().filter(callback);
    return new Datasheet(rows, this.columns);
  }

  map<T>(callback: (row: any, index: number) => T): T[] {
    const rows = this.getRows().map(callback);
    return rows;
  }

  concat(...datasheets: Datasheet[]) {
    const rows = this.getRows();
    for (const datasheet of datasheets) {
      rows.push(...datasheet.getRows());
    }

    return new Datasheet(rows, this.columns);
  }
}