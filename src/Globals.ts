import { data } from "./data/Data";
import { GeneLocation } from "./data/Gene";

export const AOEtoJ2 = new Map();
export const J2toAOE = new Map();
export const genomePtnLocDict = new Map<string, GeneLocation>();
export const genomeRnaLocDict = new Map<string, GeneLocation>();
export const Locus3A = [];

// Create a dictionary to map protein IDs to their associated metabolite IDs

export const MetPtnGenes = [
  'MMSYN1_0621', 'MMSYN1_0419',
  'MMSYN1_0513', 'MMSYN1_0512', 'MMSYN1_0117', 'MMSYN1_0139',
  'MMSYN1_0147', 'MMSYN1_0304', 'MMSYN1_0420', 'MMSYN1_0616', 'MMSYN1_0617',
  'MMSYN1_0218', 'MMSYN1_0214', 'MMSYN1_0875',
  'MMSYN1_0836', 'MMSYN1_0642', 'MMSYN1_0643', 'MMSYN1_0641'
];

export const MetLocusNums = [];

for (const gene of MetPtnGenes) {
  const locusNum = gene.split('_')[1];
  MetLocusNums.push(locusNum);
}

export const named_PTN_list = [];

// Define Transcription and Translation Rate Constants

export const ModelSpecies = new Set<string>();

export const defaultPtnCount = 10;

// Global parameters for transcription
export const rnaPolKcat = 0.155 * 187 / 493 * 20; // nt/s
export const rnaPolK0 = 1e-4; // mM
export const rnaPolKd = 0.1; // mM

export const rrnaPolKcat = 85 * 2; // nt/s
export const trnaPolKcat = 0.155 * 187 / 493 * 25; // nt/s

export const krnadeg = 0.00578 / 2; // 1/s

// Global parameter for degradation of proteins
// Derived from eLife's model, using average protein half life of 25 hours. 
export const ptnDegRate = 7.70e-6; // 1/s

export const ATPconc = 1.04; // mM
export const UTPconc = 0.68; // mM
export const CTPconc = 0.34; // mM
export const GTPconc = 0.68; // mM

// Cell radius (meters):
export const r_cell = 2.0 * (10 ** -7); // m

export const CytoVolume = (4 * Math.PI / 3) * 1000 * r_cell ** 3; // L
export const cellVolume = CytoVolume;

// Avogadro:
export const avgdr = 6.022e23; // molec/mol
export const Avognum = avgdr;

export const countToMiliMol = 1000 / (avgdr * cellVolume);

export const RnaPconc = 187 * countToMiliMol; // mM

// Global parameter for degradation of mRNAs
export const rnaDegRate = 0.00578 / 2; // 1/s

// Create a map for rna sequence to NTP concentration.
export const baseMap = new Map<string, number>([
  ["A", ATPconc],
  ["U", UTPconc],
  ["G", GTPconc],
  ["C", CTPconc]
]);

// Global parameters for translation
export const riboKcat = 12; // 1/s
export const riboK0 = 4 * 25e-6; // 
export const riboKd = 0.0001; // 

export const ribosomeConc = 503 * countToMiliMol; // mM

// Concentration of charged tRNA
export const ctRNAconc = 150 * countToMiliMol; // mM

// Ratio of volume change to surface area change for 200nm radius cell
export const volToSA = 1.5874;

export function initDicts() {
  for (const f of data.genome2.features) {
    if (f.type === "CDS") {
      const JCVSYN2_tag = f.notes['locus_tag'][0];
      // Not all entries have an AOE protein_id
      if ('protein_id' in f.notes) {
        const AOE_locus = f.notes['protein_id'][0];
        AOEtoJ2.set(AOE_locus, JCVSYN2_tag);
        J2toAOE.set(JCVSYN2_tag, AOE_locus);
      } else {
        console.log("Locus ", JCVSYN2_tag, " has no AOE id!");
      }
    }
    if (f.type === "rRNA") {
      const JCVSYN2_tag = f.notes['locus_tag'][0];
    }
    if (f.type === "tRNA") {
      const JCVSYN2_tag = f.notes['locus_tag'][0];
    }
  }

  for (const f of data.genome3A.features) {
    if (f.type === "CDS") {
      const JCVSYN3A_tag = f.notes['locus_tag'][0];
      Locus3A.push(JCVSYN3A_tag);
      // Not all entries have an AOE protein_id
      if ('protein_id' in f.notes) {
        genomePtnLocDict.set(JCVSYN3A_tag, new GeneLocation(f));
      } else {
        console.log("Locus ", JCVSYN3A_tag, " is pseudo.");
      }
    }
    if (f.type === "rRNA") {
      const JCVSYN3A_tag = f.notes['locus_tag'][0];
      Locus3A.push(JCVSYN3A_tag);
      genomeRnaLocDict.set(JCVSYN3A_tag, new GeneLocation(f));
    }
    if (f.type === "tRNA") {
      const JCVSYN3A_tag = f.notes['locus_tag'][0];
      Locus3A.push(JCVSYN3A_tag);
      genomeRnaLocDict.set(JCVSYN3A_tag, new GeneLocation(f));
    }
  }

  for (const row of data.riboPtnMetDF.rows) {
    named_PTN_list.push(row["gene"]);
  }

  for (const row of data.PtnMetDF.rows) {
    named_PTN_list.push(row["gene"]);
  }
}