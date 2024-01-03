// Define how to calculate transcription rate constants as in equation 3 for transcription reactions.
// Uses mature transcript length and proteomics for promoter strength.

import { getProteinCount } from "./Data";
import { rnaPolKcat, baseMap, rnaPolKd, ATPconc, CTPconc, UTPconc, GTPconc, rnaPolK0, RnaPconc, riboKcat, trnaPolKcat, rrnaPolKcat, CytoVolume, riboKd, riboK0, ribosomeConc } from "./Globals";
import { aaMap, aaTRNAMap } from "./RNAMaps";

function checkBases(sequence: string, baseMap: Map<string, any>) {
  for (const base of new Set(sequence)) {
    if (!baseMap.has(base)) {
      console.error(sequence, baseMap);
      throw new Error("Unknown base(s) in RNA sequence: " + base);
    }
  }
}

export function transcriptRate(rnaMetID: string, ptnMetID: string, rnasequence: string, jcvi2ID: string) {
  // Add trascription reaction

  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const { proteinCount, proteinName } = getProteinCount(ptnMetID, jcvi2ID);

  const kcat_mod = Math.min(rnaPolKcat * (proteinCount / (180)), 2 * 90);

  // Add total number of monomers to parameter dict

  const CMono1 = baseMap.get(rnasequence[0])!;

  const CMono2 = baseMap.get(rnasequence[1])!;

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const NMono_A = baseCount.get("A")!;

  const NMono_U = baseCount.get("C")!;

  const NMono_C = baseCount.get("G")!;

  const NMono_G = baseCount.get("U")!;

  const NMonoDict = [NMono_A, NMono_C, NMono_G, NMono_U];

  const NMonoSum = NMono_A * rnaPolKd / ATPconc + NMono_C * rnaPolKd / CTPconc + NMono_U * rnaPolKd / UTPconc + NMono_G * rnaPolKd / GTPconc;

  const k_transcription = kcat_mod / ((1 + rnaPolK0 / RnaPconc) * (rnaPolKd ** 2) / (CMono1 * CMono2) + NMonoSum + n_tot - 1);

  return k_transcription;
}

// Define transcription rate for ribosomal protein-coding genes.

export function riboTranscriptRate(rnaMetID: string, ptnMetID: string, rnasequence: string, jcvi2ID: string) {
  // Add trascription reaction

  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const { proteinCount, proteinName } = getProteinCount(ptnMetID, jcvi2ID);

  const kcat_mod = Math.min(riboKcat * (proteinCount / (180)), 2 * 90);

  // Add total number of monomers to parameter dict

  const CMono1 = baseMap.get(rnasequence[0])!;

  const CMono2 = baseMap.get(rnasequence[1])!;

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const NMono_A = baseCount.get("A")!;

  const NMono_U = baseCount.get("C")!;

  const NMono_C = baseCount.get("G")!;

  const NMono_G = baseCount.get("U")!;

  const NMonoDict = [NMono_A, NMono_C, NMono_G, NMono_U];

  const NMonoSum = NMono_A * rnaPolKd / ATPconc + NMono_C * rnaPolKd / CTPconc + NMono_U * rnaPolKd / UTPconc + NMono_G * rnaPolKd / GTPconc;

  const k_transcription = kcat_mod / ((1 + rnaPolK0 / RnaPconc) * (rnaPolKd ** 2) / (CMono1 * CMono2) + NMonoSum + n_tot - 1);

  return k_transcription;
}

// Define how to calculate transcription rate constants as in equation 3 for transcription reactions.
// Uses mature transcript length and proteomics for promoter strength.

export function degradationRate(rnaMetID: string, rnasequence: string) {
  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const kcat = (18 / 452) * 88; //1/s # INSTEAD OF 18 or 20

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const k_deg = kcat / n_tot;

  return k_deg;
}

// Define transcription rate for tRNA genes.

export function trnaTranscriptRate(rnaMetID: string, rnasequence: string) {
  // Add trascription reaction

  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const kcat = trnaPolKcat; //1/s # INSTEAD OF 18 or 20

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const k_transcription = kcat / n_tot;

  return k_transcription;
}

// Define transcription rate for rRNA genes.

export function rrnaTranscriptRate(rnaMetID: string, rnasequence: string) {
  // Add trascription reaction

  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const kcat = rrnaPolKcat; //1/s # INSTEAD OF 18 or 20

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const k_transcription = kcat / n_tot;

  return k_transcription;
}

// Convert particle counts to mM concentrations for the ODE Solver

export function partTomM(particles: number, pmap: Map<string, number>) {
  const NA = 6.022e23; // Avogadro's
  const r_cell = 200e-9; // 200 nm radius, 400 nm diameter
  const V = ((4 / 3) * Math.PI * (r_cell) ** 3) * (1000); // for a spherical cell

  const cellVolume = CytoVolume;

  const conc = (particles * 1000.0) / (NA * cellVolume);

  return conc;
}

// Define how to calculate translation rate constants as in equation 3 for translation reactions.

export function translateRate(rnaMetID: string, ptnMetID: string, rnasequence: string, aasequence: string, pmap: Map<string, number>) {
  // Add translation reaction

  // Check that we know all residues used in the sequence
  checkBases(aasequence, aaMap);

  // Count how many times each residue is used
  const aaCount = new Map<string, number>();
  for (const aa of new Set(aasequence)) {
    aaCount.set(aa, aasequence.split(aa).length - 1);
  }

  const NStop = aaCount.get("*")!;

  if (NStop > 1) {
    console.log("EXTRA STOP CODON: MISTAKE IN TRANSLATION");
  }

  let NMonoSum = 0;

  for (const [key, trnaID] of aaTRNAMap) {
    const aaCntPtn = aaCount.get(key) || 0;

    NMonoSum += aaCntPtn * riboKd / partTomM(Math.max(1, pmap.get(trnaID) || 0), pmap);
  }

  const n_tot = [...aaCount.values()].reduce((a, b) => a + b);

  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const transcript_length = [...baseCount.values()].reduce((a, b) => a + b);

  let ribo_num = Math.max(1, Math.round(transcript_length / 125 - 1));
  ribo_num = Math.min(15, ribo_num);

  let kcat_mod
  if (ribo_num > 1) {
    kcat_mod = ribo_num * 0.25 * riboKcat + 0.2 * riboKcat; // 503 ribosomes
  } else {
    kcat_mod = 0.45 * riboKcat;
  }

  const k_translation = kcat_mod / ((1 + riboK0 / ribosomeConc) * (riboKd ** 2) / (partTomM(Math.max(1, pmap.get("M_fmettrna_c")!), pmap) ** 2) + NMonoSum + n_tot - 1);

  return k_translation;
}

// Define how to calculate translation rate constants as in equation 3 for translation reactions.

export function translocRate(aaSequence: string) {
  // Check that we know all residues used in the sequence
  checkBases(aaSequence, aaMap);

  // Count how many times each residue is used
  const aaCount = new Map<string, number>();
  for (const aa of new Set(aaSequence)) {
    aaCount.set(aa, aaSequence.split(aa).length - 1);
  }

  const ptnLen = [...aaCount.values()].reduce((a, b) => a + b);

  const k_transloc = 50 / ptnLen;

  return k_transloc;
}