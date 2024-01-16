// Define how to calculate transcription rate constants as in equation 3 for transcription reactions.
// Uses mature transcript length and proteomics for promoter strength.

import { getProteinCount } from "../data/Data";
import { rnaPolKcat, baseMap, rnaPolKd, ATPconc, CTPconc, UTPconc, GTPconc, rnaPolK0, RnaPconc, riboKcat, trnaPolKcat, rrnaPolKcat, CytoVolume, riboKd, riboK0, ribosomeConc, ctRNAconc } from "../Globals";
import { aaMap, aaTRNAMap } from "./RNAMaps";

function checkBases(sequence: string, baseMap: Map<string, any>) {
  for (const base of new Set(sequence)) {
    if (!baseMap.has(base)) {
      console.error(sequence, baseMap);
      throw new Error("Unknown base(s) in RNA sequence: " + base);
    }
  }
}

export function transcriptRate(rnasequence: string, ptnMetID: string, jcvi2ID: string) {
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

  const CMono1 = baseMap.get(rnasequence[0]) || 0;
  const CMono2 = baseMap.get(rnasequence[1]) || 0;

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const NMono_A = baseCount.get("A") || 0;
  const NMono_U = baseCount.get("C") || 0;
  const NMono_C = baseCount.get("G") || 0;
  const NMono_G = baseCount.get("U") || 0;

  const NMonoSum = NMono_A * rnaPolKd / ATPconc + NMono_C * rnaPolKd / CTPconc + NMono_U * rnaPolKd / UTPconc + NMono_G * rnaPolKd / GTPconc;

  const k_transcription = kcat_mod / ((1 + rnaPolK0 / RnaPconc) * (rnaPolKd ** 2) / (CMono1 * CMono2) + NMonoSum + n_tot - 1);

  return k_transcription;
}

// Define transcription rate for ribosomal protein-coding genes.

export function riboTranscriptRate(rnasequence: string) {
  // Add trascription reaction

  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const kcat_mod = rnaPolKcat * 500 / 180;

  // Add total number of monomers to parameter dict

  const CMono1 = baseMap.get(rnasequence[0])!;
  const CMono2 = baseMap.get(rnasequence[1])!;

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const NMono_A = baseCount.get("A") || 0;
  const NMono_U = baseCount.get("C") || 0;
  const NMono_C = baseCount.get("G") || 0;
  const NMono_G = baseCount.get("U") || 0;

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

export function trnaTranscriptRate(rnasequence: string) {
  // Add trascription reaction

  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const NMono_A = baseCount.get("A") || 0;
  const NMono_U = baseCount.get("C") || 0;
  const NMono_C = baseCount.get("G") || 0;
  const NMono_G = baseCount.get("U") || 0;

  const NMonoDict = [NMono_A, NMono_C, NMono_G, NMono_U];

  const CMono1 = baseMap.get(rnasequence[0])!;
  const CMono2 = baseMap.get(rnasequence[1])!;

  const NMonoSum = NMono_A * rnaPolKd / ATPconc + NMono_C * rnaPolKd / CTPconc + NMono_U * rnaPolKd / UTPconc + NMono_G * rnaPolKd / GTPconc;

  const k_transcription = trnaPolKcat / ((1 + rnaPolK0 / RnaPconc) * (rnaPolKd ** 2) / (CMono1 * CMono2) + NMonoSum + n_tot - 1);

  return k_transcription;
}

// Define transcription rate for rRNA genes.

export function rrnaTranscriptRate(rnasequence: string) {
  // Add trascription reaction

  // Check that we know all bases used in the sequence
  checkBases(rnasequence, baseMap);

  // Count how many times each base is used
  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const n_tot = [...baseCount.values()].reduce((a, b) => a + b);

  const NMono_A = baseCount.get("A") || 0;
  const NMono_U = baseCount.get("C") || 0;
  const NMono_C = baseCount.get("G") || 0;
  const NMono_G = baseCount.get("U") || 0;

  const NMonoDict = [NMono_A, NMono_C, NMono_G, NMono_U];

  const CMono1 = baseMap.get(rnasequence[0])!;
  const CMono2 = baseMap.get(rnasequence[1])!;

  const NMonoSum = NMono_A * rnaPolKd / ATPconc + NMono_C * rnaPolKd / CTPconc + NMono_U * rnaPolKd / UTPconc + NMono_G * rnaPolKd / GTPconc;

  const k_transcription = rrnaPolKcat / ((1 + rnaPolK0 / RnaPconc) * (rnaPolKd ** 2) / (CMono1 * CMono2) + NMonoSum + n_tot - 1);

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

export function translateRate(rnasequence: string, aasequence: string) {
  // Add translation reaction

  // Check that we know all residues used in the sequence
  checkBases(aasequence, aaMap);

  // Count how many times each residue is used
  const aaCount = new Map<string, number>();
  for (const aa of new Set(aasequence)) {
    aaCount.set(aa, aasequence.split(aa).length - 1);
  }

  const NMono_A = aaCount.get("A") || 0;
  const NMono_R = aaCount.get("R") || 0;
  const NMono_N = aaCount.get("N") || 0;
  const NMono_D = aaCount.get("D") || 0;
  const NMono_C = aaCount.get("C") || 0;
  const NMono_E = aaCount.get("E") || 0;
  const NMono_Q = aaCount.get("Q") || 0;
  const NMono_H = aaCount.get("H") || 0;
  const NMono_I = aaCount.get("I") || 0;
  const NMono_L = aaCount.get("L") || 0;
  const NMono_K = aaCount.get("K") || 0;
  const NMono_M = aaCount.get("M") || 0;
  const NMono_P = aaCount.get("P") || 0;
  const NMono_S = aaCount.get("S") || 0;
  const NMono_T = aaCount.get("T") || 0;
  const NMono_W = aaCount.get("W") || 0;
  const NMono_Y = aaCount.get("Y") || 0;
  const NMono_G = aaCount.get("G") || 0;
  const NMono_F = aaCount.get("F") || 0;
  const NMono_V = aaCount.get("V") || 0;

  const NStop = aaCount.get("*") || 0;

  if (NStop > 1) {
    console.log("EXTRA STOP CODON: MISTAKE IN TRANSLATION");
  }

  const NMonoDict = [NMono_A, NMono_R, NMono_N, NMono_D, NMono_C, NMono_E, NMono_Q, NMono_H,
    NMono_I, NMono_L, NMono_K, NMono_M, NMono_P, NMono_S, NMono_T, NMono_W,
    NMono_Y, NMono_G, NMono_F, NMono_V];

  let NMonoSum = 0;

  for (const nmono of NMonoDict) {
    NMonoSum += nmono * riboKd / ctRNAconc;
  }

  const n_tot = [...aaCount.values()].reduce((a, b) => a + b, 0);

  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence)) {
    baseCount.set(base, rnasequence.split(base).length - 1);
  }

  const transcript_length = [...baseCount.values()].reduce((a, b) => a + b, 0);

  let ribo_num = Math.max(1, Math.round(transcript_length / 125 - 1));
  ribo_num = Math.min(15, ribo_num);

  let kcat_mod
  if (ribo_num > 1) {
    kcat_mod = ribo_num * 0.25 * riboKcat + 0.25 * riboKcat; // 503 ribosomes
  } else {
    kcat_mod = 0.45 * riboKcat;
  }

  const k_translation = kcat_mod / ((1 + riboK0 / ribosomeConc) * (riboKd ** 2) / (ctRNAconc ** 2) + NMonoSum + n_tot - 1);

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

  const ptnLen = [...aaCount.values()].reduce((a, b) => a + b, 0);

  const k_transloc = 50 / ptnLen;

  return k_transloc;
}