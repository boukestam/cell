// Create Definitions to Add RNA, Proteins, and Reactions to the Simulation

import { CME } from "./CME";
import { data, getProteinCount, getRNASequences, getSequences } from "./Data";
import { Gene } from "./Gene";
import { ModelSpecies, genomePtnLocDict, genomeRnaLocDict } from "./Globals";
import { poisson } from "./Poisson";
import { aaCostMap, aaTRNAMap, ctRNAcostMap } from "./RNAMaps";
import { transcriptRate, degradationRate, translateRate, riboTranscriptRate, trnaTranscriptRate } from "./Rates";

export const genesInModel = new Set<string>();

let unkIter = 1;

function addRNACost(
  sim: CME,
  geneMetID: string, rnaMetID: string, ptnMetID: string, jcvi2ID: string,
  rnasequence: Gene, aasequence: Gene
) {
  const TrscProd = [geneMetID, rnaMetID];
  const mRNAdegProd = [];

  for (let i = 0; i < rnasequence.length; i++) {
    TrscProd.push('ATP_trsc');
    mRNAdegProd.push('ATP_mRNAdeg');
  }

  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence.sequence)) {
    baseCount.set(base, rnasequence.sequence.split(base).length - 1);
  }

  // Add total number of monomers to parameter dict

  const N_A = baseCount.get("A")!;
  const N_U = baseCount.get("U")!;
  const N_C = baseCount.get("C")!;
  const N_G = baseCount.get("G")!;

  for (let i = 0; i < N_A; i++) {
    TrscProd.push('ATP_mRNA');
    mRNAdegProd.push('AMP_mRNAdeg');
  }

  for (let i = 0; i < N_U; i++) {
    TrscProd.push('UTP_mRNA');
    mRNAdegProd.push('UMP_mRNAdeg');
  }

  for (let i = 0; i < N_G; i++) {
    TrscProd.push('GTP_mRNA');
    mRNAdegProd.push('GMP_mRNAdeg');
  }

  for (let i = 0; i < N_C; i++) {
    TrscProd.push('CTP_mRNA');
    mRNAdegProd.push('CMP_mRNAdeg');
  }

  const TranslatProd = [rnaMetID, ptnMetID];

  for (let i = 0; i < aasequence.length; i++) {
    TranslatProd.push('ATP_translat');
    TranslatProd.push('ATP_translat');
  }

  const aaCount = new Map<string, number>();
  for (const aa of new Set(aasequence.sequence)) {
    aaCount.set(aa, aasequence.sequence.split(aa).length - 1);
  }

  const NMono_A = aaCount.get("A")!;
  const NMono_R = aaCount.get("R")!;
  const NMono_N = aaCount.get("N")!;
  const NMono_D = aaCount.get("D")!;
  const NMono_C = aaCount.get("C")!;
  const NMono_E = aaCount.get("E")!;
  const NMono_Q = aaCount.get("Q")!;
  const NMono_H = aaCount.get("H")!;
  const NMono_I = aaCount.get("I")!;
  const NMono_L = aaCount.get("L")!;
  const NMono_K = aaCount.get("K")!;
  const NMono_M = Math.max(0, aaCount.get("M")! - 1);
  const NMono_P = aaCount.get("P")!;
  const NMono_S = aaCount.get("S")!;
  const NMono_T = aaCount.get("T")!;
  const NMono_W = aaCount.get("W")!;
  const NMono_Y = aaCount.get("Y")!;
  const NMono_G = aaCount.get("G")!;
  const NMono_F = aaCount.get("F")!;
  const NMono_V = aaCount.get("V")!;
  const NMono_FM = 1;

  const NStop = aaCount.get("*")!;

  if (NStop > 1) {
    console.log("EXTRA STOP CODON: MISTAKE IN TRANSLATION");
  }

  const AaUsed: [string, number][] = [['A', NMono_A], ['R', NMono_R], ['N', NMono_N], ['D', NMono_D], ['C', NMono_C],
  ['E', NMono_E], ['Q', NMono_Q], ['H', NMono_H], ['I', NMono_I], ['L', NMono_L],
  ['K', NMono_K], ['M', NMono_M], ['P', NMono_P], ['S', NMono_S], ['T', NMono_T],
  ['W', NMono_W], ['Y', NMono_Y], ['G', NMono_G], ['F', NMono_F], ['V', NMono_V],
  ['FM', NMono_FM]];

  for (const aaCost of AaUsed) {
    const aa_ID = aaCost[0];
    const numberUsed = aaCost[1];

    const aaCostID = aaCostMap.get(aa_ID)!;

    for (let i = 0; i < numberUsed; i++) {
      TranslatProd.push(aaCostID);
    }
  }

  sim.addReaction([geneMetID], TrscProd, transcriptRate(rnaMetID, ptnMetID, rnasequence.sequence, jcvi2ID));
  sim.addReaction([rnaMetID], mRNAdegProd, degradationRate(rnaMetID, rnasequence.sequence));
  sim.addReaction([rnaMetID], TranslatProd, translateRate(rnaMetID, ptnMetID, rnasequence.sequence, aasequence.sequence, sim.species));
}

function getJCVI2ID(mmcode: string) {
  // Checks if a translation to JCVISYN2* code is available
  return data.annotatPD.findRow(5, mmcode)?.[13]?.trim() || (
    data.manGPRPD.getColumn("MM").includes(mmcode) ?
      ("JCVIman_" + mmcode) :
      ("JCVIunk_" + mmcode + "_" + (++unkIter))
  );
}

function addProtein(sim: CME, jcvi3AID: string) {
  const locusNum = jcvi3AID.split('_')[1];
  const mmcode = 'MMSYN1_' + locusNum;

  const jcvi2ID = getJCVI2ID(mmcode);

  genesInModel.add(jcvi3AID);

  const ptnMetID = 'M_PTN_' + jcvi3AID + '_c';
  ModelSpecies.add(ptnMetID);

  const { proteinCount } = getProteinCount(ptnMetID, jcvi2ID);

  const geneMetID = jcvi3AID + '_gene';
  ModelSpecies.add(geneMetID);

  // Get nucleotide and amino acid sequences, if available
  const { rnasequence, aasequence } = getSequences(jcvi3AID);

  if (!rnasequence || !aasequence) return;

  const rnaMetID = "M_RNA_" + jcvi3AID + "_c";
  ModelSpecies.add(rnaMetID);

  sim.defineSpecies([geneMetID, rnaMetID, ptnMetID]);

  sim.addParticles(geneMetID, 1);
  sim.addParticles(ptnMetID, Math.max(1, proteinCount));
  sim.addParticles(rnaMetID, 1);

  addRNACost(sim, geneMetID, rnaMetID, ptnMetID, jcvi2ID, rnasequence, aasequence);
}

function addMetabolite(sim: CME, newMetID: string, proteinCount: number, transcribe: boolean, jcvi3AID: string) {
  if (!transcribe) {
    sim.defineSpecies([newMetID]);
    sim.addParticles(newMetID, proteinCount);
    return;
  }

  const geneMetID = jcvi3AID + '_gene';
  ModelSpecies.add(geneMetID);

  // Get nucleotide and amino acid sequences, if available
  const { rnasequence, aasequence } = getSequences(jcvi3AID);

  if (!rnasequence || !aasequence) return;

  const rnaMetID = "M_RNA_" + jcvi3AID + "_c";
  ModelSpecies.add(rnaMetID);

  sim.defineSpecies([geneMetID, rnaMetID, newMetID]);

  const row = data.mRNA_counts_DF.findRow("LocusTag", jcvi3AID);
  let avg_mRNA_cnt = row[2];
  if (avg_mRNA_cnt === 0) avg_mRNA_cnt = 0.001;
  const init_mRNA_count = poisson(avg_mRNA_cnt);

  sim.addParticles(geneMetID, 1);
  sim.addParticles(newMetID, proteinCount);
  sim.addParticles(rnaMetID, init_mRNA_count);

  return { geneMetID, rnaMetID, rnasequence, aasequence };
}

function addNamedProtein(
  sim: CME,
  newMetID: string, mmcode: string, transcribe: boolean, proteinFrac: number
) {
  const locusNum = mmcode.split('_')[1];

  const jcvi3AID = 'JCVISYN3A_' + locusNum;
  genesInModel.add(jcvi3AID);

  const jcvi2ID = getJCVI2ID(mmcode);

  let { proteinCount, proteinName } = getProteinCount(newMetID, jcvi2ID);
  proteinCount = proteinFrac * proteinCount;

  ModelSpecies.add(newMetID);

  const metaboliteResult = addMetabolite(
    sim, newMetID, proteinCount, transcribe, jcvi3AID
  );
  if (!metaboliteResult) return;

  const { geneMetID, rnaMetID, rnasequence, aasequence } = metaboliteResult;
  addRNACost(sim, geneMetID, rnaMetID, newMetID, jcvi2ID, rnasequence, aasequence);
}

function addMembraneProtein(
  sim: CME,
  newMetID: string, mmcode: string, transcribe: boolean, proteinFrac: number
) {
  const locusNum = mmcode.split('_')[1];

  const jcvi3AID = 'JCVISYN3A_' + locusNum;
  genesInModel.add(jcvi3AID);

  const jcvi2ID = getJCVI2ID(mmcode);

  let { proteinCount, proteinName } = getProteinCount(newMetID, jcvi2ID);
  proteinCount = proteinFrac * proteinCount;

  ModelSpecies.add(newMetID);

  const metaboliteResult = addMetabolite(
    sim, newMetID, proteinCount, transcribe, jcvi3AID
  );
  if (!metaboliteResult) return;

  const { geneMetID, rnaMetID, rnasequence, aasequence } = metaboliteResult;
  addRNACost(sim, geneMetID, rnaMetID, newMetID, jcvi2ID, rnasequence, aasequence);
}

const RiboPtnNames = ['Ribosomal Protein'];
const RiboPtnLens: any[] = ['Gene Length (nt)'];
const RiboPtnTrscRates: any[] = ['Transcription Rate'];
const RiboPtnTranslatRates: any[] = ['Translation Rate'];

function addRiboProtein(
  sim: CME,
  newMetID: string, mmcode: string
) {
  const locusNum = mmcode.split('_')[1];

  const jcvi3AID = 'JCVISYN3A_' + locusNum;
  genesInModel.add(jcvi3AID);

  const jcvi2ID = getJCVI2ID(mmcode);

  const { proteinCount, proteinName } = getProteinCount(newMetID, jcvi2ID);

  ModelSpecies.add(newMetID);

  const metaboliteResult = addMetabolite(
    sim, newMetID, proteinCount, true, jcvi3AID
  );
  if (!metaboliteResult) return;

  const { geneMetID, rnaMetID, rnasequence, aasequence } = metaboliteResult;

  RiboPtnNames.push(newMetID);

  const trsc_rate = riboTranscriptRate(rnaMetID, newMetID, rnasequence.sequence, jcvi2ID);
  RiboPtnTrscRates.push(trsc_rate);

  const translat_rate = translateRate(rnaMetID, newMetID, rnasequence.sequence, aasequence.sequence, sim.species);
  RiboPtnTranslatRates.push(translat_rate);

  const rna_length = rnasequence.length;
  RiboPtnLens.push(rna_length)

  addRNACost(sim, geneMetID, rnaMetID, newMetID, jcvi2ID, rnasequence, aasequence);
}

const tRNAadded = new Set<string>();

function addTRNA(sim: CME, unchargedMetID: string, chargedMetID: string, mmcode: string, aminoAcid: string, synthase: string) {
  const locusNum = mmcode.split('_')[1];

  const jcvi3AID = 'JCVISYN3A_' + locusNum;

  const jcvi2ID = getJCVI2ID(mmcode);

  const geneMetID = jcvi3AID + '_gene';

  ModelSpecies.add(geneMetID);

  const rnasequence = getRNASequences(jcvi3AID);

  const gene = [geneMetID];
  sim.defineSpecies(gene);

  sim.addParticles(geneMetID, 1);

  if (!tRNAadded.has(unchargedMetID)) {
    ModelSpecies.add(unchargedMetID);
    const tRNA = [unchargedMetID];
    sim.defineSpecies(tRNA);
    tRNAadded.add(unchargedMetID);
    ModelSpecies.add(chargedMetID);
    const ctRNA = [chargedMetID];
    sim.defineSpecies(ctRNA);
    tRNAadded.add(chargedMetID);

    if (!synthase.includes('GLN')) {
      const synthase_atp = synthase + '_ATP';
      const synth_atp = [synthase_atp];
      sim.defineSpecies(synth_atp);
      sim.addParticles(synthase_atp, 1);

      const synthase_atp_aa = synthase_atp + '_AA';
      const synth_atp_aa = [synthase_atp_aa];
      sim.defineSpecies(synth_atp_aa);
      sim.addParticles(synthase_atp_aa, 1);

      const synthase_atp_aa_trna = synthase_atp_aa + '_tRNA';
      const synth_atp_aa_trna = [synthase_atp_aa_trna];
      sim.defineSpecies(synth_atp_aa_trna);
      sim.addParticles(synthase_atp_aa_trna, 1);

      const synthase_ptn = 'M_PTN_' + synthase + '_c';

      sim.addReaction(['M_atp_c', synthase_ptn], [synthase_atp], 0.01);
      sim.addReaction([aminoAcid, synthase_atp], [synthase_atp_aa], 0.01);
      sim.addReaction([synthase_atp_aa, unchargedMetID], [synthase_atp_aa_trna], 0.1);
      sim.addReaction([synthase_atp_aa_trna], [chargedMetID, 'M_amp_c', 'M_ppi_c', synthase_ptn], 30);

      const costID = ctRNAcostMap.get(chargedMetID)!;

      const cost_paid = costID + '_paid';

      sim.addReaction([chargedMetID, costID], [unchargedMetID, cost_paid], 100);
    } else {
      const synthase = 'JCVISYN3A_0126';
      const synthase_atp = synthase + '_ATP';
      const synthase_atp_aa = synthase_atp + '_AA';
      const synthase_ptn = 'M_PTN_' + synthase + '_c';

      const synthase_atp_aa_trna = synthase_atp_aa + '_tRNAgln';
      sim.defineSpecies([synthase_atp_aa, synthase_atp_aa_trna]);
      sim.addParticles(synthase_atp_aa_trna, 1);

      sim.addReaction([synthase_atp_aa, 'M_trnagln_c'], [synthase_atp_aa_trna], 0.1);
      sim.addReaction([synthase_atp_aa_trna], ['M_glutrnagln_c', 'M_amp_c', 'M_ppi_c', synthase_ptn], 30);
      const glugln = 'M_glutrnagln_c';
      const glugln_enz = 'glutrnagln_enz';
      const glugln_enz_atp = 'glutrnagln_enz_atp';
      const glugln_enz_atp_gln = 'glutrnagln_enz_atp_gln';
      const glugln_enz_atp_asn = 'glutrnagln_enz_atp_asn';
      const gge = [glugln_enz];
      sim.defineSpecies(gge);
      sim.addParticles(glugln_enz, 1);
      const ggea = [glugln_enz_atp];
      sim.defineSpecies(ggea);
      sim.addParticles(glugln_enz_atp, 1);
      const ggeagln = [glugln_enz_atp_gln];
      sim.defineSpecies(ggeagln);
      sim.addParticles(glugln_enz_atp_gln, 1);
      const ggeaasn = [glugln_enz_atp_asn];
      sim.defineSpecies(ggeaasn);
      sim.addParticles(glugln_enz_atp_asn, 1);

      sim.addReaction([glugln, 'M_PTN_JCVISYN3A_0689_c'], [glugln_enz], 0.01);
      sim.addReaction([glugln_enz, 'M_atp_c'], [glugln_enz_atp], 0.01);
      sim.addReaction([glugln_enz_atp, 'M_gln__L_c'], [glugln_enz_atp_gln], 0.01);
      sim.addReaction([glugln_enz_atp_gln], ['M_glntrna_c', 'M_glu__L_c', 'M_adp_c', 'M_pi_c', 'M_PTN_JCVISYN3A_0689_c'], 30);
      sim.addReaction([glugln_enz_atp, 'M_asn__L_c'], [glugln_enz_atp_asn], 0.001);
      sim.addReaction([glugln_enz_atp_asn], ['M_glntrna_c', 'M_asp__L_c', 'M_adp_c', 'M_pi_c', 'M_PTN_JCVISYN3A_0689_c'], 30);

      const costID = ctRNAcostMap.get('M_glntrna_c')!;

      const cost_paid = costID + '_paid';

      sim.addReaction(['M_glntrna_c', costID], ['M_trnagln_c', cost_paid], 100);
    }
  }

  sim.addParticles(unchargedMetID, Math.round(5800 / 29 * 0.2));
  sim.addParticles(chargedMetID, Math.round(5800 / 29 * 0.8));
  const TrscProd = [geneMetID, unchargedMetID];

  for (let i = 0; i < rnasequence.length; i++) {
    TrscProd.push('ATP_trsc');
  }

  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence.sequence)) {
    baseCount.set(base, rnasequence.sequence.split(base).length - 1);
  }

  // Add total number of monomers to parameter dict

  const N_A = baseCount.get("A")!;
  const N_U = baseCount.get("U")!;
  const N_C = baseCount.get("C")!;
  const N_G = baseCount.get("G")!;

  for (let i = 0; i < N_A; i++) {
    TrscProd.push('ATP_tRNA');
  }

  for (let i = 0; i < N_U; i++) {
    TrscProd.push('UTP_tRNA');
  }

  for (let i = 0; i < N_G; i++) {
    TrscProd.push('GTP_tRNA');
  }

  for (let i = 0; i < N_C; i++) {
    TrscProd.push('CTP_tRNA');
  }

  sim.addReaction([geneMetID], TrscProd, trnaTranscriptRate(unchargedMetID, rnasequence.sequence));
}

function addRRNA(sim: CME) {
  const rrnaMetDF_1 = data.rrnaMetDF_1;
  const rrnaMetDF_2 = data.rrnaMetDF_2;

  const rrna_gene_locs_1: number[] = [];
  const rrna_gene_locs_2: number[] = [];

  const rRNA_species: string[] = [];

  for (const row of rrnaMetDF_1.rows) {
    const newMetID = row[0];
    const jcvi3AID = row[1];

    const geneLocation = genomeRnaLocDict.get(jcvi3AID)!;
    rrna_gene_locs_1.push(geneLocation.start);
    rrna_gene_locs_1.push(geneLocation.end);

    rRNA_species.push(newMetID);
  }

  const rrna_pos_1 = new Gene(data.genome3A.sequence).slice(Math.min(...rrna_gene_locs_1), Math.max(...rrna_gene_locs_1));

  const rrna_operon_1 = rrna_pos_1.reverseComplement();

  for (const row of rrnaMetDF_2.rows) {
    const newMetID = row[0];
    const jcvi3AID = row[1];

    const geneLocation = genomeRnaLocDict.get(jcvi3AID)!;
    rrna_gene_locs_2.push(geneLocation.start);
    rrna_gene_locs_2.push(geneLocation.end);
  }

  const rrna_pos_2 = new Gene(data.genome3A.sequence).slice(Math.min(...rrna_gene_locs_2), Math.max(...rrna_gene_locs_2));

  const rrna_operon_2 = rrna_pos_2.reverseComplement();

  sim.defineSpecies(rRNA_species);

  for (const rRNA of rRNA_species) {
    sim.addParticles(rRNA, 1);
    ModelSpecies.add(rRNA);
  }

  for (const row of rrnaMetDF_1.rows) {
    const newMetID = row[0];
    const jcvi3AID = row[1];

    const geneMetID = jcvi3AID + '_gene';

    ModelSpecies.add(geneMetID);

    const gene = [geneMetID];
    sim.defineSpecies(gene);

    sim.addParticles(geneMetID, 1);
  }

  const rnasequence = rrna_operon_1.transcribe();

  const TrscProd = ['JCVISYN3A_0067_gene', 'M_rRNA_5S_c', 'M_rRNA_23S_c', 'M_rRNA_16S_c'];

  for (let i = 0; i < rnasequence.length; i++) {
    TrscProd.push('ATP_trsc');
  }

  const baseCount = new Map<string, number>();
  for (const base of new Set(rnasequence.sequence)) {
    baseCount.set(base, rnasequence.sequence.split(base).length - 1);
  }

  // Add total number of monomers to parameter dict

  const N_A = baseCount.get("A")!;
  const N_U = baseCount.get("U")!;
  const N_C = baseCount.get("C")!;
  const N_G = baseCount.get("G")!;

  for (let i = 0; i < N_A; i++) {
    TrscProd.push('ATP_rRNA');
  }

  for (let i = 0; i < N_U; i++) {
    TrscProd.push('UTP_rRNA');
  }

  for (let i = 0; i < N_G; i++) {
    TrscProd.push('GTP_rRNA');
  }

  for (let i = 0; i < N_C; i++) {
    TrscProd.push('CTP_rRNA');
  }

  sim.addReaction(['JCVISYN3A_0067_gene'], TrscProd, transcriptRate('M_rRNA_5S_c', 'M_rRNA_16S_c', rnasequence.sequence, 'JCVIunk_0067'));

  for (const row of rrnaMetDF_2.rows) {
    const newMetID = row[0];
    const jcvi3AID = row[1];

    const geneMetID = jcvi3AID + '_gene';

    ModelSpecies.add(geneMetID);

    const gene = [geneMetID];
    sim.defineSpecies(gene);

    sim.addParticles(geneMetID, 1);
  }

  const rnasequence2 = rrna_operon_2.transcribe();

  const TrscProd2 = ['JCVISYN3A_0532_gene', 'M_rRNA_5S_c', 'M_rRNA_23S_c', 'M_rRNA_16S_c'];

  for (let i = 0; i < rnasequence2.length; i++) {
    TrscProd2.push('ATP_trsc');
  }

  const baseCount2 = new Map<string, number>();
  for (const base of new Set(rnasequence2.sequence)) {
    baseCount2.set(base, rnasequence2.sequence.split(base).length - 1);
  }

  // Add total number of monomers to parameter dict

  const N_A2 = baseCount2.get("A")!;
  const N_U2 = baseCount2.get("U")!;
  const N_C2 = baseCount2.get("C")!;
  const N_G2 = baseCount2.get("G")!;

  for (let i = 0; i < N_A2; i++) {
    TrscProd2.push('ATP_rRNA');
  }

  for (let i = 0; i < N_U2; i++) {
    TrscProd2.push('UTP_rRNA');
  }

  for (let i = 0; i < N_G2; i++) {
    TrscProd2.push('GTP_rRNA');
  }

  for (let i = 0; i < N_C2; i++) {
    TrscProd2.push('CTP_rRNA');
  }

  sim.addReaction(['JCVISYN3A_0532_gene'], TrscProd2, transcriptRate('M_rRNA_5S_c', 'M_rRNA_16S_c', rnasequence2.sequence, 'JCVIunk_0532'));
}

export function populate(sim: CME) {
  const Nuc_counters = ['ATP_trsc', 'ATP_translat', 'ATP_mRNAdeg', 'ATP_ptndeg', 'ATP_DNArep', 'ATP_transloc',
    'ATP_mRNA', 'UTP_mRNA', 'CTP_mRNA', 'GTP_mRNA',
    'AMP_mRNAdeg', 'UMP_mRNAdeg', 'CMP_mRNAdeg', 'GMP_mRNAdeg',
    'ATP_tRNA', 'UTP_tRNA', 'CTP_tRNA', 'GTP_tRNA',
    'ATP_rRNA', 'UTP_rRNA', 'CTP_rRNA', 'GTP_rRNA',
    'dATP_DNArep', 'dTTP_DNArep', 'dCTP_DNArep', 'dGTP_DNArep'];

  sim.defineSpecies(Nuc_counters);

  const AA_counters = ["ALA_cost", "ARG_cost", "ASN_cost", "ASP_cost", "CYS_cost", "GLU_cost", "GLN_cost", "GLY_cost",
    "HIS_cost", "ILE_cost", "LEU_cost", "LYS_cost", "MET_cost", "PHE_cost", "PRO_cost", "SER_cost",
    "THR_cost", "TRP_cost", "TYR_cost", "VAL_cost", "FMET_cost"];

  sim.defineSpecies(AA_counters);

  const AA_paid: string[] = [];

  for (const cost of AA_counters) {
    const cost_paid = cost + '_paid';
    AA_paid.push(cost_paid);
  }

  sim.defineSpecies(AA_paid);

  for (const row of data.PtnMetDF.getRows()) {
    addNamedProtein(sim, row["species"], row["gene"], row["transcribe"] === '1', row["proteomics_fraction"]);
  }

  for (const row of data.memPtnMetDF.getRows()) {
    addMembraneProtein(sim, row["species"], row["gene"], row["transcribe"] === '1', row["proteomics_fraction"]);
  }

  for (const row of data.riboPtnMetDF.getRows()) {
    addRiboProtein(sim, row["species"], row["gene"]);
  }

  for (const gene of genomePtnLocDict.keys()) {
    if (!genesInModel.has(gene)) {
      addProtein(sim, gene);
    }
  }

  for (const row of data.trnaMetDF.getRows()) {
    addTRNA(sim, row["uncharged"], row["charged"], row["gene"], row["AA"], row["synthase"]);
  }

  addRRNA(sim);
}