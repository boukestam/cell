// Create Definitions to Get Protein Counts and Gene Sequences

import { ParsedGenbank } from "genbank-parser";
import { Datasheet } from "./Datasheet";
import { J2toAOE, defaultPtnCount, genomePtnLocDict, genomeRnaLocDict } from "./Globals";
import { loadExcel, loadGenbank, loadCSV, loadTSV, TSV, loadSBML } from "./Loader";
import { SBML } from "./SBML";

export const data: {
  annotatPD: Datasheet,
  genome2: ParsedGenbank,
  genome3A: ParsedGenbank,
  proteomPD: Datasheet,
  manGPRPD: Datasheet,
  PtnMetDF: Datasheet,
  memPtnMetDF: Datasheet,
  riboPtnMetDF: Datasheet,
  RNApolPtnMetDF: Datasheet,
  rrnaMetDF_1: Datasheet,
  rrnaMetDF_2: Datasheet,
  trnaMetDF: Datasheet,
  mRNA_counts_DF: Datasheet,
  ComDF_nuc: TSV,
  ComDF_cent: TSV,
  ComDF_lip: TSV,
  reconstPD: Datasheet,
  docSBML: SBML,
  nuclMetRxn: Datasheet,
  transport: TSV,
  modelParameters: Datasheet
} = {
  annotatPD: new Datasheet([]),
  genome2: undefined,
  genome3A: undefined,
  proteomPD: new Datasheet([]),
  manGPRPD: new Datasheet([]),
  PtnMetDF: new Datasheet([]),
  memPtnMetDF: new Datasheet([]),
  riboPtnMetDF: new Datasheet([]),
  RNApolPtnMetDF: new Datasheet([]),
  rrnaMetDF_1: new Datasheet([]),
  rrnaMetDF_2: new Datasheet([]),
  trnaMetDF: new Datasheet([]),
  mRNA_counts_DF: new Datasheet([]),
  ComDF_nuc: new Map(),
  ComDF_cent: new Map(),
  ComDF_lip: new Map(),
  reconstPD: new Datasheet([]),
  docSBML: undefined,
  nuclMetRxn: new Datasheet([]),
  transport: new Map(),
  modelParameters: new Datasheet([])
};

export async function loadData() {
  const [
    annotatPD, genome2, genome3A, proteomPD, manGPRPD,
    PtnMetDF, memPtnMetDF, riboPtnMetDF, RNApolPtnMetDF, rrnaMetDF_1, rrnaMetDF_2, trnaMetDF, mRNA_counts_DF,
    ComDF_nuc, ComDF_cent, ComDF_lip,
    reconstPD, docSBML, nuclMetRxn, transport, modelParameters
  ] = await Promise.all([
    // The annotation matches MMSYN1* IDs with JCVISYN3* IDs (or "locus tags").
    loadExcel("model_data/FBA/Syn3A_annotation_compilation.xlsx", "Syn3A_annotation_compilation_co", true),

    // The genome data matches "locus tags" with AOE* protein IDs.
    // It provides both the gene sequence, needed for transcription reactions in the ODE model,
    // and the protein sequence, needed for translation reactions in the model.
    // This is the NCBI Gene Bank-formated file (https://www.ncbi.nlm.nih.gov/nuccore/CP014992.1).
    loadGenbank("model_data/syn2.gb"),

    // This is the NCBI Gene Bank-formated file (https://www.ncbi.nlm.nih.gov/nuccore/CP016816.2).
    loadGenbank("model_data/syn3A.gb"),

    // The proteomics matches AOE IDs with quantitative proteomics data.
    loadExcel("model_data/proteomics.xlsx", "Proteomics", true, 1),

    // The manual GPR conversion matches proteins with known abundances (AOE IDs), with reactions that did not 
    // have associated genes in the Syn3A reconstruction.
    loadCSV("model_data/manual_GPR_conversion.csv", false, ["MM", "AOE"]),

    loadCSV("model_data/protein_metabolites_frac.csv", true),
    loadCSV("model_data/membrane_protein_metabolites.csv", true),
    loadCSV("model_data/ribo_protein_metabolites.csv", true),
    loadCSV("model_data/RNApol_proteins.csv", true),
    loadCSV("model_data/rrna_metabolites_1.csv", true),
    loadCSV("model_data/rrna_metabolites_2.csv", true),
    loadCSV("model_data/trna_metabolites_synthase.csv", true),
    loadCSV("model_data/mRNA_counts.csv", true),

    loadTSV("model_data/Nucleotide_Kinetic_Parameters.tsv"),
    loadTSV("model_data/Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv"),
    loadTSV("model_data/lipid_NoH2O_balanced_model.tsv"),

    loadExcel("model_data/reconstruction.xlsx", "Reactions", true),
    loadSBML("model_data/iMB155_NoH2O.xml"),
    loadCSV("model_data/nucleo_rxns_list.txt", false),
    loadTSV("model_data/transport_NoH2O_Zane-TB-DB.tsv"),
    loadCSV("model_data/GlobalParameters_Zane-TB-DB.csv", true)
  ]);

  data.annotatPD = annotatPD;
  data.genome2 = genome2;
  data.genome3A = genome3A;
  data.proteomPD = proteomPD;
  data.manGPRPD = manGPRPD;

  data.PtnMetDF = PtnMetDF;
  data.memPtnMetDF = memPtnMetDF;
  data.riboPtnMetDF = riboPtnMetDF;
  data.RNApolPtnMetDF = RNApolPtnMetDF;
  data.rrnaMetDF_1 = rrnaMetDF_1;
  data.rrnaMetDF_2 = rrnaMetDF_2;
  data.trnaMetDF = trnaMetDF;
  data.mRNA_counts_DF = mRNA_counts_DF;

  data.ComDF_nuc = ComDF_nuc;
  data.ComDF_cent = ComDF_cent;
  data.ComDF_lip = ComDF_lip;

  data.reconstPD = reconstPD;
  data.docSBML = docSBML;
  data.nuclMetRxn = nuclMetRxn;
  data.transport = transport;
  data.modelParameters = modelParameters;

  return data;
}

// Create list of proteins with no proteomics data
const ptnNoQuant = new Set<string>();

export function getProteinCount(newMetID: string, jcvi2ID: string): { proteinCount: number, proteinName: string } {
  // Check if protein quantification exists

  try {
    const aoeID = jcvi2ID.startsWith("JCVIman_") ?
      data.manGPRPD.findRow("MM", jcvi2ID.replace("JCVIman_", "")).AOE :
      J2toAOE.get(jcvi2ID);

    const proteinCount = Math.max(
      defaultPtnCount,
      Math.round(data.proteomPD.findRow(0, aoeID)[21])
    );

    const proteinName = data.proteomPD.findRow(0, aoeID)[1].replace(" [synthetic bacterium JCVI-Syn3.0]", "");

    return { proteinCount, proteinName };
  } catch {
    ptnNoQuant.add(newMetID);
    return { proteinCount: defaultPtnCount, proteinName: newMetID };
  }
}

export function getSequences(jcvi3AID: string) {
  // returns genomic and protein sequences
  try {
    const rnasequence = genomePtnLocDict.get(jcvi3AID).extract(data.genome3A.sequence).transcribe();

    // Using translation table 11 from NCBI: "Bacterial, Archaeal and Plant Plastid Code"
    // https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4
    const aasequence = genomePtnLocDict.get(jcvi3AID).extract(data.genome3A.sequence).transcribe().translate(4);

    return { rnasequence, aasequence };
  } catch {
    return { rnasequence: undefined, aasequence: undefined };
  }
}

export function getRNASequences(jcvi3AID: string) {
  // returns genomic and protein sequences
  try {
    return genomeRnaLocDict.get(jcvi3AID).extract(data.genome3A.sequence).transcribe();
  } catch {
    return undefined;
  }
}