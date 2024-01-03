import { volToSA } from "./Globals";
import { ODE } from "./ODE";

function addLipidParams(model: ODE) {
  // ACPPAT: Raise forward kcat, Raise substrate km to literature values
  const kcatf_ACPPAT = 20.7; // 1/s kcat propanoyl phosphate - S. enterica EC
  const km_ap_ACPPAT = 0.0179; // mM - ap - S. enterica (BRENDA)

  model.addParameter("ACPPAT", "kcatF", kcatf_ACPPAT, 0, 0, "1/s", "kcatF_ACPPAT");
  model.addParameter("ACPPAT", "KmSub2", km_ap_ACPPAT, 0, 0, "mM", "km_Sub2_ACPPAT");

  // AGPAT: Fix kcat and km values to literature values for increased flux
  const kcatf_AGPAT = 144.4; // 1/s, C. moschata (SABIO-RK) - converted from Vmax, Slabas et al

  model.addParameter("AGPAT", "kcatF", kcatf_AGPAT, 0, 0, "1/s", "kcatF_AGPAT");
  model.addParameter("AGPAT", "kcatR", 0, 0, 0, "1/s", "kcatR_AGPAT"); // Irreversible

  // PGSA: Fix kcat backwards to zero, this is an irreversible reaction, otherwise CTP issue
  model.addParameter("PGSA", "kcatR", 0, 0, 0, "1/s", "kcatR_PGSA"); // Irreversible Reaction

  // Set PG Synthesis the be irreversible
  model.addParameter("PGPP", "kcatR", 0, 0, 0, "1/s", "kcatR_PGPP");

  // UDPGALM: galactomutase rate selection from literature, to obtain proper galactolipid synthesis rate
  const kcatf_UDPGALM = 27.0; // 1/s E. coli (BRENDA)
  const kcatr_UDPGALM = 1.5; // 1/s E. coli (BRENDA)
  const kmProd1_UDPGALM = 0.45; // mM Aspergillus (BRENDA)

  model.addParameter("UDPGALM", "KmProd1", kmProd1_UDPGALM, 0, 0, "mM", "km_Prod1_UDPGALM");
  model.addParameter("UDPGALM", "kcatR", kcatr_UDPGALM, 0, 0, "1/s", "kcatR_UDPGALM");
  model.addParameter("UDPGALM", "kcatF", kcatf_UDPGALM, 0, 0, "1/s", "kcatF_UDPGALM");

  // DAGGALT: Glycolipid synthesis rate selection from literature, set galactolipid synthesis to be irreversible
  const kcatr_DAGGALT = 0.0; // !/s Irreversible reaction

  model.addParameter("DAGGALT", "kcatR", kcatr_DAGGALT, 0, 0, "1/s", "kcatR_UDPGALM");
}

function FAFBA() {
  return '($Kcat_R)';
}

function COAabcR() {
  return '($KcatR)*($Enzyme)*($Subex/($Subex+$Km1))';
}

function ActTrans() {
  return '($Kcat_R*$Enz_R)';
}

function addTransportParams(model: ODE, pmap: Map<string, number>) {
  // Add the rate forms that match FBA fluxes and account for ATP dependence
  model.addRateForm("FAFBA", FAFBA());
  model.addRateForm("COAabcR", COAabcR());

  const Enz_COAabc = Math.min(pmap.get("M_PTN_JCVISYN3A_0641_c"), pmap.get("M_PTN_JCVISYN3A_0642_c"), pmap.get("M_PTN_JCVISYN3A_0643_c"), pmap.get("M_PTN_JCVISYN3A_0836_c")); // But this has to be all COAabc genes (See M. Melo)

  const Kcat_R_COAabc = 2.0; // 1/s, Santos et al (2018), CbrT params
  const Km_R_COAabc = 2.1e-6; // mM, Santos et al (2018), CbrT params

  const M_coa_e = 0.0013; // 1.3 uM (Tables 1&2 def med v3 Table 2)

  // Define the COAabc reaction
  model.addReaction("COAabc", "COAabcR", "Coenzyme A Transport");
  model.addProduct("COAabc", "Prod1", "M_coa_c");
  model.addParameter("COAabc", "Subex", M_coa_e);
  model.addSubstrate("COAabc", "Sub2", "M_atp_c"); // Update RL
  model.addProduct("COAabc", "Prod2", "M_adp_c"); // Update RL
  model.addParameter("COAabc", "Prod3", "M_pi_c"); // Update RL
  model.addParameter("COAabc", "KcatR", Kcat_R_COAabc);
  model.addParameter("COAabc", "Km1", Km_R_COAabc);
  model.addParameter("COAabc", "Enzyme", Enz_COAabc);

  // Add Sphingomyelin uptake
  model.addMetabolite("M_sm_c", "Sphingomyelin", 0);
  const Kcat_FBA_SM = 1.1e-3 * (volToSA / 2.0) * 1.25; //*volToSA#/VolToSA # mM/s, Uptake flux from lipidomics adjusted FBA

  model.addReaction("SMuptake", "FAFBA", "Sphingomyelin Uptake");
  model.addProduct("SMuptake", "Prod1", "M_sm_c");
  model.addParameter("SMuptake", "Kcat_R", Kcat_FBA_SM);

  // Add Phosphatidlycholine uptake to the model
  model.addMetabolite("M_pc_c", "Phosphatidylcholine", 0);
  const Kcat_FBA_PC = 1.42e-4 * (volToSA / 2.0) * 1.25; //*volToSA#/VolToSA # mM/s, Uptake flux from lipidomics adjusted FBA
  model.addReaction("PCuptake", "FAFBA", "Phosphatidylcholine Uptake");
  model.addProduct("PCuptake", "Prod1", "M_pc_c");
  model.addProduct("PCuptake", "Prod2", "M_adp_c");
  model.addProduct("PCuptake", "Prod3", "M_pi_c");
  model.addSubstrate("PCuptake", "Sub1", "M_atp_c");
  model.addParameter("PCuptake", "Kcat_R", Kcat_FBA_PC);

  // Add Triacylglycerol uptake
  model.addMetabolite("M_tag_c", "triacylglycerol", 0);
  const Kcat_FBA_TAG = 0.0003095 * (volToSA / 2.0) * 1.25; //*volToSA # mM/s, Uptake Flux Triacylglycerol from FBA

  model.addReaction("TAGt", "FAFBA", "TAG uptake");
  model.addProduct("TAGt", "Prod1", "M_tag_c");
  model.addParameter("TAGt", "Kcat_R", Kcat_FBA_TAG);
}

export function lipidPatch(model: ODE, pmap: Map<string, number>) {
  addLipidParams(model);
  addTransportParams(model, pmap);
}