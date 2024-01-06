import { data } from "../data/Data";
import { Datasheet } from "../data/Datasheet";
import { ODE, ODEParameterValue } from "./ODE";
import { enzymatic, partTomM } from "./Reactions";
import { SBMLModel } from "../data/SBML";

// Default metabolite concentration (mM)
const defaultMetConcentration = 0.1;

// Default Protein concentration (mM): set to 1 micro-Mol
const defaultPtnConcentration = defaultMetConcentration / 100;

function getKcat(rxnIDs: string[], QntDF: Datasheet, ParDF: Datasheet, CentQntDF: Datasheet, RxnDF: Datasheet) {
  let selectArray = new Array(QntDF.rows.length).fill(false);

  // Turns out it isn't very easy to select rows based on a list of substrings that may appear in
  //  values of a column. We iterate over reaction IDs and create a logical index of rows to keep.
  for (let rxnIDR of rxnIDs) {
    selectArray = selectArray.map((value, index) => value || QntDF.getRow(index).Quantity.includes(rxnIDR));
  }

  QntDF = QntDF.filter((row, index) => selectArray[index]);

  // Read in the parameters in the parameter dataframe for reactions in Central metabolism.
  // This gives us the forward and reverse kcats to use later in kinetic rates.

  selectArray = new Array(ParDF.rows.length).fill(false);

  for (let rxnIDR of rxnIDs) {
    selectArray = selectArray.map((value, index) => value || ParDF.getRow(index)["Reaction:SBML:reaction:id"].includes(rxnIDR));
  }

  ParDF = ParDF.filter((row, index) => selectArray[index]);

  // Read through the parameter dataframe and pull out forward and reverse kcat values.
  // Then create a dataframe of just the kcat values.

  const kcats = [];

  for (let rxnID of rxnIDs) {

    const kcatFID = "kcatF_" + rxnID;

    const rows = ParDF.filter(row => row["Reaction:SBML:reaction:id"] === rxnID);

    const kcatF = rows.map(row => row.Mode)[3];
    const kcatFunits = rows.map(row => row.Unit)[3];

    const kF = [kcatFID, kcatF, kcatFunits];

    const kcatRID = "kcatR_" + rxnID;

    const kcatR = rows.map(row => row.Mode)[4];
    const kcatRunits = rows.map(row => row.Unit)[4];

    const kR = [kcatRID, kcatR, kcatRunits];

    kcats.push(kF);
    kcats.push(kR);
  }

  const KcatDF = new Datasheet(kcats, ["Parameter", "Value", "Units"]);

  const Kms = [];

  for (let rxnID of rxnIDs) {
    let reaction: string = RxnDF.findRow("Reaction", rxnID).ReactionFormula;
    reaction = reaction.replace(/ /g, '');

    let substrates = reaction.split("<=>")[0];

    if (substrates.includes("2.0")) {
      substrates = substrates.replace(/2\.0/g, '');
    }
    if (substrates.includes("3.0")) {
      substrates = substrates.replace(/3\.0/g, '');
    }
    if (substrates.includes("4.0")) {
      substrates = substrates.replace(/4\.0/g, '');
    }

    if (substrates.includes("+")) {
      const plusNum = substrates.split('+').length - 1;

      for (let i = 0; i < plusNum + 1; i++) {

        const substrate = substrates.split('+')[i];

        const kmID = 'kmc_' + rxnID + '_' + substrate + '_' + rxnID;

        const km_value = CentQntDF.findRow("Quantity", kmID).Value;

        if (km_value) Kms.push([kmID, km_value, "mM"]);
      }
    } else {
      const kmID = 'kmc_' + rxnID + '_' + substrates + '_' + rxnID;

      const km_value = CentQntDF.findRow("Quantity", kmID).Value;

      Kms.push([kmID, km_value, "mM"]);
    }

    let products = reaction.split("<=>")[1];

    if (products.includes("2.0")) {
      products = products.replace(/2\.0/g, '');
    }
    if (products.includes("3.0")) {
      products = products.replace(/3\.0/g, '');
    }
    if (products.includes("4.0")) {
      products = products.replace(/4\.0/g, '');
    }

    if (products.includes("+")) {
      const plusNum = products.split('+').length - 1;

      for (let i = 0; i < plusNum + 1; i++) {

        const product = products.split('+')[i];

        const kmID = 'kmc_' + rxnID + '_' + product + '_' + rxnID;

        const km_value = CentQntDF.findRow("Quantity", kmID).Value;

        Kms.push([kmID, km_value, "mM"]);
      }
    } else {
      const kmID = 'kmc_' + rxnID + '_' + products + '_' + rxnID;

      const km_value = CentQntDF.findRow("Quantity", kmID).Value;

      Kms.push([kmID, km_value, "mM"]);
    }
  }

  const KmDF = new Datasheet(Kms, ["Parameter", "Value", "Units"]);

  return { KcatDF, KmDF };
}

export function defineMetabolicReactions() {
  const modelSBML = data.docSBML.model;

  let speciesNames = modelSBML.species.map(species => species.name);
  let speciesNamesLower = speciesNames.map(name => name.toLowerCase());
  let speciesIDs = modelSBML.species.map(species => species.id);
  let rxnNamesSBML = modelSBML.reactions.map(reaction => reaction.name);

  // Read SBTab tables that will be used in the kinetic model
  let RxnDF = data.ComDF_cent.get("Reaction");
  let ParDF = data.ComDF_cent.get("Parameter");
  let QntDF = data.ComDF_cent.get("Quantity");
  let CentQntDF = QntDF;
  let ComDF = data.ComDF_cent.get("Compound");

  // Initialize dictionary for "quantity" values
  const QntDFDict = new Map();

  // Reaction IDs for central metabolism reactions
  const cntrMetRxn = ["PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK",
    "LDH_L", "PDH_acald", "PDH_E3", "PTAr", "ACKr", "NOX", "TALA", "TKT1",
    "TKT2", "RPE", "RPI", "PRPPS", "PPM", "PPM2", "DRPA", "GAPDP"]; //"PDH_E1","PDH_E2" //"PDH_E1","PDH_E2", "GAPDP",

  const cntrMetRxnID = cntrMetRxn.map(rxn => "R_" + rxn);

  const RxnDF_CentralMet = RxnDF.filter(row => cntrMetRxnID.includes(row.Reaction));

  const centralMets = ["g6p_c", "f6p_c", "man6p_c", "gam6p_c", "acgam6p_c", "acmanap_c", "pi_c", "amp_c", "adp_c", "atp_c", "acmana_c",
    "nh3_c", "pep_c", "pyr_c", "xu5p__D_c", "ru5p__D_c", "r5p_c", "prpp_c", "e4p_c", "r1p_c", "fdp_c",
    "s7p_c", "dhap_c", "g3p_c", "2dr5p_c", "2dr1p_c", "nad_c", "nadh_c", "13dpg_c", "nadp_c", "nadph_c",
    "3pg_c", "2pg_c", "lac__L_c", "co2_c", "acdhlpl_PdhC_c", "coa_c", "accoa_c", "actp_c", "ac_c", "dhlpl_PdhC_c",
    "lpl_PdhC_c", "acald_c", "o2_c", "co2_c", "10fthfglu3_c", "thfglu3_c"];

  const concentrationList = [];
  const metIDcheck = [];

  for (let met of centralMets) {
    const metID = "M_" + met;
    const conc = ComDF.getRows().find(row => row.Compound === metID).InitialConcentration;
    concentrationList.push([metID, conc]);
    metIDcheck.push(metID);
  }

  const { KcatDF: CentKcatDF, KmDF: CentKmDF } = getKcat(cntrMetRxnID, QntDF, ParDF, CentQntDF, RxnDF);

  // Amino Acid Metabolism

  const aaMetRxn = ["FMETTRS"];
  const aaMetRxnID = aaMetRxn.map(rxn => "R_" + rxn);

  const RxnDF_aaMet = RxnDF.filter(row => aaMetRxnID.includes(row.Reaction));

  RxnDF = RxnDF_CentralMet.concat(RxnDF_aaMet);

  const { KcatDF: AAKcatDF, KmDF: AAKmDF } = getKcat(aaMetRxnID, QntDF, ParDF, CentQntDF, RxnDF);

  const aaMetIDs = ["M_ala__L_c", "M_arg__L_c",
    "M_asn__L_c", "M_asp__L_c", "M_cys__L_c", "M_glu__L_c", "M_gln__L_c", "M_gly_c",
    "M_his__L_c", "M_ile__L_c", "M_leu__L_c", "M_lys__L_c", "M_met__L_c", "M_phe__L_c",
    "M_pro__L_c", "M_ser__L_c", "M_thr__L_c", "M_trp__L_c", "M_tyr__L_c", "M_val__L_c", "M_alatrna_c", "M_argtrna_c",
    "M_asntrna_c", "M_asptrna_c", "M_cystrna_c", "M_glutrna_c", "M_glntrna_c", "M_glytrna_c",
    "M_histrna_c", "M_iletrna_c", "M_leutrna_c", "M_lystrna_c", "M_mettrna_c", "M_phetrna_c",
    "M_protrna_c", "M_sertrna_c", "M_thrtrna_c", "M_trptrna_c", "M_tyrtrna_c", "M_valtrna_c",
    "M_trnaala_c", "M_trnaarg_c",
    "M_trnaasn_c", "M_trnaasp_c", "M_trnacys_c", "M_trnaglu_c", "M_trnagln_c", "M_trnagly_c",
    "M_trnahis_c", "M_trnaile_c", "M_trnaleu_c", "M_trnalys_c", "M_trnamet_c", "M_trnaphe_c",
    "M_trnapro_c", "M_trnaser_c", "M_trnathr_c", "M_trnatrp_c", "M_trnatyr_c", "M_trnaval_c",
    "M_glutrnagln_c", "M_amet_c", "M_alatrna_c", "M_argtrna_c",
    "M_asntrna_c", "M_asptrna_c", "M_cystrna_c", "M_glutrna_c", "M_glntrna_c", "M_glytrna_c",
    "M_histrna_c", "M_iletrna_c", "M_leutrna_c", "M_lystrna_c", "M_mettrna_c", "M_phetrna_c",
    "M_protrna_c", "M_sertrna_c", "M_thrtrna_c", "M_trptrna_c", "M_tyrtrna_c", "M_valtrna_c"];

  const aaMetDict = new Map();
  for (let met of aaMetIDs) {
    const metID = met;
    const conc = ComDF.findRow("Compound", metID).InitialConcentration;
    aaMetDict.set(met, conc);
    concentrationList.push([metID, conc]);
    metIDcheck.push(metID);
  }

  // Nucleotide Metabolism

  let ParDF_nuc = data.ComDF_nuc.get("Parameter");
  let QntDF_nuc = data.ComDF_nuc.get("Quantity");
  let ComDF_nuc = data.ComDF_nuc.get("Compound");

  ComDF = ComDF.concat(ComDF_nuc);
  let QntDF_Nuclt = QntDF_nuc;

  let RxnDF_Nuclt = data.ComDF_nuc.get("Reaction");

  let nuclMetRxn = data.nuclMetRxn.getRows().map(row => row[0]);
  let nuclMetRxnID = nuclMetRxn.map(x => "R_" + x);

  const nuclTurnOff = ['NTD9', 'DCDPMP', 'DCTPMP', 'DCTPDP', 'CTPDP']; // NTD1, NTD5, NTD6, NTD8
  for (let rxn of nuclTurnOff) {
    const rxnID = 'R_' + rxn;
    nuclMetRxnID = nuclMetRxnID.filter(rxnIDR => rxnIDR !== rxnID);
  }

  let indx = 0;
  for (let rxnID of nuclMetRxnID) {
    if (rxnID === 'R_PYK') {
      nuclMetRxnID.splice(indx, 1);
    } else {
      indx = indx + 1;
    }
  }
  indx = 0;
  for (let rxnID of nuclMetRxnID) {
    if (rxnID === 'R_PGK') {
      nuclMetRxnID.splice(indx, 1);
    } else {
      indx = indx + 1;
    }
  }

  RxnDF_Nuclt = RxnDF_Nuclt.filter(row => nuclMetRxnID.includes(row.Reaction));

  RxnDF = RxnDF_CentralMet.concat(RxnDF_aaMet, RxnDF_Nuclt);

  const nuclMets = ["cytd_c", "cmp_c", "cdp_c", "ctp_c", "dcyt_c", "dcmp_c", "dcdp_c", "dctp_c",
    "gsn_c", "gua_c", "gmp_c", "gdp_c", "gtp_c", "dgsn_c", "dgmp_c", "dgdp_c", "dgtp_c",
    "ura_c", "uri_c", "ump_c", "udp_c", "utp_c", "duri_c", "dump_c", "dudp_c", "dutp_c",
    "adn_c", "ade_c", "dad_2_c", "damp_c", "dadp_c", "datp_c", "thymd_c", "dtmp_c", "dtdp_c", "dttp_c",
    "glu__L_c", "gln__L_c", "ppi_c", "trdrd_c"];

  for (let met of nuclMets) {
    const metID = "M_" + met;
    const conc = ComDF_nuc.findRow("Compound", metID).InitialConcentration;
    concentrationList.push([metID, conc]);
    metIDcheck.push(metID);
  }

  let concDF = new Datasheet(concentrationList, ["Metabolite", "InitialConcentration"]);

  const { KcatDF: NucKcatDF, KmDF: NucKmDF } = getKcat(nuclMetRxnID, QntDF_nuc, ParDF_nuc, QntDF_Nuclt, RxnDF);

  let KmDF = CentKmDF.concat(NucKmDF);

  const pykRate = 3204; // Rate kcat 1/s
  const pgkRate = 220; // Rate kcat 1/s

  const editRxnsDict = {
    'R_PYK': pykRate * 1, 'R_PYK2': pykRate * .11, 'R_PYK3': pykRate * .21,
    'R_PYK4': pykRate * 0.17, 'R_PYK5': pykRate * .18, 'R_PYK6': pykRate * 0.26,
    'R_PYK7': pykRate * 0.13, 'R_PYK8': pykRate * .05, 'R_PYK9': pykRate * .01,
    'R_PGK': pgkRate * 1, 'R_PGK2': pgkRate * 0.17, 'R_PGK3': pgkRate * 0.64,
    'R_PGK4': pgkRate * 0.16
  }; //scaling values from Pollack et al. OMICS 2002

  const tmp = NucKcatDF.getColumn("Value");
  const tmpRxns = ['R_PYK2', 'R_PYK3', 'R_PYK4', 'R_PYK5', 'R_PYK6', 'R_PYK7', 'R_PYK8', 'R_PYK9',
    'R_PGK2', 'R_PGK3', 'R_PGK4'];

  for (let rxn of tmpRxns) {
    tmp[NucKcatDF.findRowIndex("Parameter", 'kcatF_' + rxn)] = editRxnsDict[rxn];
  }

  NucKcatDF.setColumn("Value", tmp);

  let KcatDF = CentKcatDF.concat(NucKcatDF);

  // Lipid Metabolism

  const ParDF_lip = data.ComDF_lip.get("Parameter");
  const QntDF_lip = data.ComDF_lip.get("Quantity");
  const ComDF_lip = data.ComDF_lip.get("Compound");

  ComDF = ComDF.concat(ComDF_lip);
  let QntDF_Lipid = QntDF_lip;
  let RxnDF_Lipid = data.ComDF_lip.get("Reaction");

  const lipMetRxn = ["GLYK", "ACPS", "BPNT", "FAKr", "ACPPAT", "APG3PAT", "AGPAT", "DASYN", "PGSA", "PGPP", "CLPNS", "PGMT", "PAPA", "GALU", "UDPG4E", "UDPGALM", "DAGGALT", "PSSYN", "DAGPST"];

  const lipTurnOff = ["DAGPST", "PSSYN"];
  for (let rxn of lipTurnOff) {
    lipMetRxn.splice(lipMetRxn.indexOf(rxn), 1);
  }

  const lipMetRxnID = lipMetRxn.map(rxn => "R_" + rxn);

  RxnDF_Lipid = RxnDF_Lipid.filter(row => lipMetRxnID.includes(row.Reaction));

  RxnDF = RxnDF_CentralMet.concat(RxnDF_aaMet, RxnDF_Nuclt, RxnDF_Lipid);

  const lipMets = ["udpgal_c", "fa_c", "ap_c", "pa_c", "glyc_c", "glyc3p_c", "apoACP_c", "coa_c", "pap_c",
    "ACP_c", "ACP_R_c", "1ag3p_c", "cdpdag_c", "pg3p_c", "pg_c", "clpn_c", "12dgr_c", "galfur12dgr_c",
    "udpgalfur_c", "glyc_c", "chsterol_c"];

  for (let met of lipMets) {
    const metID = "M_" + met;
    const conc = ComDF_lip.findRow("Compound", metID).InitialConcentration;
    concentrationList.push([metID, conc]);
    metIDcheck.push(metID);
  }

  concDF = new Datasheet(concentrationList, ["Metabolite", "InitialConcentration"]);

  const { KcatDF: LipKcatDF, KmDF: LipKmDF } = getKcat(lipMetRxnID, QntDF_lip, ParDF_lip, QntDF_Lipid, RxnDF);

  RxnDF = data.ComDF_cent.get("Reaction");
  ParDF = data.ComDF_cent.get("Parameter");
  QntDF = data.ComDF_cent.get("Quantity");
  CentQntDF = QntDF;
  ComDF = data.ComDF_cent.get("Compound");

  const cofactMetRxn = ["NCTPPRT", "NNATr", "NADS", "NADK", "RBFK", "FMNAT", "5FTHFPGS", "FTHFCL", "MTHFC", "GHMT", "MTHFD"]; //,"NADHK","GHMT2"]

  const cofactMetRxnID = cofactMetRxn.map(rxn => "R_" + rxn);

  const RxnDF_cofactMet = RxnDF.filter(row => cofactMetRxnID.includes(row.Reaction));

  RxnDF = RxnDF_CentralMet.concat(RxnDF_aaMet, RxnDF_Nuclt, RxnDF_Lipid, RxnDF_cofactMet);

  const { KcatDF: cofactKcatDF, KmDF: cofactKmDF } = getKcat(cofactMetRxnID, QntDF, ParDF, CentQntDF, RxnDF);

  const cofactMetIDs = ["M_nac_c", "M_nicrnt_c", "M_dnad_c", "M_nh3_c", "M_ribflv_c", "M_fmn_c", "M_fad_c", "M_sprm_c", "M_thmpp_c",
    "M_5fthf_c", "M_5fthfglu3_c", "M_methfglu3_c", "M_10fthfglu3_c", "M_thfglu3_c", "M_mlthfglu3_c"];

  const cofactMetDict = new Map();
  for (let met of cofactMetIDs) {
    const metID = met;
    const conc = ComDF.findRow("Compound", metID).InitialConcentration;
    cofactMetDict.set(met, conc);
    concentrationList.push([metID, conc]);
    metIDcheck.push(metID);
  }

  concDF = new Datasheet(concentrationList, ["Metabolite", "InitialConcentration"]);

  KcatDF = CentKcatDF.concat(NucKcatDF, AAKcatDF, LipKcatDF, cofactKcatDF);
  KmDF = CentKmDF.concat(NucKmDF, AAKmDF, LipKmDF, cofactKmDF);

  // Transport Reactions

  ComDF = ComDF.concat(data.transport.get("Compound"));
  QntDF = QntDF.concat(data.transport.get("Quantity"));

  let RxnDF_Transport = data.transport.get("Reaction");

  // WARNING ##
  // We are missing an ATPase reaction, so there will be no explicit tracking of proton concentration.
  // Therefore, we are removing the proton transport reaction.

  RxnDF_Transport = RxnDF_Transport.filter(row => row.Reaction !== "R_Ht");

  // Concatonate the Central and transport reaction dataframes.

  RxnDF = RxnDF_CentralMet.concat(RxnDF_Nuclt, RxnDF_Lipid, RxnDF_aaMet, RxnDF_Transport, RxnDF_cofactMet);

  // Read quantities in from the quantities dataframe.
  for (let items of QntDF.getRows()) {
    if (items.Quantity) {
      if (items.Quantity.length == 1) {
        continue;
      }
      if (!QntDFDict.has(items.Quantity)) {
        QntDFDict.set(items.Quantity, parseFloat(items.Value));
      }
    }
  }

  const AvailQnts = Array.from(QntDFDict.keys());

  // Adapt model to specialized transport rate forms for lipid and glucose transport.

  // List of reactions that need to be removed from the model in
  // order to be re-defined. The Reaction Formulas need to be redefined.

  let transpRxnsList = ["FAt", "CHOLt", "GLCpts"];

  for (let rxnName of transpRxnsList) {
    modelSBML.removeReaction(rxnName);
  }

  // List of reactions that need to be created in the SBML model.
  transpRxnsList = ["GLCpts0", "GLCpts1", "GLCpts2", "GLCpts3", "GLCpts4", "GLYCt", "FAt", "CHOLt"];

  let newRxnsDF = RxnDF_Transport.filter(row => transpRxnsList.includes(row.Name));

  let rxnCounter = 1;

  for (let items of newRxnsDF.getRows()) {
    // Create new reaction
    let rxnTmp = modelSBML.createReaction("R_" + items.Name, items.Name, true);

    // Extract reactants and products from the reaction formula
    let substratesList = items.ReactionFormula.split("<=>")[0];
    let productsList = items.ReactionFormula.split("<=>")[1];

    substratesList = substratesList.split("+");
    substratesList = substratesList.map(met => met.trim());

    productsList = productsList.split("+");
    productsList = productsList.map(met => met.trim());

    // Add substrates to the reaction object
    for (let spcTmpID of substratesList) {

      // If the species does not already exist in the SBML model, create it.
      if (!speciesIDs.includes(spcTmpID)) {
        modelSBML.createSpecies(spcTmpID, spcTmpID, 'c', false, 0, 'mmole')
      }

      rxnTmp.listOfReactants.push({
        species: spcTmpID,
        stoichiometry: 1,
        constant: false
      });
    }

    // Add products to the reaction object
    for (let spcTmpID of productsList) {

      // If the species does not already exist in the SBML model, create it.
      if (!speciesIDs.includes(spcTmpID)) {
        modelSBML.createSpecies(spcTmpID, spcTmpID, 'c', false, 0, 'mmole')
      }

      rxnTmp.listOfProducts.push({
        species: spcTmpID,
        stoichiometry: 1,
        constant: false
      });
    }

    rxnCounter += 1;
  }

  // Reset variables with new species and reactions
  speciesNames = modelSBML.species.map(species => species.name);
  speciesNamesLower = speciesNames.map(name => name.toLowerCase());
  speciesIDs = modelSBML.species.map(species => species.id);
  rxnNamesSBML = modelSBML.reactions.map(reaction => reaction.name);

  return {
    modelSBML,
    RxnDF, ComDF, concDF,
    cntrMetRxn, nuclMetRxn, cofactMetRxn, lipMetRxn, aaMetRxn,
    metIDcheck, QntDFDict, AvailQnts,
    KcatDF, KmDF
  };
}

type MetabolicData = ReturnType<typeof defineMetabolicReactions>;

const transportRxns = ['GLCpts0', 'GLCpts1', 'GLCpts2', 'GLCpts3', 'GLCpts4', 'ACt',
  'PYRt2r', 'L_LACt2r', 'GLYCt', 'FAt', 'CHOLt'];

function getSpecIDs(rxnName: string, modelSBML: SBMLModel) {
  const rxnObj = modelSBML.getReaction(rxnName);

  // Use model SBML to get IDs, names, and stoichiometries for reactants
  const specIDs = rxnObj.listOfReactants.map(x => x.species);
  const spcStoich = rxnObj.listOfReactants.map(x => -1 * x.stoichiometry);
  const spcNames = specIDs.map(spcID => modelSBML.getSpecies(spcID).name);

  if (spcStoich.some(isNaN)) {
    throw new Error(`Invalid stoichiometry for reaction: ${rxnName}`);
  }

  const returnList: [string[], string[], number[]][] = [];
  returnList.push([spcNames, specIDs, spcStoich]);

  // Now do the same for products
  const specIDs2 = rxnObj.listOfProducts.map(x => x.species);
  const spcStoich2 = rxnObj.listOfProducts.map(x => x.stoichiometry);
  const spcNames2 = specIDs2.map(spcID => modelSBML.getSpecies(spcID).name);

  if (spcStoich2.some(isNaN)) {
    throw new Error(`Invalid stoichiometry for reaction: ${rxnName}`);
  }

  returnList.push([spcNames2, specIDs2, spcStoich2]);

  return returnList;
}

// Add metabolites to ODE model.

const constMetList = [];

function addSpcToODE(
  spcID: string,
  concDF: Datasheet,
  model: ODE,
  pmap: Map<string, number>,
  constantVals: string[],
  explicitParsDict: Map<string, ODEParameterValue>
) {
  // If the species was not defined in reaction input files, return error code
  if (!concDF.getColumn("Metabolite").includes(spcID)) {
    console.log(`ERROR: species ID \"${spcID}\" not found!`);
    return 1;
  }

  // If species has already been added to model, return
  if (model.metabolites.has(spcID)) {
    console.log(`WARNING: species ID \"${spcID}\" already added to model!`);
    return 0;
  }

  const spcName = concDF.findRow("Metabolite", spcID).Metabolite;

  if (!constantVals.includes(spcID)) {

    // If it is an external concentration, not previously set, add it as a constant
    if (spcID.endsWith("_e")) {
      constantVals.push(spcID);
      const spcConc = concDF.findRow("Metabolite", spcID).InitialConcentration;
      model.addParameter("Explicit", spcID, spcConc, 0, 0, "mM", `External concentration (${spcID})`);
    } else {
      const spcConc = partTomM(pmap.get(spcID), pmap);
      model.addMetabolite(spcID, spcName, spcConc);
    }
  } else {
    // For escher maps
    constMetList.push(spcID);
    const spcConc = concDF.findRow("Metabolite", spcID).InitialConcentration;
    explicitParsDict.set(spcID, spcConc);
  }

  return 0;
}

function addExtSpcToODE(
  spcID: string,
  ComDF: Datasheet,
  model: ODE,
  pmap: Map<string, number>,
  constantVals: string[],
  explicitParsDict: Map<string, ODEParameterValue>
) {
  // If the species was not defined in reaction input files, return error code
  if (!ComDF.getColumn("Compound").includes(spcID)) {
    console.log(`ERROR: species ID \"${spcID}\" not found!`);
    return 1;
  }

  // If species has already been added to model, return
  if (model.metabolites.has(spcID)) {
    console.log(`WARNING: species ID \"${spcID}\" already added to model!`);
    return 0;
  }

  const spcName = ComDF.findRow("Compound", spcID).Name;

  if (!constantVals.includes(spcID)) {

    // If it is an external concentration, not previously set, add it as a constant
    if (spcID.endsWith("_e")) {
      constantVals.push(spcID);
      const spcConc = ComDF.findRow("Compound", spcID).InitialConcentration;
      model.addParameter("Explicit", spcID, spcConc, 0, 0, "mM", `External concentration (${spcID})`);
    } else {
      const spcConc = partTomM(pmap.get(spcID), pmap);
      model.addMetabolite(spcID, spcName, spcConc);
    }
  } else {
    // For escher maps
    constMetList.push(spcID);
    const spcConc = partTomM(pmap.get(spcID), pmap);
    explicitParsDict.set(spcID, spcConc);
  }

  return 0;
}

function createGeneExpression(rxnID: string, pmap: Map<string, number>) {
  // Creates protein species and auxiliary functions for effective protein levels.
  // Creates necessary trnscription and translation reactions, when info is available.

  // List of reaction IDs in *metabolic* reconstruction
  const reconstRxnIDList = data.reconstPD.getColumn("Reaction ID");

  if (!reconstRxnIDList.includes(rxnID)) {
    console.log(`WARNING!! Reaction ${rxnID} not found in reconstruction!`);
    return { enzLvlVar: '', enzConc: defaultPtnConcentration };
  }

  const gprStr = data.reconstPD.findRow("Reaction ID", rxnID)["GPR rule"];

  // Checks if GPR tule is available.
  if (gprStr === undefined) {
    //         print("\n\tReaction", rxnID, "has no GPR!!\n")
    // Keeps default value of 0.001 mM enzyme concentration.
    return { enzLvlVar: '', enzConc: defaultPtnConcentration };
  }

  // Get list of genes in MM* code.
  const genes = gprStr.match(/\bMM[A-Z0-9_]{1,}\b/g);

  // Loop over genes associated with the reaction
  const rxnPtns = [];
  for (let mmcode of genes) {
    const locusNum = mmcode.split('SYN1_')[1];
    const jcvi3AID = 'JCVISYN3A_' + locusNum;

    const newMetID = "M_PTN_" + jcvi3AID + "_c";

    // Keep list of protein species associated with the current Rxn.
    rxnPtns.push(newMetID);
  }

  let enzLvlVar = "---";
  let enzConc = 0;

  // Check if this reaction has multiple genes.
  if (rxnPtns.length === 1) {
    // If there is just one protein, this is used directly in the reaction's rate law
    enzLvlVar = rxnPtns[0];
    enzConc = partTomM(pmap.get(rxnPtns[0]), pmap);
  } else {
    // IF there are multiple proteins and a logical rule, we need an extra variable
    // which adds the proteins' concentrations (OR rule) or get the minimal concentration (AND rule).

    // To create an "effective" enzyme level, the model will calculate, at every step, the "sum" or "min" 
    //   of multiple protein concentrations.
    enzLvlVar = rxnID + "_effecEnz";

    if (gprStr.toLowerCase().includes(" or ")) {
      let enzCnts = 0;
      for (let i = 0; i < rxnPtns.length; i++) {
        const cnt = pmap.get(rxnPtns[i]);
        enzCnts = enzCnts + cnt;
      }
      enzConc = partTomM(enzCnts, pmap);
    } else {
      const enzCnts = [];
      for (let i = 0; i < rxnPtns.length; i++) {
        const cnt = pmap.get(rxnPtns[i]);
        enzCnts.push(cnt);
      }
      const enzCnt = Math.min(...enzCnts);
      enzConc = partTomM(enzCnt, pmap);
    }
  }

  return { enzLvlVar, enzConc };
}

function addGlobalParams(model: ODE) {
  // Parse global parameters to model
  model.parseParameters(data.modelParameters);

  const explPar = model.explicitParameters;

  // Store explicit and global parameters (fixed values used throughout the model)
  // They may contain metabolites with constant concentration, from the external media
  const explicitParsDict = new Map<string, ODEParameterValue>();
  for (let par of explPar.values()) {
    explicitParsDict.set(par.formKey, par.value);
  }

  const constantVals = Array.from(explicitParsDict.keys());

  return { constantVals, explicitParsDict };
}

export function addReactionsToModel(model: ODE, pmap: Map<string, number>, data: MetabolicData) {
  const {
    modelSBML,
    RxnDF, ComDF, concDF,
    cntrMetRxn, nuclMetRxn, cofactMetRxn, lipMetRxn, aaMetRxn,
    metIDcheck, QntDFDict, AvailQnts,
    KcatDF, KmDF
  } = data;

  // To store the rate forms and their names
  const rateFormsDict = {};
  const knownForms = new Set<string>();

  // Diagnostics:
  // parameters with *zero/inf/nan* value:
  const errorPar = [];

  let rxnCounter = 1;

  const { constantVals, explicitParsDict } = addGlobalParams(model);

  for (const items of RxnDF.getRows()) {

    if (nuclMetRxn.includes(items.Name) || cntrMetRxn.includes(items.Name) || lipMetRxn.includes(items.Name) || aaMetRxn.includes(items.Name) || cofactMetRxn.includes(items.Name)) {

      // For clarity
      const rxnID = items.Name;

      const { enzConc } = createGeneExpression(rxnID, pmap);

      const specLists = getSpecIDs(rxnID, modelSBML);

      const reactantsDict = {};

      // Matches the substrates in the rate form to standardized IDs.
      const substratesList = [];
      let substrates = 0;
      let spcCounter = 1;

      for (let i = 0; i < specLists[0][1].length; i++) {

        const spc = specLists[0][1][i];
        const spcStoich = specLists[0][2][i];

        let subsStr = "";

        for (let j = 0; j < Math.floor(Math.abs(spcStoich)); j++) {
          subsStr = "Sub" + spcCounter;
          substratesList.push([spc, subsStr, spcStoich]);
          spcCounter += 1;
        }

        // we create this lookup dict to standardize parameter names later on
        reactantsDict[spc] = subsStr;

        for (let j = 0; j < Math.floor(Math.abs(spcStoich)); j++) {
          substrates = substrates + 1;
        }
      }

      // Matches the products in the rate form to standardized IDs.
      const productsList = [];

      spcCounter = 1;
      let products = 0;

      for (let i = 0; i < specLists[1][1].length; i++) {

        const spc = specLists[1][1][i];
        const spcStoich = specLists[1][2][i];

        let prodStr = "";

        for (let j = 0; j < Math.floor(Math.abs(spcStoich)); j++) {
          prodStr = "Prod" + spcCounter;
          productsList.push([spc, prodStr, spcStoich]);
          spcCounter += 1;
        }

        reactantsDict[spc] = prodStr;

        for (let j = 0; j < Math.floor(Math.abs(spcStoich)); j++) {
          products = products + 1;
        }
      }

      const rxnName = 'R_' + rxnID;

      const rateLaw = enzymatic(substrates, products);

      const RateName = rxnID + '_Rate';

      model.addRateForm(RateName, rateLaw);

      model.addReaction(rxnID, RateName, "Reaction " + rxnID);

      for (let spc of Object.keys(reactantsDict)) {
        if (model.metabolites.has(spc)) {
          // In case we already know about this metabolite.
          continue;
        }
        if (metIDcheck.includes(spc)) {
          if (addSpcToODE(spc, concDF, model, pmap, constantVals, explicitParsDict)) {
            // In case the metabolite ID referenced by a rate form is not found
            //  among defined metabolites.
            break;
          }
        } else {
          if (addExtSpcToODE(spc, ComDF, model, pmap, constantVals, explicitParsDict)) {
            // In case the metabolite ID referenced by a rate form is not found
            //  among defined metabolites.
            break;
          }
        }
      }

      const rxnMetsAdded = [];

      for (let [modelID, stdID, spcStoich] of substratesList) {
        if (constantVals.includes(modelID)) {
          model.addParameter(rxnID, stdID, modelID);
        } else if (rxnMetsAdded.includes(modelID)) {
          model.addParameter(rxnID, stdID, modelID);
        } else {
          model.addSubstrate(rxnID, stdID, modelID, spcStoich);
          rxnMetsAdded.push(modelID);
        }

        const km_ID = 'kmc_R_' + rxnID + '_' + modelID + '_R_' + rxnID;

        const km_Value = KmDF.findRow("Parameter", km_ID).Value;

        const KM_rxnID = "Km" + stdID.replace('$', '');

        model.addParameter(rxnID, KM_rxnID, km_Value);
      }

      for (let [modelID, stdID, spcStoich] of productsList) {
        if (constantVals.includes(modelID)) {
          model.addParameter(rxnID, stdID, modelID);
        } else if (rxnMetsAdded.includes(modelID)) {
          model.addParameter(rxnID, stdID, modelID);
        } else {
          model.addProduct(rxnID, stdID, modelID, spcStoich);
          rxnMetsAdded.push(modelID);
        }

        const km_ID = 'kmc_R_' + rxnID + '_' + modelID + '_R_' + rxnID;

        const km_Value = KmDF.findRow("Parameter", km_ID).Value;

        const KM_rxnID = "Km" + stdID.replace('$', '');

        model.addParameter(rxnID, KM_rxnID, km_Value);
      }

      const kcatFID = "kcatF_" + rxnName;
      const kcatF = KcatDF.findRow("Parameter", kcatFID).Value;

      model.addParameter(rxnID, 'kcatF', kcatF);

      const kcatRID = "kcatR_" + rxnName;
      const kcatR = KcatDF.findRow("Parameter", kcatRID).Value;

      model.addParameter(rxnID, 'kcatR', kcatR);

      // Add enzyme level parameter.
      if (model.getReaction(rxnID).getKeys().has("Enzyme")) {
        model.addParameter(rxnID, "Enzyme", enzConc, 0, 0, "mM", "Enzyme concentration");
      }

      if (model.getReaction(rxnID).getKeys().has("onoff")) {
        model.addParameter(rxnID, "onoff", 1, 0, 1, "mM", "Debug On/Off switch");
      }

      rxnCounter += 1;
    } else if (transportRxns.includes(items.Name)) {
      let kl = items.KineticLaw.replace(/\s{2,}/g, ' ');

      // For clarity
      const rxnID = items.Name;

      const specLists = getSpecIDs(rxnID, modelSBML);

      const reactantsDict = {};

      // Matches the substrates in the rate form to standardized IDs.
      const substratesList: [string, string, number][] = [];
      let spcCounter = 1;
      for (let i = 0; i < specLists[0][1].length; i++) {

        const spc = specLists[0][1][i];
        const spcStoich = specLists[0][2][i];

        const subsStr = "Sub" + spcCounter;

        substratesList.push([spc, "$" + subsStr, spcStoich]);

        // we create this lookup dict to standardize parameter names later on
        reactantsDict[spc] = subsStr;

        // We use regular expression syntax to assure that no substrings will
        // be matched. For example, when substituting the species M_abc_c for $Sub1, 
        // we would "sub-match" the parameter kmc_R_XYZ_M_abc_c, and create kmc_R_XYZ_Sub1,
        // which would be bad...
        kl = kl.replace(new RegExp(`\\b${spc}\\b`, 'g'), "$" + subsStr);

        spcCounter += 1;
      }

      // Matches the products in the rate form to standardized IDs.
      const productsList: [string, string, number][] = [];
      spcCounter = 1;
      for (let i = 0; i < specLists[1][1].length; i++) {

        const spc = specLists[1][1][i];
        const spcStoich = specLists[1][2][i];

        const prodStr = "Prod" + spcCounter;

        productsList.push([spc, "$" + prodStr, spcStoich]);
        reactantsDict[spc] = prodStr;

        kl = kl.replace(new RegExp(`\\b${spc}\\b`, 'g'), "$" + prodStr);

        spcCounter += 1;
      }

      // Break main loop
      let metError = false;

      // Sanity Checks that the metabolites found here are in the ODECell model.
      for (let spc of Object.keys(reactantsDict)) {
        if (model.metabolites.has(spc)) {
          // In case we already know about this metabolite.
          continue;
        }
        if (metIDcheck.includes(spc)) {
          if (addSpcToODE(spc, concDF, model, pmap, constantVals, explicitParsDict)) {
            // In case the metabolite ID referenced by a rate form is not found
            //  among defined metabolites.
            metError = true;
            break;
          }
        } else {
          if (addExtSpcToODE(spc, ComDF, model, pmap, constantVals, explicitParsDict)) {
            // In case the metabolite ID referenced by a rate form is not found
            //  among defined metabolites.
            metError = true;
            break;
          }
        }
      }

      if (metError) {
        break;
      }

      const paramsDict = {};
      const paramCount = {};

      // Substitute parameter values
      for (let qnt of AvailQnts) {
        // Checks that the parameter refers to the reaction
        // Each parameter has the "R_XXX" reaction code in its name.
        // Needed because "R_AB1" parameters would be found in "R_AB12" laws.
        // Check if the quantity is an explicit paramter (a global constant)
        //   shared between rate laws.
        if (!qnt.includes(items.Reaction) || constantVals.includes(qnt)) {
          continue;
        }

        // Sanity Checks if the parameter is found in the kinetic law
        if (!kl.includes(qnt)) {
          continue;
        }

        let paramType = "";

        // Classify the parameter, store ID, name and value.
        if (qnt.includes("hco")) {
          paramType = "hco";
          paramCount["hco"] = (paramCount["hco"] || 0) + 1;
          paramsDict[qnt] = ["$hco" + paramCount["hco"], "reaction cooperativity"];
        } else if (qnt.includes("keq")) {
          paramType = "keq";
          paramCount["keq"] = (paramCount["keq"] || 0) + 1;
          paramsDict[qnt] = ["$keq" + paramCount["keq"], "equilibrium constant"];
        } else if (qnt.includes("kcr")) {
          paramType = "kcr";
          paramCount["kcr"] = (paramCount["kcr"] || 0) + 1;
          paramsDict[qnt] = ["$kcat" + paramCount["kcr"], "catalitic rate constant"];
        } else if (qnt.includes("km")) {
          paramType = "km";
          paramCount["km"] = (paramCount["km"] || 0) + 1;

          // Find species related to the Km value
          const relatedSpcs = new Set<string>();
          for (let spc of [...substratesList, ...productsList].map(x => x[0])) {
            if (spc.includes(qnt)) {
              relatedSpcs.add(spc);
            }
          }

          if (!relatedSpcs.size) {
            console.log("WARNING: no species found for parameter", qnt);
          } else {
            paramsDict[qnt] = ["$km_" + [...relatedSpcs].map(x => reactantsDict[x]).join("_"), "species constant"];
          }
        } else {
          paramType = "OTHER";
          paramCount["OTHER"] = (paramCount["OTHER"] || 0) + 1;

          paramsDict[qnt] = ["$Const" + paramCount["OTHER"], "rate law-specific constant"];
        }

        // Sanity check:
        if (QntDFDict.get(qnt) === 0 || isNaN(QntDFDict.get(qnt)) || isFinite(QntDFDict.get(qnt))) {
          errorPar.push([qnt, paramType, QntDFDict.get(qnt)]);
          if (qnt.includes('keq')) {
            QntDFDict.set(qnt, 10 ** -7);
          }
        }

        // Use regular expression to extract *only* the parameter ID and substitute it
        // with a standardized ID.
        kl = kl.replace(new RegExp(`\\b${qnt}\\b`, 'g'), paramsDict[qnt][0]);
      }

      let enzLevlPar = "";
      // Create argument for protein level.
      if (kl.startsWith("( 0.001 / 1.0)")) {
        // If the search returns anything but a "None", it found a match.
        // We need the test to create a reaction-specific parameter later on.
        // Creates key for rate law.
        enzLevlPar = "$enzLvl";

        // Substitute hard-coded default value for key in rate law.
        kl = kl.replace("( 0.001 / 1.0)", enzLevlPar);
      }

      // DEBUG -- create an ON/OFF switch for all metabolic reactions:
      kl = "$onoff * (" + kl + " )";

      if (!knownForms.has(kl)) {
        // If we found a new rate form, report on it and store it for easy search.
        rateFormsDict[kl] = ["RateForm" + knownForms.size, `${substratesList.length}s_${productsList.length}p`];
        knownForms.add(kl);

        // Now add the rate form to the model.
        model.addRateForm(rateFormsDict[kl][0], kl);
      }

      let enzLvlVar = "---";

      if (enzLevlPar) {
        // In case we are explicitly tracking enzyme levels:
        // Add gene expression for necessary enzyme(s) for the current reaction.
        // This has to run *before*  the actual reaction is added in case a GPR rule 
        //   needs an auxiliary reaction to determine effective enzyme levels.
        enzLvlVar = createGeneExpression(rxnID, pmap).enzLvlVar;
      }

      // Add reaction to model, using the rate form defined above
      model.addReaction(rxnID, rateFormsDict[kl][0], "Reaction " + rxnID);

      // Add substrate species
      for (let [modelID, stdID, spcStoich] of substratesList) {
        if (constantVals.includes(modelID)) {
          model.addParameter(rxnID, stdID.slice(1), modelID);
        } else {
          model.addSubstrate(rxnID, stdID.slice(1), modelID, spcStoich);
        }
      }

      // Add product species
      for (let [modelID, stdID, spcStoich] of productsList) {
        if (constantVals.includes(modelID)) {
          model.addParameter(rxnID, stdID.slice(1), modelID);
        } else {
          model.addProduct(rxnID, stdID.slice(1), modelID, spcStoich);
        }
      }

      // Add reaction-specific parameter values
      for (let [modelID, val] of Object.entries(paramsDict)) {
        model.addParameter(rxnID, val[0].slice(1), QntDFDict.get(modelID), 0, 0, "", val[1]);
      }

      // Add enzyme level parameter.
      if (model.getReaction(rxnID).getKeys().has("enzLvl")) {
        model.addParameter(rxnID, "enzLvl", enzLvlVar, 0, 0, "mM", "Enzyme concentration");
      }

      // DEBUG
      if (model.getReaction(rxnID).getKeys().has("onoff")) {
        model.addParameter(rxnID, "onoff", 1, 0, 1, "mM", "Debug On/Off switch");
      }

      rxnCounter += 1;
    }
  }

  const addReaction = (
    name: string, rate: string,
    enzConc: number, kcatF: number, kcatR: number,
    subNames: ODEParameterValue[], subKms: ODEParameterValue[], prodNames: ODEParameterValue[], prodKms: ODEParameterValue[],
  ) => {
    model.addRateForm(name + "_Rate", rate);

    model.addReaction(name, name + "_Rate", name + " transport");
    model.addParameter(name, "Enzyme", enzConc);

    model.addParameter(name, "kcatF", kcatF);
    model.addParameter(name, "kcatR", kcatR);

    for (let i = 0; i < Math.max(subNames.length, subKms.length); i++) {
      if (i < subNames.length) {
        if (typeof subNames[i] === "string") model.addSubstrate(name, "Sub" + (i + 1), subNames[i] as string);
        else if (typeof subNames[i] === "number") model.addParameter(name, "Sub" + (i + 1), subNames[i] as number);
      }
      if (i < subKms.length) model.addParameter(name, "KmSub" + (i + 1), subKms[i]);
    }

    for (let i = 0; i < Math.max(prodNames.length, prodKms.length); i++) {
      if (i < prodNames.length) {
        if (typeof prodNames[i] === "string") model.addProduct(name, "Prod" + (i + 1), prodNames[i] as string);
        else if (typeof prodNames[i] === "number") model.addParameter(name, "Prod" + (i + 1), prodNames[i] as number);
      }
      if (i < prodKms.length) model.addParameter(name, "KmProd" + (i + 1), prodKms[i]);
    }

    model.addParameter(name, "onoff", 1, 0, 1, "mM", "Debug On/Off switch");
  };

  // PUNP5  --- Explicit Addition
  const punp5RateLaw = "$onoff * $Enzyme * (($kcatF*($Sub1/$KmSub1)*($Sub2/$KmSub2)-$kcatR*($Prod1/$KmProd1)*($Prod2/$KmProd2)) / ((1+$Sub1/$KmSub1)*(1+$Sub2/$KmSub2)+(1+$Prod1/$KmProd1)*(1+$Prod2/$KmProd2)-1))";
  const enzConc_PUNP5 = partTomM(pmap.get("M_PTN_JCVISYN3A_0747_c"), pmap); // mM
  model.addMetabolite("M_uri_c", "uridine", partTomM(pmap.get("M_uri_c"), pmap)); // mM - from pb file

  addReaction("PUNP5", punp5RateLaw, enzConc_PUNP5, 26.123325, 109.8933, ["M_uri_c", "M_pi_c"], [0.115, 2.1566], ["M_ura_c", "M_r1p_c"], [2.8972, 0.0137]);

  // ADD ATPase
  const enzConc_ATPase = partTomM(pmap.get("M_PTN_JCVISYN3A_0789_c"), pmap);
  const ATPaseRateLaw = enzymatic(2, 1);

  addReaction("ATPase", ATPaseRateLaw, enzConc_ATPase, 20, 217 / 3, ["M_adp_c", "M_pi_c"], [0.1, 4.2], ["M_atp_c"], [0.6]);

  // ADD Piabc
  const enzConc_PIabc = partTomM(pmap.get("M_PTN_JCVISYN3A_0427_c"), pmap);
  const PitRateLaw = enzymatic(2, 3);

  addReaction("PIabc", PitRateLaw, enzConc_PIabc, 25, 0, [134, "M_atp_c"], [0.0031, 0.023], ["M_pi_c", "M_pi_c", "M_adp_c"], [0.02, 0.385, 0.654]);

  // ADD ADNabc
  const enzConc_NUCabc = partTomM(Math.min(pmap.get("M_PTN_JCVISYN3A_0010_c"), pmap.get("M_PTN_JCVISYN3A_0009_c")), pmap);
  const AdnRateLaw = enzymatic(2, 3);

  addReaction("ADNabc", AdnRateLaw, enzConc_NUCabc, 2, 0, [0.15, "M_atp_c"], [1.9 / 10 ** 3, 0.6], ["M_adn_c", "M_pi_c", "M_adp_c"], [0.02, 10, 2.8]);

  // ADD DADNabc
  const DAdnRateLaw = enzymatic(2, 3);

  addReaction("DADNabc", DAdnRateLaw, enzConc_NUCabc, 0.5, 0, [0.02, "M_atp_c"], [1.9 / 10 ** 3, 0.6], ["M_dad_2_c", "M_pi_c", "M_adp_c"], [0.02, 10, 2.8]);

  // ADD GSNabc
  const GsnRateLaw = enzymatic(2, 3);

  addReaction("GSNabc", GsnRateLaw, enzConc_NUCabc, 1, 0, [0.13, "M_atp_c"], [1.9 / 10 ** 3, 0.6], ["M_gsn_c", "M_pi_c", "M_adp_c"], [0.02, 10, 2.8]);

  // ADD DGSNabc
  const DGsnRateLaw = enzymatic(2, 3);

  addReaction("DGSNabc", DGsnRateLaw, enzConc_NUCabc, 0.7, 0, [0.02, "M_atp_c"], [1.9 / 10 ** 3, 0.6], ["M_dgsn_c", "M_pi_c", "M_adp_c"], [0.02, 10, 2.8]);

  // ADD URIabc
  const UriRateLaw = enzymatic(2, 3);

  addReaction("URIabc", UriRateLaw, enzConc_NUCabc, 2, 0, [0.18, "M_atp_c"], [1.7 / 10 ** 3, 0.6], ["M_uri_c", "M_pi_c", "M_adp_c"], [0.02, 10, 2.8]);

  // ADD DCYTabc
  const DCytRateLaw = enzymatic(2, 3);

  addReaction("DCYTabc", DCytRateLaw, enzConc_NUCabc, 0.5, 0, [0.02, "M_atp_c"], [1.1 / 10 ** 3, 0.6], ["M_dcyt_c", "M_pi_c", "M_adp_c"], [0.02, 10, 2.8]);

  // ADD THMDabc
  const ThmdRateLaw = enzymatic(2, 3);

  addReaction("THMDabc", ThmdRateLaw, enzConc_NUCabc, 0.75, 0, [0.08, "M_atp_c"], [1.7 / 10 ** 3, 0.6], ["M_thymd_c", "M_pi_c", "M_adp_c"], [0.02, 10, 2.8]);

  // ADD SPRMabc
  const enzConc_SPRMabc = partTomM(Math.min(pmap.get("M_PTN_JCVISYN3A_0195_c"), pmap.get("M_PTN_JCVISYN3A_0196_c"), pmap.get("M_PTN_JCVISYN3A_0197_c")), pmap);
  const SprmabcRateLaw = enzymatic(1, 1);
  model.addMetabolite("M_sprm_c", "spermidine", partTomM(pmap.get("M_sprm_c"), pmap));

  addReaction("SPRMabc", SprmabcRateLaw, enzConc_SPRMabc, 3, 0, [0.1, "M_atp_c"], [2, 0.6], ["M_sprm_c", "M_pi_c", "M_adp_c"], [2, 2, 2]);

  // ADD NACabc
  const enzConc_NACabc = partTomM(pmap.get("M_PTN_JCVISYN3A_0314_c"), pmap);
  const NacabcRateLaw = enzymatic(1, 1);

  addReaction("NACabc", NacabcRateLaw, enzConc_NACabc, 1.66, 0, [0.016, "M_atp_c"], [0.0019, 0.023], ["M_nac_c", "M_pi_c", "M_adp_c"], [0.0019, 0.385, 0.654]);

  // ADD RIBFLVabc
  const enzConc_RIBFLVabc = partTomM(Math.min(pmap.get("M_PTN_JCVISYN3A_0641_c"), pmap.get("M_PTN_JCVISYN3A_0642_c"), pmap.get("M_PTN_JCVISYN3A_0643_c"), pmap.get("M_PTN_JCVISYN3A_0877_c")), pmap);
  const RIBFLVabcRateLaw = enzymatic(1, 1);

  addReaction("RIBFLVabc", RIBFLVabcRateLaw, enzConc_RIBFLVabc, 1.66, 0, [0.003, "M_atp_c"], [0.0019, 0.023], ["M_ribflv_c", "M_pi_c", "M_adp_c"], [0.0019, 0.385, 0.654]);

  // ADD FTHFabc
  const enzConc_FTHFabc = partTomM(Math.min(pmap.get("M_PTN_JCVISYN3A_0641_c"), pmap.get("M_PTN_JCVISYN3A_0642_c"), pmap.get("M_PTN_JCVISYN3A_0643_c"), pmap.get("M_PTN_JCVISYN3A_0822_c")), pmap);
  const FTHFabcRateLaw = enzymatic(1, 1);

  addReaction("FTHFabc", FTHFabcRateLaw, enzConc_FTHFabc, 1.66, 0, [0.05, "M_atp_c"], [0.0019, 0.023], ["M_5fthf_c", "M_pi_c", "M_adp_c"], [0.0019, 0.385, 0.654]);

  // ADD THMPPabc
  const enzConc_THMPPabc = partTomM(Math.min(pmap.get("M_PTN_JCVISYN3A_0706_c"), pmap.get("M_PTN_JCVISYN3A_0707_c"), pmap.get("M_PTN_JCVISYN3A_0708_c")), pmap);
  const THMPPabcRateLaw = enzymatic(1, 1);
  model.addMetabolite("M_thmpp_c", "thmpp", partTomM(pmap.get("M_thmpp_c"), pmap));

  addReaction("THMPPabc", THMPPabcRateLaw, enzConc_THMPPabc, 1.66, 0, [0.008, "M_atp_c"], [0.0019, 0.023], ["M_thmpp_c", "M_pi_c", "M_adp_c"], [0.0019, 0.385, 0.654]);

  // ADD P5Pabc
  const enzConc_P5Pabc = partTomM(Math.min(pmap.get("M_PTN_JCVISYN3A_0641_c"), pmap.get("M_PTN_JCVISYN3A_0642_c"), pmap.get("M_PTN_JCVISYN3A_0643_c"), pmap.get("M_PTN_JCVISYN3A_0345_c")), pmap);
  const P5PabcRateLaw = enzymatic(1, 1);
  model.addMetabolite("M_pydx5p_c", "pydx5p", partTomM(pmap.get("M_pydx5p_c"), pmap));

  addReaction("P5Pabc", P5PabcRateLaw, enzConc_P5Pabc, 1.66, 0, [0.006, "M_atp_c"], [0.0019, 0.023], ["M_pydx5p_c", "M_pi_c", "M_adp_c"], [0.0019, 0.385, 0.654]);

  // Ion Transporters

  const Kt6RateLaw = enzymatic(3, 4);
  const enzConc_Kt6 = partTomM(pmap.get("M_PTN_JCVISYN3A_0686_c"), pmap);
  model.addMetabolite("M_k_c", "potassium", partTomM(pmap.get("M_k_c"), pmap));

  addReaction("Kt6", Kt6RateLaw, enzConc_Kt6, 3, 0, [10, "M_atp_c", 10], [0.46, 0.03, 7.47], ["M_k_c", "M_adp_c", "M_pi_c", 134.0], [1.9, 1.5, 10, 12.7]);

  // CA2abc
  const CA2tRateLaw = enzymatic(2, 3);
  const enzConc_cat2 = partTomM(pmap.get("M_PTN_JCVISYN3A_0879_c"), pmap);
  model.addMetabolite("M_ca2_c", "calcium", partTomM(pmap.get("M_ca2_c"), pmap));

  addReaction("CA2abc", CA2tRateLaw, enzConc_cat2, 9.5, 0, [0.68, "M_atp_c"], [0.0075, 0.075], ["M_ca2_c", "M_adp_c", "M_pi_c"], [0.0075, 0.1, 10]);

  // MG2abc
  const MG2tRateLaw = enzymatic(2, 3);
  const enzConc_mg2 = partTomM(pmap.get("M_PTN_JCVISYN3A_0787_c"), pmap);
  model.addMetabolite("M_mg2_c", "magnesium", partTomM(pmap.get("M_mg2_c"), pmap));

  addReaction("MG2abc", MG2tRateLaw, enzConc_mg2, 22, 0, [1, "M_atp_c"], [0.05, 2.8], ["M_mg2_c", "M_adp_c", "M_pi_c"], [0.05, 2.8, 10]);

  // ADD Amino Acid Transporters
  const enzConc_GluTrans = partTomM(pmap.get("M_PTN_JCVISYN3A_0886_c"), pmap);
  const enzConc_AaTransSymp = partTomM(pmap.get("M_PTN_JCVISYN3A_0876_c"), pmap) + partTomM(pmap.get("M_PTN_JCVISYN3A_0878_c"), pmap);
  const enzConc_AaTransAbc = partTomM(pmap.get("M_PTN_JCVISYN3A_0168_c"), pmap); //0165

  const kcat_aa_abc = 0.2; //1/s

  const aaTransport: [string, string, string, number, number, number, number][] = [
    ["R_ARGt2r", "R_ARG4abc", "M_arg__L_c", 0.03, 3.4, kcat_aa_abc, 4],
    ["R_ASPt2pr", "R_ASP4abc", "M_asp__L_c", 0.03, 9.1, kcat_aa_abc, 0.1],
    ["R_GLYt2r", "R_GLY4abc", "M_gly_c", 0.03, 5.1, kcat_aa_abc, 4],
    ["R_ISOt2r", "R_ISO4abc", "M_ile__L_c", 0.03, 9.1, kcat_aa_abc, 4],
    ["R_ALAt2r", "R_ALA4abc", "M_ala__L_c", 0.03, 5.1, kcat_aa_abc, 4],
    ["R_ASNt2r", "R_ASN4abc", "M_asn__L_c", 0.03, 6.1, kcat_aa_abc, 4],
    ["R_LEUt2r", "R_LEU4abc", "M_leu__L_c", 0.03, 9.1, kcat_aa_abc, 4.2],
    ["R_HISt2r", "R_HIS4abc", "M_his__L_c", 0.03, 3.4, kcat_aa_abc, 4],
    ["R_LYSt2r", "R_LYS4abc", "M_lys__L_c", 0.03, 9.1, kcat_aa_abc, 4.2],
    ["R_PROt2r", "R_PRO4abc", "M_pro__L_c", 0.03, 3.4, kcat_aa_abc, 4],
    ["R_PHEt2r", "R_PHE4abc", "M_phe__L_c", 0.03, 3.4, kcat_aa_abc, 1],
    ["R_THRt2r", "R_THR4abc", "M_thr__L_c", 0.03, 5.1, kcat_aa_abc, 1],
    ["R_TRPt2r", "R_TRP4abc", "M_trp__L_c", 0.03, 3.4, kcat_aa_abc, 0.5],
    ["R_TYRt2r", "R_TYR4abc", "M_tyr__L_c", 0.03, 3.4, kcat_aa_abc, 4],
    ["R_VALt2r", "R_VAL4abc", "M_val__L_c", 0.03, 5.1, kcat_aa_abc, 8],
    ["R_SERt2r", "R_SER4abc", "M_ser__L_c", 0.03, 6.1, kcat_aa_abc, 1.1],
    ["R_METt2r", "R_MET4abc", "M_met__L_c", 0.03, 3.4, kcat_aa_abc, 8],
    ["R_CYSt2r", "R_CYS4abc", "M_cys__L_c", 0.03, 3.4, kcat_aa_abc, 2.3],
    ["R_GLUt2pr", "R_GLU4abc", "M_glu__L_c", 0.03, 3.4, kcat_aa_abc, 4.3],
    ["R_GLNt2r", "R_GLN4abc", "M_gln__L_c", 0.03, 7.1, kcat_aa_abc, 2],
  ];

  for (const rxn of aaTransport) {
    const [rxnIDsymp, rxnIDabc, metID, Km, KcatSymp, KcatAbc, ExtConc] = rxn;

    if (!metID.includes("glu") && !metID.includes("gln") && !metID.includes("gly") && !metID.includes("ser")) {
      model.addMetabolite(metID, metID, partTomM(pmap.get(metID), pmap));
    }

    const AaTransSympRateLaw = enzymatic(1, 1);
    addReaction(rxnIDsymp, AaTransSympRateLaw, metID.includes("glu") ? enzConc_GluTrans : enzConc_AaTransSymp, KcatSymp, 0, [ExtConc], [Km], [metID], [Km]);

    const AaTransAbcRateLaw = enzymatic(1, 1);
    addReaction(rxnIDabc, AaTransAbcRateLaw, enzConc_AaTransAbc, KcatAbc, 0, [ExtConc, "M_atp_c"], [Km, 2.8], [metID, "M_adp_c", "M_pi_c"], [Km, 2.8, 10]);
  }

  // Missing Metabolic Reactions
  // GHMT2: h2o_c + methfglu3_c --> 5fthfglu3_c + h_c
  // No exisiting data other than Keq
  model.addRateForm("GHMT2_Rate", enzymatic(1, 1));

  model.addReaction("GHMT2", "GHMT2_Rate", "enzymatic reaction GHMT2");
  model.addParameter("GHMT2", "Enzyme", partTomM(pmap.get("M_PTN_JCVISYN3A_0799_c"), pmap));
  model.addParameter("GHMT2", "kcatF", 640); // Taken from GHMT reaction
  model.addParameter("GHMT2", "kcatR", 23.173); // Taken from GHMT reaction

  model.addSubstrate("GHMT2", "Sub1", "M_methfglu3_c");
  model.addParameter("GHMT2", "KmSub1", 0.0684); // Taken from GHMT reaction, related substrate
  model.addProduct("GHMT2", "Prod1", "M_5fthfglu3_c");
  model.addParameter("GHMT2", "KmProd1", 0.1463); // Taken from GHMT reaction, related substrate
  model.addParameter("GHMT2", "onoff", 1, 0, 1, "mM", "Debug On/Off switch");

  // NADHK: atp_c + nadh_c --> adp_c + h_c + nadph_c
  // exists in the TSV file but not read, code below should not be needed
  model.addRateForm("NADHK_Rate", enzymatic(2, 2));

  model.addReaction("NADHK", "NADHK_Rate", "enzymatic reaction NADHK");
  model.addParameter("NADHK", "Enzyme", partTomM(pmap.get("M_PTN_JCVISYN3A_0259_c"), pmap));
  model.addParameter("NADHK", "kcatF", 31.3405); // From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
  model.addParameter("NADHK", "kcatR", 104.2685); // From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv

  model.addSubstrate("NADHK", "Sub1", "M_nadh_c");
  model.addParameter("NADHK", "KmSub1", 2.0); // 0.1251);//From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
  model.addSubstrate("NADHK", "Sub2", "M_atp_c");
  model.addParameter("NADHK", "KmSub2", 0.4054); // From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv

  model.addProduct("NADHK", "Prod1", "M_nadph_c");
  model.addParameter("NADHK", "KmProd1", 0.3); // 7.4674);//From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
  model.addProduct("NADHK", "Prod2", "M_adp_c");
  model.addParameter("NADHK", "KmProd2", 4.624); // From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
  model.addParameter("NADHK", "onoff", 1, 0, 1, "mM", "Debug On/Off switch");
}