import { CME } from "./CME";
import { data } from "./Data";
import { Avognum, cellVolume } from "./Globals";

function mMtoPart(mM: number) {
  return (mM / 1000) * Avognum * cellVolume;
}

function initNucConcs(sim: CME) {
  const ComDF_nuc = data.ComDF_nuc.get("Compound");

  const nuclMets = ["cytd_c", "cmp_c", "cdp_c", "ctp_c", "dcyt_c", "dcmp_c", "dcdp_c", "dctp_c",
    "gsn_c", "gua_c", "gmp_c", "gdp_c", "gtp_c", "dgsn_c", "dgmp_c", "dgdp_c", "dgtp_c",
    "ura_c", "uri_c", "ump_c", "udp_c", "utp_c", "duri_c", "dump_c", "dudp_c", "dutp_c",
    "adn_c", "ade_c", "dad_2_c", "damp_c", "dadp_c", "datp_c", "thymd_c", "dtmp_c", "dtdp_c", "dttp_c",
    "glu__L_c", "gln__L_c", "ppi_c"]

  for (const met of nuclMets) {
    const metID = "M_" + met;
    const conc = ComDF_nuc.findRow("Compound", metID).InitialConcentration;
    const parts = Math.floor(Math.round(mMtoPart(parseFloat(conc))));
    sim.defineSpecies([metID]);
    sim.addParticles(metID, parts);
  }
}

function initCentsConcs(sim: CME) {
  const ComDF_cent = data.ComDF_cent.get("Compound");

  const centralMets = ["g6p_c", "g1p_c", "f6p_c", "man6p_c", "gam6p_c", "acgam6p_c", "acmanap_c", "pi_c", "amp_c", "adp_c", "atp_c", "acmana_c", "thfglu3_c", "10fthfglu3_c", "nh3_c", "pep_c", "pyr_c", "xu5p__D_c", "ru5p__D_c", "r5p_c", "prpp_c", "e4p_c", "r1p_c", "fdp_c", "s7p_c", "dhap_c", "g3p_c", "2dr5p_c", "2dr1p_c", "nad_c", "nadh_c", "13dpg_c", "nadp_c", "nadph_c", "3pg_c", "2pg_c", "lac__L_c", "co2_c", "accoa_c", "actp_c", "ac_c", "acald_c", "o2_c",
    "nac_c", "nicrnt_c", "dnad_c", "ribflv_c", "fmn_c", "fad_c", "pydx5p_c", "thmpp_c", "5fthf_c", "5fthfglu3_c", "methfglu3_c", "mlthfglu3_c", "sprm_c"];

  for (const met of centralMets) {
    const metID = "M_" + met;
    const conc = ComDF_cent.findRow("Compound", metID).InitialConcentration;
    const parts = Math.floor(Math.round(mMtoPart(parseFloat(conc))));
    sim.defineSpecies([metID]);
    sim.addParticles(metID, parts);
  }
}

function initLipConcs(sim: CME) {
  const ComDF_lip = data.ComDF_lip.get("Compound");

  const lipMets = ["udpg_c", "udpgal_c", "ap_c", "pa_c", "glyc_c", "glyc3p_c", "coa_c", "pap_c", "1ag3p_c", "cdpdag_c", "pg3p_c", "pg_c", "clpn_c", "12dgr_c", "galfur12dgr_c",
    "udpgalfur_c", "fa_c"]; //,"ACP_c","ACP_R_c"] # "fa_c" #"chsterol_c", ACP_R_c, ACP_c

  for (const met of lipMets) {
    const metID = "M_" + met;
    const conc = ComDF_lip.findRow("Compound", metID).InitialConcentration;
    const parts = Math.floor(Math.round(mMtoPart(parseFloat(conc))));
    sim.defineSpecies([metID]);
    sim.addParticles(metID, parts);
  }

  const lipDict_additional = {
    "M_tag_c": 1.95,
    "M_pc_c": 1.20,
    "M_sm_c": 9.22,
    "M_chsterol_c": 23.29,
    'M_atpUsed_GLYK_c': 0.0,
    'M_adpMade_GLYK_c': 0.0,
    'M_adpMade_FAKr_c': 0.0,
    'M_ampMade_BPNT_c': 0.0,
    'M_piMade_BPNT_c': 0.0
  };

  for (const key in lipDict_additional) {
    sim.defineSpecies([key]);
  }

  for (const key in lipDict_additional) {
    const val = Math.floor(Math.round(mMtoPart(lipDict_additional[key])));
    sim.addParticles(key, val);
  }
}

function initAAConcs(sim: CME) {
  const ComDF_aa = data.ComDF_cent.get("Compound");

  const aaMetIDs = ["M_ala__L_c", "M_arg__L_c",
    "M_asn__L_c", "M_asp__L_c", "M_cys__L_c", "M_gly_c",
    "M_his__L_c", "M_ile__L_c", "M_leu__L_c", "M_lys__L_c", "M_met__L_c", "M_phe__L_c",
    "M_pro__L_c", "M_ser__L_c", "M_thr__L_c", "M_trp__L_c", "M_tyr__L_c", "M_val__L_c", "M_amet_c"];

  for (const metID of aaMetIDs) {
    const conc = ComDF_aa.findRow("Compound", metID).InitialConcentration;
    const parts = Math.floor(Math.round(mMtoPart(parseFloat(conc))));
    sim.defineSpecies([metID]);
    sim.addParticles(metID, parts);
  }
}

function initTRNAConcs(sim: CME) {
  sim.defineSpecies(['M_fmettrna_c']);
  sim.addParticles('M_fmettrna_c', 10);

  sim.defineSpecies(['M_glutrnagln_c']);
  sim.addParticles('M_glutrnagln_c', 10);
}

function initTransportConcs(sim: CME) {
  const transport_Dict = {
    "M_ac_e": 0.01, "M_pyr_e": 0.0, "M_lac__L_e": 0.0, "M_h_e": 0.0,
    "M_cytd_e": 0.0, "M_ura_e": 0.0, "M_adn_e": 0.0, "M_dad_2_e": 0.0, "M_gsn_e": 0.0,
    "M_dgsn_e": 0.0, "M_dcyt_e": 0.0, "M_duri_e": 0.0, "M_thymd_e": 0.0, "M_uri_e": 0.0, //"M_co2_c":0.0234,"M_o2_c":0.1,
    "M_glc__D_e": 40, "M_coa_e": 0.0013, "M_glyc_e": 5.0, "M_fa_e": 1.4,
    'M_chsterol_e': 0.005, 'M_nac_e': 0.0016, 'M_arg__L_e': 3.32, 'M_ala__L_e': 2.81,
    'M_asp__L_e': 2.25, 'M_asn__L_e': 0.0, 'M_cys__L_e': 16.5, 'M_glu__L_e': 5.1,
    'M_his__L_e': 0.95, 'M_ile__L_e': 1.53, 'M_leu__L_e': 4.58, 'M_lys__L_e': 3.83,
    'M_met__L_e': 1.01, 'M_phe__L_e': 1.52, 'M_pro__L_e': 3.48, 'M_ser__L_e': 2.38,
    'M_thr__L_e': 2.52, 'M_trp__L_e': 0.49, 'M_tyr__L_e': 2.21, 'M_val__L_e': 2.14,
    'M_gly_e': 6.67
  };

  for (const key in transport_Dict) {
    sim.defineSpecies([key]);
  }

  for (const key in transport_Dict) {
    const val = Math.floor(Math.round(mMtoPart(transport_Dict[key])));
    sim.addParticles(key, val);
  }
}

export function initIC(sim: CME) {
  initLipConcs(sim);
  initNucConcs(sim);
  initCentsConcs(sim);
  initAAConcs(sim);
  initTRNAConcs(sim);
  initTransportConcs(sim);

  sim.defineSpecies(['CellSA']);
  sim.addParticles('CellSA', 231875 + 270956); // Rescale FA, 54% protein

  sim.defineSpecies(['M_k_c']);
  sim.addParticles('M_k_c', mMtoPart(10));

  sim.defineSpecies(['M_ca2_c']);
  sim.addParticles('M_ca2_c', mMtoPart(1.41));

  sim.defineSpecies(['M_mg2_c']);
  sim.addParticles('M_mg2_c', mMtoPart(2.36));

  sim.defineSpecies(['CellSA_Lip']);
  sim.addParticles('CellSA_Lip', 231875); // Rescale FA, 46% Lipid

  // Add a cell surface area covered by protein species, equal to amount covered by prot from SA estimate
  sim.defineSpecies(['CellSA_Prot']);
  sim.addParticles('CellSA_Prot', 270956); // Rescale FA, 54% protein 
  sim.defineSpecies(['CellV']);

  sim.addParticles('CellV', 335); // Cell Volume x 10^19 L
}