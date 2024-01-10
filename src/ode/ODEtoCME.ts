import { ODE } from "./ODE";
import { mMtoPart } from "./Reactions";

export function writeODEtoCME(pmap: Map<string, number>, res: Map<string, number>) {
  // For loop iterating over the species and updating their counts with ODE results
  for (const metabolite of res.keys()) {
    if (metabolite === 'CellSA' || metabolite === 'CellSA_Prot' || metabolite === 'CellSA_Lip') {
      continue;
    } else {
      pmap.set(metabolite, mMtoPart(res.get(metabolite), pmap));
    }
  }

  const ATP_hydro_counters = ['ATP_translat', 'ATP_trsc', 'ATP_mRNAdeg', 'ATP_DNArep', 'ATP_transloc'];

  const NTP_counters = [['ATP_mRNA', 'M_atp_c'], ['CTP_mRNA', 'M_ctp_c'], ['UTP_mRNA', 'M_utp_c'],
  ['GTP_mRNA', 'M_gtp_c'], ['ATP_tRNA', 'M_atp_c'], ['CTP_tRNA', 'M_ctp_c'],
  ['UTP_tRNA', 'M_utp_c'], ['GTP_tRNA', 'M_gtp_c'], ['ATP_rRNA', 'M_atp_c'],
  ['CTP_rRNA', 'M_ctp_c'], ['UTP_rRNA', 'M_utp_c'], ['GTP_rRNA', 'M_gtp_c']];

  const NMP_counters = [['AMP_mRNAdeg', 'M_amp_c'], ['UMP_mRNAdeg', 'M_ump_c'],
  ['CMP_mRNAdeg', 'M_cmp_c'], ['GMP_mRNAdeg', 'M_gmp_c']];

  const dNTP_counters = [['dATP_DNArep', 'M_datp_c'], ['dTTP_DNArep', 'M_dttp_c'],
  ['dCTP_DNArep', 'M_dctp_c'], ['dGTP_DNArep', 'M_dgtp_c']];

  const aatrna_counters = [["FMET_cost", "M_fmettrna_c", "M_trnamet_c"]];

  for (const costID of ATP_hydro_counters) {
    if (costID.includes('translat')) {
      const costCnt = pmap.get(costID);
      const gtpCnt = pmap.get('M_gtp_c');

      if (costCnt > gtpCnt) {
        pmap.set(costID, costCnt - gtpCnt);
        pmap.set('M_gdp_c', pmap.get('M_gdp_c') + gtpCnt);
        pmap.set('M_pi_c', pmap.get('M_pi_c') + gtpCnt);
        pmap.set('M_gtp_c', 0);
      } else {
        pmap.set('M_gtp_c', pmap.get('M_gtp_c') - pmap.get(costID));
        pmap.set('M_gdp_c', pmap.get('M_gdp_c') + pmap.get(costID));
        pmap.set('M_pi_c', pmap.get('M_pi_c') + pmap.get(costID));
        pmap.set(costID, 0);
      }
    } else {
      const costCnt = pmap.get(costID);
      const atpCnt = pmap.get('M_atp_c');

      if (costCnt > atpCnt) {
        pmap.set(costID, costCnt - atpCnt);
        pmap.set('M_adp_c', pmap.get('M_adp_c') + atpCnt);
        pmap.set('M_pi_c', pmap.get('M_pi_c') + atpCnt);
        pmap.set('M_atp_c', 0);
      } else {
        pmap.set('M_atp_c', pmap.get('M_atp_c') - pmap.get(costID));
        pmap.set('M_adp_c', pmap.get('M_adp_c') + pmap.get(costID));
        pmap.set('M_pi_c', pmap.get('M_pi_c') + pmap.get(costID));
        pmap.set(costID, 0);
      }
    }
  }

  for (const cost of NTP_counters) {
    const costID = cost[0];
    const metID = cost[1];

    const cost_count = pmap.get(costID);
    const met_count = pmap.get(metID);

    if (cost_count > met_count) {
      pmap.set(costID, cost_count - met_count);
      pmap.set(metID, 0);
      pmap.set('M_ppi_c', pmap.get('M_ppi_c') + met_count);
    } else {
      pmap.set(metID, pmap.get(metID) - pmap.get(costID));
      pmap.set('M_ppi_c', pmap.get('M_ppi_c') + pmap.get(costID));
      pmap.set(costID, 0);
    }
  }

  for (const cost of dNTP_counters) {
    const costID = cost[0];
    const metID = cost[1];

    const cost_count = pmap.get(costID);
    const met_count = pmap.get(metID);

    if (cost_count > met_count) {
      pmap.set(costID, cost_count - met_count);
      pmap.set(metID, 0);
      pmap.set('M_ppi_c', pmap.get('M_ppi_c') + met_count);
    } else {
      pmap.set(metID, pmap.get(metID) - pmap.get(costID));
      pmap.set('M_ppi_c', pmap.get('M_ppi_c') + pmap.get(costID));
      pmap.set(costID, 0);
    }
  }

  for (const recycle of NMP_counters) {
    const recycleID = recycle[0];
    const metID = recycle[1];

    pmap.set(metID, pmap.get(metID) + pmap.get(recycleID));
    pmap.set(recycleID, 0);
  }

  for (const cost of aatrna_counters) {
    const costID = cost[0];
    const chargedID = cost[1];
    const unchargedID = cost[2];

    const cost_count = pmap.get(costID);
    const charged_pool = pmap.get(chargedID);

    if (cost_count > charged_pool) {
      pmap.set(costID, cost_count - charged_pool);
      pmap.set(unchargedID, pmap.get(unchargedID) + charged_pool);
      pmap.set(chargedID, 0);
    } else {
      pmap.set(chargedID, charged_pool - cost_count);
      pmap.set(unchargedID, pmap.get(unchargedID) + cost_count);
      pmap.set(costID, 0);
    }
  }

  // Call the surface area calculation every write step, every comm. timestep

  const saDict = {
    'M_clpn_c': 0.4,
    'M_chsterol_c': 0.35,
    'M_sm_c': 0.45,
    'M_pc_c': 0.55,
    'M_pg_c': 0.6,
    'M_galfur12dgr_c': 0.6,
    'M_12dgr_c': 0.5,
    'M_pa_c': 0.5,
    'M_cdpdag_c': 0.5,
  }

  let saFlt = 0.0; // surface area float (nm^2)
  for (const key of Object.keys(saDict)) {
    saFlt += pmap.get(key) * saDict[key];
  }

  // Here we partition between inner and outer membrane surface area
  const fracOutMem = 0.513; // 51.3% of the membrane surface is the outer layer of the membrane
  const saFltOut = saFlt * fracOutMem;

  const saInt = Math.round(saFltOut); // Obtain the surface area as an integer nm^2 value for easy storage

  // Calculate Protein Surface Area

  const avgProtSA = 28.0; // nm^2, average protein surface area to produce expected 54% coverage for 9.6K membrane proteins

  // Add code that was used to do this parsing, recalculate Area/protein, remove non-needed ATPase subunits: delta, a, b
  const memProtList = ['JCVISYN3A_0005', 'JCVISYN3A_0008', 'JCVISYN3A_0009', 'JCVISYN3A_0010', 'JCVISYN3A_0011', 'JCVISYN3A_0030', 'JCVISYN3A_0034', 'JCVISYN3A_0060', 'JCVISYN3A_0095', 'JCVISYN3A_0113', 'JCVISYN3A_0114', 'JCVISYN3A_0116', 'JCVISYN3A_0117', 'JCVISYN3A_0132', 'JCVISYN3A_0143', 'JCVISYN3A_0146', 'JCVISYN3A_0164', 'JCVISYN3A_0165', 'JCVISYN3A_0166', 'JCVISYN3A_0167', 'JCVISYN3A_0168', 'JCVISYN3A_0169', 'JCVISYN3A_0195', 'JCVISYN3A_0196', 'JCVISYN3A_0197', 'JCVISYN3A_0235', 'JCVISYN3A_0239', 'JCVISYN3A_0248', 'JCVISYN3A_0249', 'JCVISYN3A_0296', 'JCVISYN3A_0304', 'JCVISYN3A_0314', 'JCVISYN3A_0317', 'JCVISYN3A_0326', 'JCVISYN3A_0332', 'JCVISYN3A_0338', 'JCVISYN3A_0345', 'JCVISYN3A_0346', 'JCVISYN3A_0371', 'JCVISYN3A_0372', 'JCVISYN3A_0379', 'JCVISYN3A_0388', 'JCVISYN3A_0398', 'JCVISYN3A_0399', 'JCVISYN3A_0411', 'JCVISYN3A_0425', 'JCVISYN3A_0426', 'JCVISYN3A_0427', 'JCVISYN3A_0428', 'JCVISYN3A_0439', 'JCVISYN3A_0440', 'JCVISYN3A_0478', 'JCVISYN3A_0481', 'JCVISYN3A_0505', 'JCVISYN3A_0516', 'JCVISYN3A_0601', 'JCVISYN3A_0639', 'JCVISYN3A_0641', 'JCVISYN3A_0642', 'JCVISYN3A_0643', 'JCVISYN3A_0652', 'JCVISYN3A_0685', 'JCVISYN3A_0686', 'JCVISYN3A_0691', 'JCVISYN3A_0696', 'JCVISYN3A_0706', 'JCVISYN3A_0707', 'JCVISYN3A_0708', 'JCVISYN3A_0774', 'JCVISYN3A_0777', 'JCVISYN3A_0778', 'JCVISYN3A_0779', 'JCVISYN3A_0787', 'JCVISYN3A_0789', 'JCVISYN3A_0790', 'JCVISYN3A_0791', 'JCVISYN3A_0792', 'JCVISYN3A_0795', 'JCVISYN3A_0797', 'JCVISYN3A_0822', 'JCVISYN3A_0827', 'JCVISYN3A_0830', 'JCVISYN3A_0835', 'JCVISYN3A_0836', 'JCVISYN3A_0839', 'JCVISYN3A_0852', 'JCVISYN3A_0870', 'JCVISYN3A_0872', 'JCVISYN3A_0876', 'JCVISYN3A_0878', 'JCVISYN3A_0879', 'JCVISYN3A_0881', 'JCVISYN3A_0908'] // _0793, _0794, _0796 = ATPase not in transmembrane

  // Proteins with special naming conventions in the model
  const otherNamesDict = { 'JCVISYN3A_0779': ['ptsg', 'ptsg_P'] };

  let count = 0; // Count of number of membrane proteins
  for (const elementID of memProtList) {
    if (elementID in otherNamesDict) {
      const elementIDList = otherNamesDict[elementID];
      for (const obj of elementIDList) {
        count += pmap.get(obj);
      }
    } else {
      count += pmap.get('M_PTN_' + elementID + '_c');
    }
  }

  const protSA = Math.round(count * avgProtSA); // nm^2

  // Assign the values
  pmap.set('CellSA_Lip', saInt);
  pmap.set('CellSA', saInt + protSA);
  pmap.set('CellSA_Prot', protSA);
}
