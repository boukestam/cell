import { ParsedGenbank } from "genbank-parser";
import { CME } from "./CME";
import { Avognum, cellVolume, countToMiliMol } from "../Globals";

export function initReplication(sim: CME) {
  const k_high = 7800 * 1000 / Avognum / cellVolume;
  const k_low = 35 * 1000 / Avognum / cellVolume;
  const k_on = 100 * 1000 / Avognum / cellVolume; //molecule^-1 sec^-1
  const k_on_2 = k_on / 2;
  const k_off = 0.55; //sec^-1

  const helicase_removal_rate = 600; //s^-1

  //Chromosome 1
  const nonOricSpec = ['High_Affinity_Site', 'High_Affinity_Bound', 'High_Affinity_Site_oriC',
    'High_Affinity_Bound_oriC', 'Low_Affinity_Site_1', 'Low_Affinity_Site_2',
    'Low_Affinity_Bound_1', 'Low_Affinity_Bound_2', 'ssDNAunboundSite_1'];

  sim.defineSpecies(nonOricSpec);
  sim.addParticles('High_Affinity_Site', 16);
  sim.addParticles('High_Affinity_Bound', 0);
  sim.addReaction(['High_Affinity_Site', 'M_DnaA_c'], ['High_Affinity_Bound'], k_high);

  sim.addParticles('High_Affinity_Site_oriC', 1);
  sim.addParticles('Low_Affinity_Site_1', 0);
  sim.addParticles('Low_Affinity_Site_2', 0);

  sim.addReaction(['High_Affinity_Site_oriC', 'M_DnaA_c'], ['High_Affinity_Bound_oriC', 'Low_Affinity_Site_1'], k_high);
  sim.addReaction(['Low_Affinity_Site_1', 'M_DnaA_c'], ['Low_Affinity_Bound_1', 'Low_Affinity_Site_2'], k_low);
  sim.addReaction(['Low_Affinity_Site_2', 'M_DnaA_c'], ['Low_Affinity_Bound_2', 'ssDNAunboundSite_1'], k_low);

  let species: string[] = [];

  for (let i = 0; i < 30; i++) {         //loop adds 30 terms for unbound sites
    let term = 'ssDNAunboundSite_';
    term = term + (i + 1);
    species.push(term);
  }

  for (let j = 0; j < 30; j++) {         //adds 30 more terms of Bound sites
    let bnd = 'ssDNABoundSite_';
    bnd = bnd + (j + 1);
    species.push(bnd);
  }

  for (let k = 0; k < 30; k++) {
    let unbnd = 'ssDNA_Unbinding_';
    unbnd = unbnd + (k + 1);
    species.push(unbnd);
  }

  species.push('Initiator_C');
  species.push('Initiator_CC');

  sim.defineSpecies(species);

  for (const specy of species) {
    sim.addParticles(specy, 0); //gives all particles count of 0 
  }

  sim.addParticles('ssDNAunboundSite_1', 0); //starts at 1(only site existing @ start

  for (let i = 1; i < 30; i++) {
    sim.addReaction(['ssDNAunboundSite_' + i, 'M_DnaA_c'], ['ssDNAunboundSite_' + (i + 1), 'ssDNABoundSite_' + i], k_on);
    sim.addReaction(['ssDNAunboundSite_' + (i + 1), 'ssDNABoundSite_' + i], ['ssDNAunboundSite_' + i, 'M_DnaA_c'], k_off);
  }

  sim.addReaction(['ssDNAunboundSite_30', 'M_DnaA_c'], ['ssDNABoundSite_30', 'ssDNA_Unbinding_30', 'Initiator_C', 'Initiator_CC'], k_on);

  //Chromosme 2

  const nonOricSpec2 = ['High_Affinity_Site2', 'High_Affinity_Bound2', 'High_Affinity_Site_oriC2',
    'High_Affinity_Bound_oriC2', 'Low_Affinity_Site2_1', 'Low_Affinity_Site2_2',
    'Low_Affinity_Bound2_1', 'Low_Affinity_Bound2_2', 'ssDNAunboundSite2_1'];

  sim.defineSpecies(nonOricSpec2);
  sim.addParticles('High_Affinity_Site2', 0);
  sim.addParticles('High_Affinity_Bound2', 0);
  sim.addReaction(['High_Affinity_Site2', 'M_DnaA_c'], ['High_Affinity_Bound2'], k_high);

  sim.addParticles('High_Affinity_Site_oriC2', 0);
  sim.addParticles('Low_Affinity_Site2_1', 0);
  sim.addParticles('Low_Affinity_Site2_2', 0);

  sim.addReaction(['High_Affinity_Site_oriC2', 'M_DnaA_c'], ['High_Affinity_Bound_oriC2', 'Low_Affinity_Site2_1'], k_high);
  sim.addReaction(['Low_Affinity_Site2_1', 'M_DnaA_c'], ['Low_Affinity_Bound2_1', 'Low_Affinity_Site2_2'], k_low);
  sim.addReaction(['Low_Affinity_Site2_2', 'M_DnaA_c'], ['Low_Affinity_Bound2_2', 'ssDNAunboundSite2_1'], k_low);

  species = [];

  for (let i = 0; i < 30; i++) {         //loop adds 30 terms for unbound sites
    let term = 'ssDNAunboundSite2_';
    term = term + (i + 1);
    species.push(term);
  }

  for (let j = 0; j < 30; j++) {         //adds 30 more terms of Bound sites
    let bnd = 'ssDNABoundSite2_';
    bnd = bnd + (j + 1);
    species.push(bnd);
  }

  for (let k = 0; k < 30; k++) {
    let unbnd = 'ssDNA_Unbinding2_';
    unbnd = unbnd + (k + 1);
    species.push(unbnd);
  }

  species.push('Initiator_C2');
  species.push('Initiator_CC2');

  sim.defineSpecies(species);

  for (const specy of species) {
    sim.addParticles(specy, 0); //gives all particles count of 0 
  }

  sim.addParticles('ssDNAunboundSite2_1', 0); //starts at 1(only site existing @ start

  // Bidning and unbinding reactions are constructed so that DnaA in the middle of the filament cannot unbind and
  // DnaA only bind next to the last DnaA in the filament.

  // Add the binding reactions for each of the 30 binding events in filament formation.
  for (let i = 1; i < 30; i++) {
    sim.addReaction(['ssDNAunboundSite2_' + i, 'M_DnaA_c'], ['ssDNAunboundSite2_' + (i + 1), 'ssDNABoundSite2_' + i], k_on_2);
    sim.addReaction(['ssDNAunboundSite2_' + (i + 1), 'ssDNABoundSite2_' + i], ['ssDNAunboundSite2_' + i, 'M_DnaA_c'], k_off);
  }

  let repInitProds = ['ssDNABoundSite2_30', 'ssDNA_Unbinding2_30', 'Initiator_C', 'Initiator_CC'];

  sim.addReaction(['ssDNAunboundSite2_30', 'M_DnaA_c'], repInitProds, k_on_2);

  //Chromosme 3
  const nonOricSpec3 = ['High_Affinity_Site3', 'High_Affinity_Bound3', 'High_Affinity_Site_oriC3',
    'High_Affinity_Bound_oriC3', 'Low_Affinity_Site3_1', 'Low_Affinity_Site3_2',
    'Low_Affinity_Bound3_1', 'Low_Affinity_Bound3_2', 'ssDNAunboundSite3_1'];
  sim.defineSpecies(nonOricSpec3);
  sim.addParticles('High_Affinity_Site3', 0);
  sim.addParticles('High_Affinity_Bound3', 0);
  sim.addReaction(['High_Affinity_Site3', 'M_DnaA_c'], ['High_Affinity_Bound3'], k_high);

  sim.addParticles('High_Affinity_Site_oriC3', 0);
  sim.addParticles('Low_Affinity_Site3_1', 0);
  sim.addParticles('Low_Affinity_Site3_2', 0);

  sim.addReaction(['High_Affinity_Site_oriC3', 'M_DnaA_c'], ['High_Affinity_Bound_oriC3', 'Low_Affinity_Site3_1'], k_high);
  sim.addReaction(['Low_Affinity_Site3_1', 'M_DnaA_c'], ['Low_Affinity_Bound3_1', 'Low_Affinity_Site3_2'], k_low);
  sim.addReaction(['Low_Affinity_Site3_2', 'M_DnaA_c'], ['Low_Affinity_Bound3_2', 'ssDNAunboundSite3_1'], k_low);

  species = [];

  for (let i = 0; i < 30; i++) {         //loop adds 30 terms for unbound sites
    let term = 'ssDNAunboundSite3_';
    term = term + (i + 1);
    species.push(term);
  }

  for (let j = 0; j < 30; j++) {         //adds 30 more terms of Bound sites
    let bnd = 'ssDNABoundSite3_';
    bnd = bnd + (j + 1);
    species.push(bnd);
  }

  for (let k = 0; k < 30; k++) {
    let unbnd = 'ssDNA_Unbinding3_';
    unbnd = unbnd + (k + 1);
    species.push(unbnd);
  }

  species.push('Initiator_C3');
  species.push('Initiator_CC3');

  sim.defineSpecies(species);

  for (const specy of species) {
    sim.addParticles(specy, 0); //gives all particles count of 0 
  }

  sim.addParticles('ssDNAunboundSite3_1', 0); //starts at 1(only site existing @ start

  // Bidning and unbinding reactions are constructed so that DnaA in the middle of the filament cannot unbind and
  // DnaA only bind next to the last DnaA in the filament.

  // Add the binding reactions for each of the 30 binding events in filament formation.
  for (let i = 1; i < 30; i++) {
    sim.addReaction(['ssDNAunboundSite3_' + i, 'M_DnaA_c'], ['ssDNAunboundSite3_' + (i + 1), 'ssDNABoundSite3_' + i], k_on_2);
    sim.addReaction(['ssDNAunboundSite3_' + (i + 1), 'ssDNABoundSite3_' + i], ['ssDNAunboundSite3_' + i, 'M_DnaA_c'], k_off);
  }

  repInitProds = ['ssDNABoundSite3_30', 'ssDNA_Unbinding3_30', 'Initiator_C', 'Initiator_CC'];

  sim.addReaction(['ssDNAunboundSite3_30', 'M_DnaA_c'], repInitProds, k_on_2);

  // Add the unbinding reactions for each of the 30 possible unbinding events in filament formation.
  for (let i = 2; i < 31; i++) {
    sim.addReaction(['ssDNA_Unbinding_' + i, 'ssDNABoundSite_' + i], ['ssDNA_Unbinding_' + (i - 1), 'M_DnaA_c'], helicase_removal_rate);
  }

  sim.addReaction(['ssDNA_Unbinding_1', 'ssDNABoundSite_1'], ['M_DnaA_c'], helicase_removal_rate);

  for (let i = 2; i < 31; i++) {
    sim.addReaction(['ssDNA_Unbinding2_' + i, 'ssDNABoundSite2_' + i], ['ssDNA_Unbinding2_' + (i - 1), 'M_DnaA_c'], helicase_removal_rate);
  }

  sim.addReaction(['ssDNA_Unbinding2_1', 'ssDNABoundSite2_1'], ['M_DnaA_c'], helicase_removal_rate);

  for (let i = 2; i < 31; i++) {
    sim.addReaction(['ssDNA_Unbinding3_' + i, 'ssDNABoundSite3_' + i], ['ssDNA_Unbinding3_' + (i - 1), 'M_DnaA_c'], helicase_removal_rate);
  }

  sim.addReaction(['ssDNA_Unbinding3_1', 'ssDNABoundSite3_1'], ['M_DnaA_c'], helicase_removal_rate);
}

export function addReplication(sim: CME, dna: ParsedGenbank, ModelSpecies: Set<string>) {
  const K0rep = 0.26e-3;
  const KDrep = 0.001;
  const kcatrep = 100;
  const datp = 0.009 * 2; // mM
  const dttp = 0.011 * 2; // mM
  const dctp = 0.006 * 2; // mM
  const dgtp = 0.0035 * 2; // mM
  const DNApol3 = 35 * countToMiliMol; // mM

  const chromosome_C = ['chromosome_C'];
  sim.defineSpecies(chromosome_C);
  sim.addParticles('chromosome_C', 1);

  const chromosome_CC = ['chromosome_CC'];
  sim.defineSpecies(chromosome_CC);
  sim.addParticles('chromosome_CC', 1);

  const gene_list = [];

  for (let i = 0; i < dna.features.length; i++) {
    if ('product' in dna.features[i].notes) {
      if (dna.features[i].notes['product'][0]) { // Figure out how to sort out for ribosomal operons?
        gene_list.push(i);
      }
    }
  }

  // console.log(gene_list);

  const intergenic_list = [];

  let C_bp = 0;
  let CC_bp = 0;

  let position = 1;

  const CC_genes = [];

  for (const gene of gene_list) {
    const locusTag = dna.features[gene].notes['locus_tag'][0];
    const start = dna.features[gene].start;
    const end = dna.features[gene].end;

    if (start < dna.sequence.length / 2) {
      const geneSeq = dna.sequence.slice(position, end).toUpperCase();

      if (start === 1) {
        const baseCount = new Map<string, number>();
        for (const base of geneSeq) {
          baseCount.set(base, (baseCount.get(base) || 0) + 1);
        }

        const n_tot = [...baseCount.values()].reduce((a, b) => a + b, 0);

        C_bp += n_tot;

        const NMono_A = baseCount.get("A") || 0;
        const NMono_C = baseCount.get("C") || 0;
        const NMono_G = baseCount.get("G") || 0;
        const NMono_T = baseCount.get("T") || 0;

        const NMonoSum = NMono_A * KDrep / datp + NMono_C * KDrep / dctp + NMono_T * KDrep / dttp + NMono_G * KDrep / dgtp;

        const k_gene_dup = kcatrep / ((1 + K0rep / DNApol3) * (KDrep ** 2 / datp / dttp) + NMonoSum + n_tot - 1);

        const intergenic_region = locusTag + '_inter';

        const geneID = locusTag + '_gene';

        if (!ModelSpecies.has(geneID)) {
          ModelSpecies.add(geneID);
          sim.defineSpecies([geneID]);
          sim.addParticles(geneID, 1);
        }

        sim.defineSpecies([intergenic_region]);
        sim.addParticles(intergenic_region, 0);

        const RepProd = [intergenic_region, geneID];

        for (let i = 0; i < geneSeq.length; i++) {
          RepProd.push('ATP_DNArep');
        }

        for (let i = 0; i < NMono_A; i++) {
          RepProd.push('dATP_DNArep');
          RepProd.push('dTTP_DNArep');
        }

        for (let i = 0; i < NMono_T; i++) {
          RepProd.push('dATP_DNArep');
          RepProd.push('dTTP_DNArep');
        }

        for (let i = 0; i < NMono_G; i++) {
          RepProd.push('dGTP_DNArep');
          RepProd.push('dCTP_DNArep');
        }

        for (let i = 0; i < NMono_C; i++) {
          RepProd.push('dCTP_DNArep');
          RepProd.push('dGTP_DNArep');
        }

        sim.addReaction(['chromosome_C', 'Initiator_C'], RepProd, k_gene_dup);

        intergenic_list.push(intergenic_region);

        position = end;
      } else {
        const baseCount = new Map<string, number>();
        for (const base of geneSeq) {
          baseCount.set(base, (baseCount.get(base) || 0) + 1);
        }

        const n_tot = [...baseCount.values()].reduce((a, b) => a + b, 0);

        C_bp += n_tot;

        const NMono_A = baseCount.get("A") || 0;
        const NMono_C = baseCount.get("C") || 0;
        const NMono_G = baseCount.get("G") || 0;
        const NMono_T = baseCount.get("T") || 0;

        const NMonoSum = NMono_A * KDrep / datp + NMono_C * KDrep / dctp + NMono_T * KDrep / dttp + NMono_G * KDrep / dgtp;

        const k_gene_dup = kcatrep / (NMonoSum + n_tot - 1);

        const intergenic_region = locusTag + '_inter';

        const geneID = locusTag + '_gene';

        if (!ModelSpecies.has(geneID)) {
          ModelSpecies.add(geneID);
          sim.defineSpecies([geneID]);
          sim.addParticles(geneID, 1);
        }

        sim.defineSpecies([intergenic_region]);
        sim.addParticles(intergenic_region, 0);

        const RepProd = [intergenic_region, geneID];

        for (let i = 0; i < geneSeq.length; i++) {
          RepProd.push('ATP_DNArep');
        }

        for (let i = 0; i < NMono_A; i++) {
          RepProd.push('dATP_DNArep');
          RepProd.push('dTTP_DNArep');
        }

        for (let i = 0; i < NMono_T; i++) {
          RepProd.push('dATP_DNArep');
          RepProd.push('dTTP_DNArep');
        }

        for (let i = 0; i < NMono_G; i++) {
          RepProd.push('dGTP_DNArep');
          RepProd.push('dCTP_DNArep');
        }

        for (let i = 0; i < NMono_C; i++) {
          RepProd.push('dCTP_DNArep');
          RepProd.push('dGTP_DNArep');
        }

        sim.addReaction(intergenic_list.slice(-1), RepProd, k_gene_dup);

        intergenic_list.push(intergenic_region);

        position = end;
      }
    }
  }

  position = 543086;

  for (const gene of gene_list) {
    if (dna.features[gene].start > dna.sequence.length / 2) {
      CC_genes.push(gene);
    }
  }

  CC_genes.reverse();

  for (const gene of CC_genes) {

    const locusTag = dna.features[gene].notes['locus_tag'][0];
    const start = dna.features[gene].start;
    const end = dna.features[gene].end;

    const geneSeq = dna.sequence.slice(start, position).toUpperCase();

    if (end === 543086) {
      const baseCount = new Map<string, number>();
      for (const base of geneSeq) {
        baseCount.set(base, (baseCount.get(base) || 0) + 1);
      }

      const n_tot = Array.from(baseCount.values()).reduce((a, b) => a + b);

      CC_bp += n_tot;

      const NMono_A = baseCount.get("A") || 0;
      const NMono_C = baseCount.get("C") || 0;
      const NMono_G = baseCount.get("G") || 0;
      const NMono_T = baseCount.get("T") || 0;

      const NMonoSum = NMono_A * KDrep / datp + NMono_C * KDrep / dctp + NMono_T * KDrep / dttp + NMono_G * KDrep / dgtp;

      const k_gene_dup = kcatrep / ((1 + K0rep / DNApol3) * (KDrep ** 2 / datp / dttp) + NMonoSum + n_tot - 1);

      const intergenic_region = locusTag + '_inter';

      const geneID = locusTag + '_gene';

      if (!ModelSpecies.has(geneID)) {
        ModelSpecies.add(geneID);
        sim.defineSpecies([geneID]);
        sim.addParticles(geneID, 1);
      }

      const intergenic = [intergenic_region];
      sim.defineSpecies(intergenic);
      sim.addParticles(intergenic_region, 0);

      const RepProd = [intergenic_region, geneID];

      for (let i = 0; i < geneSeq.length; i++) {
        RepProd.push('ATP_DNArep');
      }

      for (let i = 0; i < NMono_A; i++) {
        RepProd.push('dATP_DNArep');
        RepProd.push('dTTP_DNArep');
      }

      for (let i = 0; i < NMono_T; i++) {
        RepProd.push('dATP_DNArep');
        RepProd.push('dTTP_DNArep');
      }

      for (let i = 0; i < NMono_G; i++) {
        RepProd.push('dGTP_DNArep');
        RepProd.push('dCTP_DNArep');
      }

      sim.addReaction(['chromosome_CC', 'Initiator_CC'], RepProd, k_gene_dup);

      intergenic_list.push(intergenic_region);

      position = start;
    } else if (dna.features[gene].notes['locus_tag'][0] === 'JCVISYN3A_0421') {
      const baseCount = new Map<string, number>();
      for (const base of geneSeq) {
        baseCount.set(base, (baseCount.get(base) || 0) + 1);
      }

      const n_tot = Array.from(baseCount.values()).reduce((a, b) => a + b, 0);

      CC_bp += n_tot;

      const NMono_A = baseCount.get("A") || 0;
      const NMono_C = baseCount.get("C") || 0;
      const NMono_G = baseCount.get("G") || 0;
      const NMono_T = baseCount.get("T") || 0;

      const NMonoSum = NMono_A * KDrep / datp + NMono_C * KDrep / dctp + NMono_T * KDrep / dttp + NMono_G * KDrep / dgtp;

      const k_gene_dup = kcatrep / (NMonoSum + n_tot - 1);

      const intergenic_region = locusTag + '_inter';

      const geneID = locusTag + '_gene';

      if (!ModelSpecies.has(geneID)) {
        ModelSpecies.add(geneID);
        sim.defineSpecies([geneID]);
        sim.addParticles(geneID, 1);
      }

      sim.defineSpecies([intergenic_region]);

      const gene_rep_end_products = ['High_Affinity_Site_oriC2', 'High_Affinity_Site_oriC3',
        'chromosome_C', 'chromosome_CC', 'chromosome_C', 'chromosome_CC',
        intergenic_region, geneID];

      for (let i = 0; i < 16; i++) {
        gene_rep_end_products.push('High_Affinity_Site2');
      }

      for (let i = 0; i < geneSeq.length; i++) {
        gene_rep_end_products.push('ATP_DNArep');
      }

      for (let i = 0; i < NMono_A; i++) {
        gene_rep_end_products.push('dATP_DNArep');
        gene_rep_end_products.push('dTTP_DNArep');
      }

      for (let i = 0; i < NMono_T; i++) {
        gene_rep_end_products.push('dATP_DNArep');
        gene_rep_end_products.push('dTTP_DNArep');
      }

      for (let i = 0; i < NMono_G; i++) {
        gene_rep_end_products.push('dGTP_DNArep');
        gene_rep_end_products.push('dCTP_DNArep');
      }

      for (let i = 0; i < NMono_C; i++) {
        gene_rep_end_products.push('dCTP_DNArep');
        gene_rep_end_products.push('dGTP_DNArep');
      }

      sim.addParticles(intergenic_region, 0);
      sim.addReaction(intergenic_list.slice(-1), gene_rep_end_products, k_gene_dup);

      intergenic_list.push(intergenic_region);

      position = start;
    } else {
      const baseCount = new Map<string, number>();
      for (const base of geneSeq) {
        baseCount.set(base, (baseCount.get(base) || 0) + 1);
      }

      const n_tot = Array.from(baseCount.values()).reduce((a, b) => a + b);

      CC_bp = CC_bp + n_tot;

      const NMono_A = baseCount.get("A") || 0;
      const NMono_C = baseCount.get("C") || 0;
      const NMono_G = baseCount.get("G") || 0;
      const NMono_T = baseCount.get("T") || 0;

      const NMonoSum = NMono_A * KDrep / datp + NMono_C * KDrep / dctp + NMono_T * KDrep / dttp + NMono_G * KDrep / dgtp;

      const k_gene_dup = kcatrep / (NMonoSum + n_tot - 1);

      const intergenic_region = locusTag + '_inter';

      const geneID = locusTag + '_gene';

      if (!ModelSpecies.has(geneID)) {
        ModelSpecies.add(geneID);
        sim.defineSpecies([geneID]);
        sim.addParticles(geneID, 1);
      }

      const intergenic = [intergenic_region];
      sim.defineSpecies(intergenic);
      sim.addParticles(intergenic_region, 0);
      const RepProd = [intergenic_region, geneID];

      for (let i = 0; i < geneSeq.length; i++) {
        RepProd.push('ATP_DNArep');
      }

      for (let i = 0; i < NMono_A; i++) {
        RepProd.push('dATP_DNArep');
        RepProd.push('dTTP_DNArep');
      }

      for (let i = 0; i < NMono_T; i++) {
        RepProd.push('dATP_DNArep');
        RepProd.push('dTTP_DNArep');
      }

      for (let i = 0; i < NMono_G; i++) {
        RepProd.push('dGTP_DNArep');
        RepProd.push('dCTP_DNArep');
      }

      for (let i = 0; i < NMono_C; i++) {
        RepProd.push('dCTP_DNArep');
        RepProd.push('dGTP_DNArep');
      }

      sim.addReaction(intergenic_list.slice(-1), RepProd, k_gene_dup);

      intergenic_list.push(intergenic_region);

      position = start;
    }
  }
}