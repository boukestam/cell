import { ODE } from "./ODE";

export function pppPatch(model: ODE) {
  // GAPDP Parameters
  model.addParameter('GAPDP', 'KmSub2', 0.385); // nadp
  model.addParameter('GAPDP', 'KmProd2', 0.202); // nadph
  model.addParameter('GAPDP', 'kcatF', 2.8);
  model.addParameter('GAPDP', 'kcatR', 0);

  // FMETTRS Parameters
  model.addParameter('FMETTRS', 'kcatF', 0.45);

  // MTHFC Parameters
  model.addParameter('MTHFC', 'kcatF', 185);

  // GHMT2 Paramters
  model.addParameter('GHMT2', 'kcatF', 0.0);
  model.addParameter('GHMT2', 'kcatR', 0.0);

  // TKT1 Parameters
  model.addParameter('TKT1', 'kcatF', 20.58);
  model.addParameter('TKT1', 'kcatR', 0.8);

  model.addParameter('TKT1', 'KmSub1', 0.743); // g3p
  model.addParameter('TKT1', 'KmSub2', 3.7298); // s7p
  model.addParameter('TKT1', 'KmProd1', 0.4717); // r5p
  model.addParameter('TKT1', 'KmProd2', 0.134); // xu5p

  // TKT2 Parameters
  model.addParameter('TKT2', 'kcatF', 26.87);
  model.addParameter('TKT2', 'kcatR', 1.4);

  model.addParameter('TKT2', 'KmSub1', 0.25); // f6p
  model.addParameter('TKT2', 'KmSub2', 0.743); // g3p
  model.addParameter('TKT2', 'KmProd1', 0.0227); // e4p
  model.addParameter('TKT2', 'KmProd2', 0.134); // xu5p

  // TALA Parameters
  model.addParameter('TALA', 'kcatF', 22.3);
  model.addParameter('TALA', 'kcatR', 0.54);

  model.addParameter('TALA', 'KmSub1', 0.0401); // e4p
  model.addParameter('TALA', 'KmSub2', 0.6688); // f6p
  model.addParameter('TALA', 'KmProd1', 1.9); // g3p
  model.addParameter('TALA', 'KmProd2', 0.285); // s7p

  // Speed up DGSN Pathway
  model.addParameter('DGSNK', 'kcatF', 2.25);

  // Speed up DADN pathway
  model.addParameter('PUNP2', 'kcatF', 13.3);

  // Speed up FBA rxn
  // model.addParameter('FBA',rxnFormKey='kcatF',value=64.5)

  model.addParameter('RNDR2', 'KmSub1', 0.24);

  // RPI Parameters
  model.addParameter('RPI', 'kcatF', 10.0);
  model.addParameter('RPI', 'kcatR', 1.0);

  model.addParameter('FBA', 'KmSub1', 0.12);
  model.addParameter('FBA', 'KmProd2', 0.05);

  model.addParameter('GAPD', 'kcatF', 442.0);
  model.addParameter('GAPD', 'kcatR', 73.6);

  model.addParameter('FBA', 'kcatR', 12.6);

  model.addParameter('TPI', 'kcatR', 67);
  model.addParameter('TPI', 'KmSub1', 0.077);
  model.addParameter('TPI', 'KmProd1', 0.084);

  model.addParameter('FBA', 'kcatF', 21.0);

  model.addParameter('PGK', 'kcatR', 3.4);

  model.addParameter('PGM', 'KmSub1', 3.6);
  model.addParameter('PGM', 'KmProd1', 0.2);

  model.addParameter('PGK', 'KmSub1', 0.01);
  model.addParameter('PGK', 'KmProd1', 0.1);

  model.addParameter('GAPD', 'KmProd1', 0.47);
  model.addParameter('GAPD', 'KmProd2', 0.061);

  model.addParameter('DRPA', 'kcatR', 34.0);

  model.addParameter('DRPA', 'KmProd1', 0.267);
  model.addParameter('DRPA', 'KmProd2', 0.2);

  model.addParameter('PPM2', 'kcatF', 173);

  model.addParameter('PPM2', 'KmSub1', 0.013);
  model.addParameter('PPM2', 'KmProd1', 1.2);
}