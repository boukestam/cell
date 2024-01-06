import { Avognum } from "../Globals";

export function calcCellVolume(pmap: Map<string, number>): number {
  const surfaceArea = pmap.get('CellSA');
  const cellRadius = Math.sqrt(surfaceArea / (4 * Math.PI)) * 1e-9;
  const cellVolume = (4 / 3) * Math.PI * Math.pow(cellRadius, 3) * 1000;
  if (cellVolume > 6.70e-17) {
    pmap.set("CellV", 670);
    return 6.70e-17;
  } else {
    pmap.set("CellV", Math.round(cellVolume * 1e19));
    return cellVolume;
  }
}

export function partTomM(particles: number, pmap: Map<string, number>): number {
  const cellVolume = calcCellVolume(pmap);
  const concentration = (particles * 1000) / (Avognum * cellVolume);

  if (isNaN(concentration)) {
    console.error('NaN concentration', particles, cellVolume);
    throw new Error('NaN concentration');
  }

  return concentration;
}

export function mMtoPart(conc: number, pmap: Map<string, number>): number {
  const cellVolume = calcCellVolume(pmap);
  const particle = Math.floor(Math.round((conc / 1000) * Avognum * cellVolume));

  if (isNaN(particle)) {
    console.error('NaN particle', conc, cellVolume);
    throw new Error('NaN particle');
  }

  return particle;
}

export function enzymatic(subs: number, prods: number): string {
  const numerator = (subs: number, prods: number) => {
    const subterm = [...Array(subs).keys()].map(i => `( $Sub${i + 1} / $KmSub${i + 1} )`);
    const subNumer = subterm.join(' * ');

    const prodterm = [...Array(prods).keys()].map(i => `( $Prod${i + 1} / $KmProd${i + 1} )`);
    const prodNumer = prodterm.join(' * ');

    const numerator = `( $kcatF * ${subNumer} - $kcatR * ${prodNumer} )`;
    return numerator;
  }

  const denominator = (subs: number, prods: number) => {
    const subterm = [...Array(subs).keys()].map(i => `( 1 + $Sub${i + 1} / $KmSub${i + 1} )`);
    const subDenom = subterm.join(' * ');

    const prodterm = [...Array(prods).keys()].map(i => `( 1 + $Prod${i + 1} / $KmProd${i + 1} )`);
    const prodDenom = prodterm.join(' * ');

    const denominator = `( ${subDenom} + ${prodDenom} - 1 )`;
    return denominator;
  }

  const rate = `$onoff * $Enzyme * ( ${numerator(subs, prods)} / ${denominator(subs, prods)} )`;

  return rate;
}
