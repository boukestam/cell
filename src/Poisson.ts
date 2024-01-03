// Draw a random number from a Poisson distribution
export function poisson(lam: number = 1.0) {
  const L = Math.exp(-lam);
  let k = 0;
  let p = 1;

  while (p > L) {
    k++;
    p *= Math.random();
  }

  return k - 1;
}