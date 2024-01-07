import { DerivativeFunction } from "./Integrator";

declare global {
  interface Window {
    Module: EmscriptenModule
  }
}

function lsodaCall(
  y: Map<string, number>,
  fn: DerivativeFunction,
  h: number,
  hook: (time: number, d: Map<string, number>) => void
): Map<string, number> {
  const Module: EmscriptenModule = window.Module;
  if (!Module['calledRun']) throw new Error('Module not initialized');

  const yKeys = Array.from(y.keys());
  const yValues = Array.from(y.values());

  const yPointer = Module._malloc(y.size * 8);
  const yArray = Module.HEAPF64.subarray(yPointer / 8, yPointer / 8 + y.size);
  for (let i = 0; i < y.size; i++) {
    yArray[i] = yValues[i];
  }

  const derivePointer = addFunction((t, yPointer, ydotPointer, voidPointer) => {
    // For some reason ydot pointer is offset by 4 bytes, so we need to handle this manually
    const ydot = new Float64Array(y.size);

    // Update y map
    for (let i = 0; i < y.size; i++) {
      y.set(yKeys[i], yArray[i]);
    }

    const derivatives = fn(y);

    for (let i = 0; i < y.size; i++) {
      ydot[i] = derivatives.get(yKeys[i]);
    }

    // Write ydot to memory
    Module.HEAPU8.set(new Uint8Array(ydot.buffer), ydotPointer);

    hook(t, derivatives);

    return 0;
  }, "idiii");

  const lsoda = cwrap("lsoda_run", "number", [
    "number", // int(*derive)(double, double*, double*, void*)
    "number", // double *y
    "number", // int neq
    "number", // double atol
    "number", // double rtol
    "number"  // double time
  ]);

  const status = lsoda(derivePointer, yPointer, y.size, 1e-12, 1e-6, h);

  if (status <= 0) {
    throw new Error(`LSODA failed with status ${status}`);
  }

  const result = new Map<string, number>();

  for (let i = 0; i < y.size; i++) {
    result.set(yKeys[i], yArray[i]);
  }

  removeFunction(derivePointer);

  Module._free(yPointer);

  return result;
}

export function lsoda(
  y: Map<string, number>,
  fn: DerivativeFunction,
  totalTime: number,
  h: number,
  hook: (time: number, d: Map<string, number>) => void
): Map<string, number> {
  const steps = Math.round(totalTime / h);

  for (let i = 0; i < steps; i++) {
    y = lsodaCall(y, fn, h, hook);
  }

  return y;
}