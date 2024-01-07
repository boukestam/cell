import { ODE, ODEReaction } from "./ODE";
import { lsoda } from "./lsoda";

export type DerivativeFunction = (y: Map<string, number>) => Map<string, number>;

function generateFunction(model: ODE, reactionsByKey: Map<string, ODEReaction[]>): DerivativeFunction {
  let code = '(y) => {\n';

  for (const key of model.metabolites.keys()) {
    code += `const ${key} = y.get("${key}");\n`;
  }

  for (const key of model.explicitParameters.keys()) {
    code += `const ${key} = ${model.explicitParameters.get(key).value};\n`;
  }

  code += 'const result = new Map();\n';
  code += 'const debug = new Map();\n';

  for (const reaction of model.reactions.values()) {
    code += `const Rate_${reaction.id} = ${reaction.compiledRateForm};\n`;
    code += `debug.set("Rate_${reaction.id}", Rate_${reaction.id});\n`;
  }

  for (const key of reactionsByKey.keys()) {
    const reactions = reactionsByKey.get(key);

    code += `result.set("${key}", ${reactions.map(r => `${r.getStoichiometry(key)} * Rate_${r.id}`).join(' + ')});\n`;
  }

  code += 'return result;\n';
  code += '}';

  return eval(code);
}

function euler(
  y: Map<string, number>,
  fn: DerivativeFunction,
  totalTime: number,
  h: number,
  hook: (time: number, d: Map<string, number>) => void
) {
  let time = 0;

  while (time < totalTime) {
    const d = fn(y);

    for (const key of y.keys()) {
      y.set(key, y.get(key) + d.get(key) * h);
    }

    hook(time, d);

    time += h;
  }

  return y;
}

function adams(
  y: Map<string, number>,
  fn: DerivativeFunction,
  totalTime: number,
  h: number,
  hook: (time: number, d: Map<string, number>) => void
): Map<string, number> {
  // Adams Moulton method

  const N = Math.floor(totalTime / h);
  let time = 0;

  for (let n = 0; n < N; n++) {
    // Predictor step
    const dPredict = fn(y);
    const yPredict = new Map<string, number>();
    for (const key of dPredict.keys()) {
      yPredict.set(key, y.get(key) + dPredict.get(key) * h);
    }

    // Corrector step
    let yCorrect = new Map(yPredict);
    let iteration = 0;
    let maxIteration = 100;
    let difference = 0;
    let tolerance = 1 * h;

    do {
      const yPrevious = new Map(yCorrect);
      const dCorrect = fn(yPrevious);

      for (const key of dCorrect.keys()) {
        yCorrect.set(key, y.get(key) + (h / 2) * (dPredict.get(key) + dCorrect.get(key)));
      }

      // Calculate difference
      difference = 0;
      for (const key of yCorrect.keys()) {
        difference += Math.abs(yCorrect.get(key) - yPrevious.get(key));
      }

      iteration++;
    } while (difference > tolerance && iteration < maxIteration);

    // Update y
    for (const key of y.keys()) {
      y.set(key, yCorrect.get(key));
    }

    hook(time, dPredict);

    time += h;
  }

  return y;
}

export function integrate(
  model: ODE,
  totalTime: number,
  h: number,
  method: 'euler' | 'adams' | 'lsoda'
): Map<string, number> {
  const y = new Map<string, number>();
  for (const metabolite of model.metabolites.keys()) {
    y.set(metabolite, model.metabolites.get(metabolite).initialValue);
  }

  const fn = generateFunction(model, model.build());

  const hook = (time: number, d: Map<string, number>) => {
    for (const key of y.keys()) {
      if (isNaN(y.get(key))) {
        console.log(y);
        throw new Error(`NaN for ${key}`);
      }
    }
  };

  if (method === 'euler') {
    return euler(y, fn, totalTime, h, hook);
  } else if (method === 'adams') {
    return adams(y, fn, totalTime, h, hook);
  } else if (method === 'lsoda') {
    return lsoda(y, fn, totalTime, h, hook);
  }
}