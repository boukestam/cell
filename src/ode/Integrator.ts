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


function generatePythonCode(model: ODE, reactionsByKey: Map<string, ODEReaction[]>) {
  let code = 'def compiledCall(self, t,y):\n\n';

  const ident = '    ';

  for (const key of model.explicitParameters.keys()) {
    code += ident + `# ${model.explicitParameters.get(key).name}\n`;
    code += ident + `${key} = ${model.explicitParameters.get(key).value}\n\n`;
  }

  code += '\n\n';

  [...model.metabolites.keys()].forEach((key, i) => {
    code += ident + `# (${key}) ${model.metabolites.get(key).name}\n`;
    code += ident + `${key} = y[${i}]\n\n`;
  });

  code += '\n\n';

  const reactionIdToIndex = new Map<string, number>();

  [...model.reactions.values()].forEach((reaction, i) => {
    reactionIdToIndex.set(reaction.id, i + 1);
    code += ident + `# (${reaction.id}) ${reaction.name}\n`;
    code += ident + `V${i + 1} = ${reaction.compiledRateForm}\n\n`;
  });

  code += '\n\n';

  for (const key of reactionsByKey.keys()) {
    const reactions = reactionsByKey.get(key);

    // Split reactions into groups of 10
    const reactionGroups = [];
    for (let i = 0; i < reactions.length; i += 10) {
      reactionGroups.push(reactions.slice(i, i + 10));
    }

    for (let i = 0; i < reactionGroups.length; i++) {
      const reactions = reactionGroups[i];

      code += ident + `Delta${key} ${i === 0 ? '=' : '+='} ${reactions.map(r => {
        const stoichiometry = r.getStoichiometry(key);
        let stoichiometryString = stoichiometry.toString();
        if (stoichiometry >= 0) stoichiometryString = '+' + stoichiometryString;

        return `${stoichiometryString}.0*V${reactionIdToIndex.get(r.id)}`;
      }).join(' ')}${i === reactionGroups.length - 1 ? ' + 0' : ''}\n`;
    }
  }

  code += '\n';
  code += ident + `return [${[...reactionsByKey.keys()].map(key => 'Delta' + key).join(', ')}]\n`;
  code += '\n';

  return code;
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

  // const yArray = [...y.values()];

  // const compareArray = `3.40230918e+00 1.18630544e+00 2.74180498e+00 1.06056537e+00
  // 7.56497180e+00 6.94942204e-03 6.29177119e-03 2.19367188e+00
  // 1.74347057e+01 6.61030903e-01 1.48575880e-02 1.21272592e-02
  // 1.26963090e+00 3.21009641e-02 1.07182936e-01 4.07228638e-02
  // 2.44481399e+00 1.04157433e+00 1.83467811e-01 3.76532901e-01
  // 7.11082532e-03 8.28898731e-01 1.90330343e-03 4.57405971e-02
  // 1.05382232e-01 1.82884484e-02 4.59880432e-01 7.20673610e-02
  // 7.67747588e-02 7.43514834e-02 1.66592507e-01 6.19934338e-02
  // 1.54231795e-01 2.02331917e-04 1.06962355e-01 3.52178763e-03
  // 7.66123320e-01 3.97781499e-03 3.21810567e-03 1.53645748e-02
  // 1.78313708e-02 5.32386833e-03 1.77536419e-02 2.10632188e-02
  // 2.11786046e-02 5.81598116e-01 8.33273778e-03 6.47015851e-01
  // 5.16453348e-03 8.37336045e-01 4.51467579e-03 7.97434633e-03
  // 8.22026353e-03 1.60831671e-01 6.76062088e-01 6.78996760e-01
  // 6.94142111e-01 2.17301987e-03 3.37348011e-01 3.34232386e-01
  // 3.48215992e-01 1.19390497e-02 9.94362810e-02 1.20705774e-02
  // 1.24879009e-02 8.95480323e-02 1.28738745e-02 8.56864885e-02
  // 3.57199149e-02 2.25860985e-02 2.23311096e-02 2.13775219e-03
  // 2.61758984e-01 3.27886469e+00 1.00204265e-01 1.45510155e-03
  // 1.95077730e+00 9.55983965e-04 5.71067634e-03 6.98917054e-03
  // 6.07269914e+00 3.23411149e-02 6.91721692e-05 1.02027109e-01
  // 1.00523582e-01 9.95639805e-02 9.91230463e-02 1.22001320e+00
  // 5.20754386e+00 6.98630177e+00 1.09087913e-01 2.15339273e+00
  // 6.56481523e-02 1.01693090e-01 5.14735675e+00 1.09549795e-01
  // 2.36022065e-02 6.66601634e-04 9.52538058e-02 1.70087457e-02
  // 4.25173640e-04 1.06802280e-02 3.63341596e-03 1.13297698e-02
  // 4.17255048e-03 2.85363974e-02 1.25719920e-02 2.32780403e+01
  // 1.00086357e-01 9.92902951e-02 1.00616334e-01 1.97945224e-02
  // 5.42102888e-02 6.39906438e-02 9.99738155e-02 9.95141456e-02
  // 5.48503031e-02 3.68655336e-01 7.20801302e-01 1.04614949e-01
  // 1.62648942e+00 9.99606997e-02 9.99538546e-02 9.99524411e-02
  // 9.99703611e+00 1.41060682e+00 2.35936949e+00 4.93008626e-01
  // 2.12968542e+00 1.18214072e-01 1.49041136e+00 1.76225076e-01
  // 1.82813750e-01 6.95021099e-02 2.56963473e-01 3.71297271e-01
  // 2.86922896e-01 3.42503768e-01 3.28478608e-02 3.13267263e-02
  // 1.81877956e+00 8.30354428e-02 3.51258128e-02 9.21527276e+00
  // 1.19938849e+00 1.94901411e+00`.split(/\s+/).map(v => parseFloat(v));

  // const percentageDifferences = yArray.map((v, i) => Math.abs(v - compareArray[i]) / compareArray[i] * 100);

  // const yStringArray = percentageDifferences.map(v => v.toExponential(8));
  // const keys = [...y.keys()];
  // // Split into groups of 4
  // const yArrayGroups = [];
  // const keyGroups = [];
  // for (let i = 0; i < yStringArray.length; i += 4) {
  //   yArrayGroups.push(yStringArray.slice(i, i + 4));
  //   keyGroups.push(keys.slice(i, i + 4));
  // }

  // let out = yArrayGroups.map(array => array.join(" ")).join('\n');
  // // replace e[+-]n with e[+-]0n
  // out = out.replace(/e([+-])(\d)/g, 'e$10$2');

  // console.log(out);

  // console.log(keyGroups.map(array => array.join(" ")).join('\n'));

  // console.log(generatePythonCode(model, model.build()));

  // Replace y with compareArray
  // for (let i = 0; i < yArray.length; i++) {
  //   y.set(y.keys()[i], compareArray[i]);
  // }

  const fn = generateFunction(model, model.build());

  const hook = (time: number, d: Map<string, number>) => {
    for (const key of y.keys()) {
      if (isNaN(y.get(key))) {
        console.log(y);
        throw new Error(`NaN for ${key}`);
      }
    }
  };

  try {
    if (method === 'euler') {
      return euler(y, fn, totalTime, h, hook);
    } else if (method === 'adams') {
      return adams(y, fn, totalTime, h, hook);
    } else if (method === 'lsoda') {
      return lsoda(y, fn, totalTime, h, hook);
    }
  } catch (e) {
    console.error(e);
    return y;
  }
}