import { ODE, ODEReaction } from "./ODE";
import { lsoda } from "./lsoda";

export type DerivativeFunction = (y: Float64Array, ydot: Float64Array) => void;

function generateFunction(model: ODE, reactionsByKey: Map<string, ODEReaction[]>): DerivativeFunction {
  let code = '(y, ydot) => {\n';

  for (const key of model.explicitParameters.keys()) {
    code += `const ${key} = ${model.explicitParameters.get(key).value};\n`;
  }

  code += '\n';

  [...model.metabolites.keys()].forEach((key, i) => {
    code += `const ${key} = y[${i}];\n`;
  });

  code += '\n';

  const reactionIdToIndex = new Map<string, number>();

  [...model.reactions.values()].forEach((reaction, i) => {
    reactionIdToIndex.set(reaction.id, i + 1);
    code += `const V${i + 1} = ${reaction.compiledRateForm};\n`;
  });

  code += '\n\n';

  [...reactionsByKey.keys()].forEach((key, i) => {
    const reactions = reactionsByKey.get(key);

    code += `ydot[${i}] = ${reactions.map(r => {
      const stoichiometry = r.getStoichiometry(key);
      let stoichiometryString = stoichiometry.toString();
      if (stoichiometry >= 0) stoichiometryString = '+' + stoichiometryString;

      return `${stoichiometryString}*V${reactionIdToIndex.get(r.id)}`;
    }).join(' ')};\n`;
  });

  code += '}';

  return eval(code);
}

export function integrate(
  model: ODE,
  totalTime: number,
  h: number,
  method: 'lsoda'
): Map<string, number> {
  const y = new Map<string, number>();
  for (const metabolite of model.metabolites.keys()) {
    y.set(metabolite, model.metabolites.get(metabolite).initialValue);
  }

  const fn = generateFunction(model, model.build());

  if (method === 'lsoda') {
    return lsoda(y, fn, totalTime, h);
  }
}