export function integrate(
  values: Map<string, number>,
  derivativesByKey: Map<string, string[]>,
  totalTime: number,
  timeStep: number,
  method: 'euler' | 'rk4'
) {
  values = new Map(values);

  let code = '(values) => {\n';

  for (const key of values.keys()) {
    code += `const ${key} = values.get("${key}");\n`;
  }

  code += 'const result = new Map();\n';
  code += 'const debug = new Map();\n';

  for (const key of derivativesByKey.keys()) {
    const derivatives = derivativesByKey.get(key);
    const variables = [];

    for (let i = 0; i < derivatives.length; i++) {
      code += `const d_${key}_${i} = ${derivatives[i]};\n`;
      code += `debug.set("${key}_${i}", d_${key}_${i});\n`;
      variables.push(`d_${key}_${i}`);
    }

    code += `const d_${key} = ${variables.join(' + ')};\n`;
    code += `result.set("${key}", d_${key});\n`;
  }

  code += 'return {result, debug};\n';
  code += '}';

  const fn: ((values: Map<string, number>) => ({ result: Map<string, number>, debug: Map<string, number> })) = eval(code);

  let time = 0;

  const updateValues = (values: Map<string, number>, derivatives: Map<string, number>, h: number) => {
    for (const key of values.keys()) {
      values.set(key, Math.max(0, values.get(key) + derivatives.get(key) * h));
    }
  };

  while (time < totalTime) {
    if (method === "euler") {
      const { result, debug } = fn(values);
      console.log(result, debug);
      throw new Error("stop");
      updateValues(values, result, timeStep);
    } else if (method === "rk4") {
      const k1 = fn(values).result;

      const k2Values = new Map(values);
      updateValues(k2Values, k1, timeStep / 2);
      const k2 = fn(k2Values).result;

      const k3Values = new Map(values);
      updateValues(k3Values, k2, timeStep / 2);
      const k3 = fn(k3Values).result;

      const k4Values = new Map(values);
      updateValues(k4Values, k3, timeStep);
      const k4 = fn(k4Values).result;

      const derivatives = new Map();
      for (const key of values.keys()) {
        derivatives.set(key, (k1.get(key) + 2 * k2.get(key) + 2 * k3.get(key) + k4.get(key)) / 6);
      }

      updateValues(values, derivatives, timeStep);
    }

    for (const key of values.keys()) {
      if (isNaN(values.get(key))) {
        throw new Error(`NaN for ${key}`);
      }
    }

    time += timeStep;
  }

  return values;
}