import * as math from "mathjs";
import { Datasheet } from "./Datasheet";

// expression + '_' + metaboliteId => derivative
const derivativeCache = new Map<string, math.MathNode>();

export class ODERateForm {
  template: string;
  keySet: Set<string>;

  constructor(template: string = "$Vmax * $Sub1 /($Km + $Sub1)") {
    this.template = template;

    this.keySet = new Set<string>();
    for (const key of template.match(/\$[a-zA-Z0-9_]+/g)) {
      this.keySet.add(key.replace("$", ""));
    }
  }

  compile(parameters: Map<string, ODEParameterValue>) {
    let rateForm = this.template;
    for (const key of this.keySet) {
      if (!parameters.has(key)) {
        throw new Error(`Parameter ${key} not found in rate form`);
      }

      rateForm = rateForm.replaceAll(`$${key}`, parameters.get(key).toString());
    }

    return rateForm;
  }
}

export class ODEReaction {
  id: string;
  name: string;
  substrates: Map<string, string>;
  products: Map<string, string>;
  parameters: Map<string, ODEParameter>;
  stochiometry: Map<string, number>;

  rateForm: ODERateForm;

  constructor(id: string, name: string) {
    this.id = id;
    this.name = name;
    this.substrates = new Map();
    this.products = new Map();
    this.parameters = new Map();
    this.stochiometry = new Map();
  }

  setRateForm(rateForm: ODERateForm) {
    this.rateForm = rateForm;
  }

  getKeys() {
    if (this.rateForm === undefined) {
      throw new Error("Rate form not set on reaction: " + this.id);
    }

    return this.rateForm.keySet;
  }

  addSubstrate(formKey: string, metaboliteId: string, stoichiometry: number = 0) {
    this.substrates.set(formKey, metaboliteId);
    if (stoichiometry !== 0) {
      this.stochiometry.set(metaboliteId, stoichiometry);
    }
  }

  addProduct(formKey: string, metaboliteId: string, stoichiometry: number = 0) {
    this.products.set(formKey, metaboliteId);
    if (stoichiometry !== 0) {
      this.stochiometry.set(metaboliteId, stoichiometry);
    }
  }

  getStoichiometry(metaboliteId: string) {
    if (this.stochiometry.has(metaboliteId)) {
      return this.stochiometry.get(metaboliteId)
    } else {
      const subStoich = Array.from(this.substrates.values()).filter(substrate => substrate === metaboliteId).length;
      const prodStoich = Array.from(this.products.values()).filter(product => product === metaboliteId).length;
      return prodStoich - subStoich;
    }
  }
}

export class ODEMetabolite {
  id: string;
  name: string;
  initialValue: number;

  constructor(id: string, name: string, initialValue: number) {
    this.id = id;
    this.name = name;
    this.initialValue = initialValue;
  }
}

export type ODEParameterValue = number | string;

export class ODEParameter {
  name: string;
  value: ODEParameterValue;
  unit: string;
  formKey: string;
  lowerBound: number;
  upperBound: number;

  constructor(name: string, value: ODEParameterValue, unit: string, formKey: string, lowerBound: number, upperBound: number) {
    this.name = name;
    this.value = value;
    this.unit = unit;
    this.formKey = formKey;
    this.lowerBound = lowerBound;
    this.upperBound = upperBound;
  }
}

export class ODE {
  rateForms: Map<string, ODERateForm>;
  reactions: Map<string, ODEReaction>;
  metabolites: Map<string, ODEMetabolite>;

  explicitParameters: Map<string, ODEParameter>;
  globalParameters: Map<string, ODEParameter>;
  reactionParameters: Map<string, ODEParameter>;

  constructor() {
    this.rateForms = new Map();
    this.reactions = new Map();
    this.metabolites = new Map();

    this.explicitParameters = new Map();
    this.globalParameters = new Map();
    this.reactionParameters = new Map();
  }

  addRateForm(id: string, template: string) {
    const rateForm = new ODERateForm(template);
    this.rateForms.set(id, rateForm);
  }

  getReaction(id: string) {
    return this.reactions.get(id);
  }

  addReaction(id: string, rateFormId: string, name: string) {
    if (this.reactions.has(id)) {
      throw new Error(`Reaction ${id} already exists`);
    }

    if (!this.rateForms.has(rateFormId)) {
      throw new Error(`Rate form ${rateFormId} not found`);
    }

    const reaction = new ODEReaction(id, name);
    reaction.setRateForm(this.rateForms.get(rateFormId));
    this.reactions.set(id, reaction);
  }

  addSubstrate(reactionId: string, reactionFormKey, metaboliteId: string, stoichiometry: number = 0) {
    if (!this.reactions.has(reactionId)) {
      throw new Error(`Reaction ${reactionId} not found`);
    }

    if (!this.metabolites.has(metaboliteId)) {
      throw new Error(`Metabolite ${metaboliteId} not found`);
    }

    this.reactions.get(reactionId).addSubstrate(reactionFormKey, metaboliteId, stoichiometry);
  }

  addProduct(reactionId: string, reactionFormKey, metaboliteId: string, stoichiometry: number = 0) {
    if (!this.reactions.has(reactionId)) {
      throw new Error(`Reaction ${reactionId} not found`);
    }

    if (!this.metabolites.has(metaboliteId)) {
      throw new Error(`Metabolite ${metaboliteId} not found`);
    }

    this.reactions.get(reactionId).addProduct(reactionFormKey, metaboliteId, stoichiometry);
  }

  addMetabolite(id: string, name: string, initialValue: number) {
    const metabolite = new ODEMetabolite(id, name, initialValue);
    this.metabolites.set(id, metabolite);
  }

  addParameter(reactionId: string, formKey: string, value: ODEParameterValue, lowerBound: number = 0, upperBound: number = 0, unit: string = "", name: string = "") {
    const parameter = new ODEParameter(name, value, unit, formKey, lowerBound, upperBound);

    if (reactionId === 'Explicit') {
      this.explicitParameters.set(formKey, parameter);
    } else if (reactionId === 'Global') {
      this.globalParameters.set(formKey, parameter);
    } else {
      this.reactionParameters.set(formKey, parameter);
      this.reactions.get(reactionId).parameters.set(formKey, parameter);
    }
  }

  parseParameters(parameters: Datasheet) {
    for (const row of parameters.getRows()) {
      const rxnID = row["ReactionID"];
      const paramKey = row["FormKey"];
      const parVal = row["parVal"];
      const parName = row["parName"];
      const parUnit = row["unit"];

      let parLB = 0;
      let parUB = 0;

      if ("optimizeMin" in row && "optimizeMax" in row) {
        if (row["optimizeMin"] == "") {
          parLB = 0;
        } else {
          parLB = row["optimizeMin"];
        }

        if (row["optimizeMax"] == "") {
          parUB = 0;
        } else {
          parUB = row["optimizeMax"];
        }
      }

      this.addParameter(rxnID, paramKey, parVal, parLB, parUB, parUnit, parName);
    }
  }

  build() {
    // Check that all metabolites in reactions are defined
    for (const reaction of this.reactions.values()) {
      for (const substrate of reaction.substrates.values()) {
        if (this.metabolites.has(substrate)) continue;
        if (this.explicitParameters.has(substrate)) continue;

        throw new Error(`Metabolite ${substrate} not found in reaction ${reaction.id}`);
      }

      for (const product of reaction.products.values()) {
        if (this.metabolites.has(product)) continue;
        if (this.explicitParameters.has(product)) continue;

        throw new Error(`Metabolite ${product} not found in reaction ${reaction.id}`);
      }
    }

    const parseValue = (value: ODEParameterValue): ODEParameterValue => {
      if (typeof value === "string") {
        if (this.metabolites.has(value)) {
          return value;
        } else if (this.explicitParameters.has(value)) {
          return this.explicitParameters.get(value).value as number;
        } else {
          throw new Error(`Parameter not found in reaction ${value}`);
        }
      }

      return value;
    };

    const metaboliteDerivatives = new Map<string, string[]>();
    for (const metabolite of this.metabolites.values()) {
      metaboliteDerivatives.set(metabolite.id, []);
    }

    // Check that all reaction parameters are defined
    for (const reaction of this.reactions.values()) {
      if (reaction.rateForm === undefined) {
        throw new Error(`Rate form not set for reaction ${reaction.id}`);
      }

      const values = new Map<string, ODEParameterValue>();

      for (const parameter of reaction.rateForm.keySet) {
        try {
          if (reaction.parameters.has(parameter)) {
            values.set(parameter, parseValue(reaction.parameters.get(parameter).value));
            continue;
          }
          if (reaction.substrates.has(parameter)) {
            values.set(parameter, parseValue(reaction.substrates.get(parameter)));
            continue;
          }
          if (reaction.products.has(parameter)) {
            values.set(parameter, parseValue(reaction.products.get(parameter)));
            continue;
          }
          if (this.explicitParameters.has(parameter)) {
            values.set(parameter, parseValue(this.explicitParameters.get(parameter).value));
            continue;
          }
        } catch (e) {
          console.error(e);
        }

        console.log(reaction);
        throw new Error(`Parameter ${parameter} not found in reaction ${reaction.id}`);
      }

      const REPLACE_STRING = "xzx";
      const compiled = reaction.rateForm.compile(values).replaceAll("_", REPLACE_STRING);

      const root = math.parse(compiled);

      let xi = 0;
      const xToReal = new Map<string, string | number>();
      const realToX = new Map<string | number, string>();

      const mapToX = (node: math.MathNode): math.MathNode => {
        if (node.type === "ConstantNode" || node.type === "SymbolNode") {
          if (node.type === "SymbolNode" && realToX.has((node as math.SymbolNode).name)) {
            return new math.SymbolNode(realToX.get((node as math.SymbolNode).name));
          }

          const xName = "x" + (xi++);

          if (node.type === "ConstantNode") {
            xToReal.set(xName, (node as math.ConstantNode).value);
          }

          if (node.type === "SymbolNode") {
            xToReal.set(xName, (node as math.SymbolNode).name);
            realToX.set((node as math.SymbolNode).name, xName);
          }

          return new math.SymbolNode(xName);
        } else {
          return node;
        }
      };

      const mapToReal = (node: math.MathNode): math.MathNode => {
        if (node.type === "SymbolNode") {
          const name = (node as math.SymbolNode).name;
          if (!xToReal.has(name)) throw new Error(`Symbol ${name} not found in xToReal`);

          const value = xToReal.get(name);
          if (typeof value === "number") {
            return new math.ConstantNode(value);
          } else {
            return new math.SymbolNode(value);
          }
        } else {
          return node;
        }
      };

      const xRoot = root.transform(mapToX);

      const getDerivative = (metaboliteId: string): math.MathNode => {
        const replacedMetaboliteId = metaboliteId.replaceAll("_", REPLACE_STRING);

        if (!realToX.has(replacedMetaboliteId)) {
          return xRoot;
        }

        const xRootString = xRoot.toString();
        const xName = realToX.get(replacedMetaboliteId);
        const key = xRootString + "_" + xName;

        if (derivativeCache.has(key)) return derivativeCache.get(key);

        const derivative = math.derivative(xRoot, xName);
        derivativeCache.set(key, derivative);

        return derivative;
      };

      for (const metabolite of [...reaction.substrates.values(), ...reaction.products.values()]) {
        const d = getDerivative(metabolite).transform(mapToReal).toString().replaceAll(REPLACE_STRING, "_");
        metaboliteDerivatives.get(metabolite).push(`${reaction.getStoichiometry(metabolite)} * (${d})`);
      }
    }

    return metaboliteDerivatives;
  }
}