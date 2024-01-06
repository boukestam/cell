import { XMLParser } from "fast-xml-parser"

export interface SBMLUnitDefinition {
  id: string,
  units: {
    kind: string,
    exponent: number,
    scale: number,
    multiplier: number
  }[]
}

export interface SBMLCompartment {
  id: string,
  constant: boolean,
}

export interface SBMLSpecies {
  id: string,
  name: string,
  compartment: string,
  hasOnlySubstanceUnits?: boolean,
  boundaryCondition?: boolean,
  constant: boolean,
  initialAmount?: number,
  substanceUnits?: string,
  'fbc:charge'?: number,
  'fbc:chemicalFormula'?: string,
}

export interface SBMLReaction {
  id: string,
  name: string,
  reversible: boolean,
  'fbc:lowerFluxBound': number,
  'fbc:upperFluxBound': number,
  listOfReactants: {
    species: string,
    stoichiometry: number,
    constant: boolean
  }[],
  listOfProducts: {
    species: string,
    stoichiometry: number,
    constant: boolean
  }[]
}

export interface SBMLParameter {
  id: string,
  value: number,
  units: string,
  constant: boolean,
  sboTerm: string,
}

export interface FBCObjective {
  id: string,
  type: string,
  listOfFluxObjectives: {
    reaction: string,
    coefficient: number,
  }[]
}

export interface FBCGeneProduct {
  id: string,
  name: string,
  label: string
}

export class SBMLModel {
  id: string;
  'fbc:strict': boolean;
  unitDefinitions: SBMLUnitDefinition[];
  compartments: SBMLCompartment[];
  species: SBMLSpecies[];
  reactions: SBMLReaction[];
  parameters: SBMLParameter[];
  'fbc:objectives': FBCObjective[];
  'fbc:geneProducts': FBCGeneProduct[];

  constructor(data: any) {
    this.id = data.id;
    this['fbc:strict'] = data['fbc:strict'] === 'true';

    this.unitDefinitions = data.listOfUnitDefinitions.unitDefinition.map(unitDefinition => ({
      id: unitDefinition.id,
      units: unitDefinition.listOfUnits.unit.map(unit => ({
        kind: unit.kind,
        exponent: parseFloat(unit.exponent),
        scale: parseFloat(unit.scale),
        multiplier: parseFloat(unit.multiplier)
      }))
    }));

    this.parameters = data.listOfParameters.parameter.map(parameter => ({
      id: parameter.id,
      value: parseFloat(parameter.value),
      units: parameter.units,
      constant: parameter.constant === 'true',
      sboTerm: parameter.sboTerm
    }));

    this.compartments = data.listOfCompartments.compartment.map(compartment => ({
      id: compartment.id,
      constant: compartment.constant === 'true'
    }));

    this.species = data.listOfSpecies.species.map(species => ({
      id: species.id,
      name: species.name,
      compartment: species.compartment,
      hasOnlySubstanceUnits: species.hasOnlySubstanceUnits === 'true',
      boundaryCondition: species.boundaryCondition === 'true',
      constant: species.constant === 'true',
      'fbc:charge': this.parseValue(species['fbc:charge']),
      'fbc:chemicalFormula': species['fbc:chemicalFormula']
    }));

    this.reactions = data.listOfReactions.reaction.map(reaction => ({
      id: reaction.id,
      name: reaction.name,
      reversible: reaction.reversible === 'true',
      'fbc:lowerFluxBound': this.parseValue(reaction['fbc:lowerFluxBound']),
      'fbc:upperFluxBound': this.parseValue(reaction['fbc:upperFluxBound']),
      listOfReactants: reaction.listOfReactants?.speciesReference.map(speciesReference => ({
        species: speciesReference.species,
        stoichiometry: this.parseValue(speciesReference.stoichiometry),
        constant: speciesReference.constant === 'true'
      })) || [],
      listOfProducts: reaction.listOfProducts?.speciesReference.map(speciesReference => ({
        species: speciesReference.species,
        stoichiometry: this.parseValue(speciesReference.stoichiometry),
        constant: speciesReference.constant === 'true'
      })) || []
    }));

    this['fbc:objectives'] = data['fbc:listOfObjectives']['fbc:objective'].map(objective => ({
      id: objective["fbc:id"],
      type: objective["fbc:type"],
      listOfFluxObjectives: objective["fbc:listOfFluxObjectives"]["fbc:fluxObjective"].map(fluxObjective => ({
        reaction: fluxObjective["fbc:reaction"],
        coefficient: this.parseValue(fluxObjective["fbc:coefficient"])
      }))
    }));

    this['fbc:geneProducts'] = data['fbc:listOfGeneProducts']['fbc:geneProduct'].map(geneProduct => ({
      id: geneProduct["fbc:id"],
      name: geneProduct["fbc:name"],
      label: geneProduct["fbc:label"]
    }));
  }

  parseValue(value: string) {
    if (isNaN(parseFloat(value))) {
      const parameter = this.parameters.find(parameter => parameter.id === value);
      if (!parameter) throw new Error(`Parameter ${value} not found`);

      const unit = this.unitDefinitions.find(unitDefinition => unitDefinition.id === parameter.units);
      if (!unit) throw new Error(`Unit ${parameter.units} not found`);

      let n = parameter.value;

      for (const u of unit.units) {
        n *= Math.pow(u.multiplier, u.exponent);
      }

      return n;
    } else {
      return parseFloat(value);
    }
  }

  getReaction(name: string) {
    return this.reactions.find(reaction => reaction.name === name);
  }

  createReaction(id: string, name: string, reversible: boolean): SBMLReaction {
    const reaction = {
      id,
      name,
      reversible,
      'fbc:lowerFluxBound': 0,
      'fbc:upperFluxBound': 0,
      listOfReactants: [],
      listOfProducts: []
    };

    this.reactions.push(reaction);

    return reaction;
  }

  removeReaction(name: string) {
    const index = this.reactions.findIndex(reaction => reaction.name === name);
    if (index === -1) {
      //throw new Error(`Reaction ${name} not found`);
      return;
    }

    this.reactions.splice(index, 1);
  }

  getSpecies(id: string) {
    return this.species.find(species => species.id === id);
  }

  createSpecies(id: string, name: string, compartment: string, constant: boolean, initialAmount: number, substanceUnits: string): SBMLSpecies {
    const species = {
      id,
      name,
      compartment,
      constant,
      initialAmount,
      substanceUnits
    };

    this.species.push(species);

    return species;
  }
}

export interface SBML {
  version: string;
  sboTerm: string;
  'fbc:required': boolean;
  level: string;
  model: SBMLModel;
}

export function parseSBML(text: string): SBML {
  const parser = new XMLParser({
    ignoreAttributes: false, attributeNamePrefix: '', isArray: (tagName, jPath) => {
      const parts = jPath.split('.');
      return parts.length >= 2 && parts[parts.length - 2].includes('listOf');
    }
  });
  const result = parser.parse(text);

  const model = new SBMLModel(result.sbml.model);

  const sbml = {
    version: result.sbml.version,
    sboTerm: result.sbml.sboTerm,
    'fbc:required': result.sbml['fbc:required'] === 'true',
    level: result.sbml.level,
    model
  };

  return sbml;
}