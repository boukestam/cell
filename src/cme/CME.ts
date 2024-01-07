class Reaction {
  reactants: string[];
  products: string[];
  rate: number;
  propensity: number;

  constructor(reactants: string[], products: string[], rate: number) {
    this.reactants = reactants;
    this.products = products;
    this.rate = rate;
    this.propensity = 0;
  }

  updatePropensity(species: Map<string, number>) {
    this.propensity = this.reactants.reduce((acc, reactant) => {
      return acc * (species.get(reactant) || 0);
    }, this.rate);
  }
}

export class CME {

  species: Map<string, number>;
  reactions: Reaction[];
  speciesReactions: Map<string, Reaction[]>;

  constructor() {
    this.species = new Map<string, number>();
    this.reactions = [];
    this.speciesReactions = new Map<string, Reaction[]>();
  }

  defineSpecies(species: string | string[]) {
    if (typeof species === 'string') {
      this.species.set(species, 0);
    } else {
      for (const s of species) {
        this.species.set(s, 0);
      }
    }
  }

  addReaction(reactants: string[], products: string[], rate: number) {
    if (isNaN(rate)) {
      throw new Error(`Rate ${rate} is not a number`);
    }

    if (!Array.isArray(reactants)) {
      throw new Error(`Reactants ${reactants} is not an array`);
    }

    for (const species of [...reactants, ...products]) {
      if (!this.species.has(species)) {
        throw new Error(`Species ${species} does not exist`);
      }
    }

    const reaction = new Reaction(reactants, products, rate);
    this.reactions.push(reaction);

    for (const species of reactants) {
      const reactions = this.speciesReactions.get(species) || [];
      this.speciesReactions.set(species, [...reactions, reaction]);
    }
  }

  addParticles(species: string, amount: number) {
    let count = this.species.get(species);

    if (count === undefined) {
      throw new Error(`Species ${species} does not exist`);
    }

    this.species.set(species, count + amount);
  }

  solve(totalTime: number, hookInterval: number, hook: (time: number) => void) {
    let time = 0;

    // Calculate propensities
    for (const reaction of this.reactions) {
      reaction.updatePropensity(this.species);
    }

    let lastHookTime = 0;

    while (time < totalTime) {
      const result = this.solveStep(time);
      time = result.time;

      if (time - lastHookTime >= hookInterval) {
        hook(time);
        lastHookTime = time;
      }
    }
  }

  solveStep(time: number) {
    const totalPropensity = this.reactions.reduce((acc, reaction) => acc + reaction.propensity, 0);

    // Calculate time to next reaction
    const r1 = Math.random();
    const timeToNextReaction = -Math.log(r1) / totalPropensity;
    time += timeToNextReaction;

    // Choose a reaction
    const r2 = Math.random();
    let reactionIndex = 0;
    let cumulativePropensity = 0;

    while (cumulativePropensity < r2 * totalPropensity) {
      cumulativePropensity += this.reactions[reactionIndex].propensity;
      reactionIndex++;
    }

    // Update counts
    const { reactants, products } = this.reactions[reactionIndex - 1];

    reactants.forEach(reactant => {
      const count = this.species.get(reactant) || 0;
      this.species.set(reactant, count - 1);
    });

    products.forEach(product => {
      const count = this.species.get(product) || 0;
      this.species.set(product, count + 1);
    });

    // Update propensities
    for (const species of [...reactants, ...products]) {
      for (const reaction of this.speciesReactions.get(species) || []) {
        reaction.updatePropensity(this.species);
      }
    }

    return { time, reactants, products };
  }

  toString() {
    // Write the header
    let output = "# CME Simulation File\n";

    // Write the species and counts
    output += "# Species\n";
    for (const [species, count] of this.species.entries()) {
      output += `species ${species}: ${count}\n`;
    }

    // Write the reactions
    output += "# Reactions\n";
    for (const reaction of this.reactions) {
      output += `reaction ${reaction.reactants.join(" + ")} -> ${reaction.products.join(" + ")}: ${reaction.rate}\n`;
    }

    return output;
  }
}