class Reaction {
  reactants: Map<string, number>;
  products: Map<string, number>;
  rate: number;
  propensity: number;

  constructor(reactants: string[], products: string[], rate: number) {
    const reactantCounts: Map<string, number> = new Map();
    const productCounts: Map<string, number> = new Map();

    for (const reactant of reactants) {
      const count = reactantCounts.get(reactant) || 0;
      reactantCounts.set(reactant, count + 1);
    }

    for (const product of products) {
      const count = productCounts.get(product) || 0;
      productCounts.set(product, count + 1);
    }

    this.reactants = reactantCounts;
    this.products = productCounts;
    this.rate = rate;
    this.propensity = 0;
  }

  updatePropensity(species: Map<string, number>) {
    let propensity = this.rate;

    for (const entry of this.reactants) {
      propensity *= Math.pow(species.get(entry[0]), entry[1]);
    }

    this.propensity = propensity;

    return propensity;
  }

  toString() {
    // 1A + 2B + 3C -> 2D + 1E (K = 0.1)

    return [...this.reactants.entries()].map(([reactant, count]) => `${count}${reactant}`).join(" + ") +
      " -> " +
      [...this.products.entries()].map(([product, count]) => `${count}${product}`).join(" + ") +
      ` (K = ${this.rate.toFixed(6)})`;
  }
}

export class CME {

  species: Map<string, number>;
  reactions: Reaction[];
  speciesReactions: Map<string, Reaction[]>;

  totalPropensity: number;

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

    for (const species of reaction.reactants.keys()) {
      if (this.speciesReactions.has(species)) {
        this.speciesReactions.get(species).push(reaction);
      } else {
        this.speciesReactions.set(species, [reaction]);
      }
    }
  }

  addParticles(species: string, amount: number) {
    let count = this.species.get(species);

    if (count === undefined) {
      throw new Error(`Species ${species} does not exist`);
    }

    this.species.set(species, count + Math.floor(amount));
  }

  async solve(totalTime: number, hookInterval: number, hook: (time: number) => Promise<void>) {
    let time = 0;

    // Calculate propensities
    for (const reaction of this.reactions) {
      reaction.updatePropensity(this.species);
    }

    this.totalPropensity = this.reactions.reduce((acc, reaction) => acc + reaction.propensity, 0);

    let lastHookTime = 0;

    while (time < totalTime) {
      time = this.solveStep(time);

      if (time - lastHookTime >= hookInterval) {
        await hook(time);
        lastHookTime = time;
      }
    }
  }

  solveStep(time: number) {
    // Calculate time to next reaction
    const r1 = Math.random();
    const timeToNextReaction = (1 / this.totalPropensity) * Math.log(1 / r1);
    time += timeToNextReaction;

    // Choose a reaction
    const r2 = Math.random();
    let reactionIndex = 0;
    let cumulativePropensity = 0;

    while (cumulativePropensity < r2 * this.totalPropensity) {
      cumulativePropensity += this.reactions[reactionIndex].propensity;
      reactionIndex++;
    }

    // Update counts
    const { reactants, products } = this.reactions[reactionIndex - 1];

    for (const reactant of reactants) {
      const count = this.species.get(reactant[0]) || 0;
      this.species.set(reactant[0], count - reactant[1]);
    }

    for (const product of products) {
      const count = this.species.get(product[0]) || 0;
      this.species.set(product[0], count + product[1]);
    }

    // Update propensities
    for (const species of [...reactants, ...products]) {
      for (const reaction of this.speciesReactions.get(species[0]) || []) {
        const oldPropensity = reaction.propensity;
        reaction.updatePropensity(this.species);

        // Update total propensity
        this.totalPropensity += reaction.propensity - oldPropensity;
      }
    }

    return time;
  }

  toString() {
    // Write the header
    let output = "# CME Simulation File\n";

    // Write the species and counts
    output += "# Species\n";
    for (const [species, count] of this.species.entries()) {
      output += `${species}: ${count}\n`;
    }

    // Write the reactions
    output += "# Reactions\n";
    for (const reaction of this.reactions) {
      output += reaction.toString() + "\n";
    }

    return output;
  }
}