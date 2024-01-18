import { CME } from "./cme/CME";
import { loadData } from "./data/Data";
import { ModelSpecies, initDicts } from "./Globals";
import { initIC } from "./cme/IC";
import { lipidPatch } from "./ode/LipidPatch";
import { addReactionsToModel, defineMetabolicReactions } from "./ode/MetabolicReactions";
import { ODE } from "./ode/ODE";
import { pppPatch } from "./ode/PppPatch";
import { populate } from "./cme/Populate";
import { addReplication, initReplication } from "./cme/Replication";
import { integrate } from "./ode/Integrator";
import { writeODEtoCME } from "./ode/ODEtoCME";
import { initUI, updateUI } from "./ui";
import { mMtoPart } from "./ode/Reactions";

async function sleep(ms: number) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function run() {
  initUI();

  const data = await loadData();

  console.log("Loaded data");

  initDicts();

  // Create our simulation object
  const sim = new CME();

  initIC(sim);

  populate(sim);

  initReplication(sim);
  addReplication(sim, data.genome3A, ModelSpecies);

  //console.log([...sim.species.entries()].map(([key, value]) => key + ": " + value).join("\n"));
  //console.log(sim.reactions.map(r => r.toString()).join("\n"));

  const metabolicData = defineMetabolicReactions();

  // Run simulation

  console.log("Running simulation");
  const startTime = Date.now();

  let i = 0;

  await sim.solve(120 * 60, 1, async (time) => {
    const uiData = new Map(sim.species);

    const model = new ODE();
    model.addRateForm("zeroOrderOnOff", "$onoff * $K");

    addReactionsToModel(model, sim.species, metabolicData);
    lipidPatch(model, sim.species);
    pppPatch(model);

    const newValues = integrate(model, 1, 0.1, "lsoda");

    for (const [key, value] of newValues.entries()) {
      if (!uiData.has(key)) continue;
      uiData.set("ODE_" + key, mMtoPart(value, sim.species) - uiData.get(key));
    }

    writeODEtoCME(sim.species, newValues);

    if ((i++) > 0) {
      if (updateUI(time / 60, uiData)) await sleep(0);
    }
  });

  console.log("Finished simulation");
  const endTime = Date.now();
  console.log("Simulation took", endTime - startTime, "ms");
}

run().catch(console.error);