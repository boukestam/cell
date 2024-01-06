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

// Create our simulation object
const sim = new CME();

loadData().then((data) => {
  console.log("Loaded", data);

  initDicts();

  initIC(sim);

  populate(sim);

  initReplication(sim);
  addReplication(sim, data.genome3A, ModelSpecies);

  const metabolicData = defineMetabolicReactions();
  console.log(metabolicData);

  // Run simulation

  console.log("Running simulation");
  const startTime = Date.now();

  console.log("ATP", sim.species.get("M_atp_c"));

  console.log([...sim.species.keys()].map(k => `${k}: ${Math.floor(sim.species.get(k))}`).join("\n"));

  sim.solve(10, 1, (time) => {
    console.log("ATP", sim.species.get("M_atp_c"));

    const start = Date.now();

    const model = new ODE();

    model.addRateForm("zeroOrderOnOff", "$onoff * $K");

    addReactionsToModel(model, sim.species, metabolicData);

    lipidPatch(model, sim.species);
    pppPatch(model);

    console.log("ODE model", model);

    const integrateStart = Date.now();
    const newValues = integrate(model, 1, 0.1, "lsoda");
    console.log("Integrated ODE in", Date.now() - integrateStart, "ms");

    writeODEtoCME(sim.species, newValues);

    const end = Date.now();
    console.log("Finished ODE at time", time, "in", end - start, "ms");
  });

  console.log("Finished simulation");
  const endTime = Date.now();
  console.log("Simulation took", endTime - startTime, "ms");
});