import { CME } from "./CME";
import { loadData } from "./Data";
import { ModelSpecies, initDicts } from "./Globals";
import { initIC } from "./IC";
import { lipidPatch } from "./LipidPatch";
import { addReactionsToModel, defineMetabolicReactions } from "./MetabolicReactions";
import { ODE } from "./ODE";
import { pppPatch } from "./PppPatch";
import { populate } from "./Populate";
import { addReplication, initReplication } from "./Replication";
import { integrate } from "./Integrator";
import { writeODEtoCME } from "./ODEtoCME";

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

  sim.solve(10, 1, (time) => {
    console.log("ATP", sim.species.get("M_atp_c"));

    const start = Date.now();

    const model = new ODE();

    model.addRateForm("zeroOrderOnOff", "$onoff * $K");

    addReactionsToModel(model, sim.species, metabolicData);

    lipidPatch(model, sim.species);
    pppPatch(model);

    const buildStart = Date.now();
    const derivatives = model.build();

    const values = new Map<string, number>();
    for (const metabolite of model.metabolites.keys()) {
      values.set(metabolite, model.metabolites.get(metabolite).initialValue);
    }

    console.log("Built ODE in", Date.now() - buildStart, "ms", model, derivatives, values);

    const integrateStart = Date.now();
    const newValues = integrate(values, derivatives, 1, 0.001, "euler");
    console.log("Integrated ODE in", Date.now() - integrateStart, "ms");

    writeODEtoCME(sim.species, newValues);

    const end = Date.now();
    console.log("Finished ODE at time", time, "in", end - start, "ms");
  });

  console.log("Finished simulation");
  const endTime = Date.now();
  console.log("Simulation took", endTime - startTime, "ms");
});