import { CME } from "./cme/CME";
import Chart from 'chart.js/auto'
import { partTomM } from "./ode/Reactions";

const DECIMATION = 100;

const times = { data: [], view: [] };
let lastUpdateTime = 0;

let volumeChart: Chart;
const volumeData = { data: [], view: [] };

let metaboliteChart: Chart;
const metaboliteData = {
  'M_atp_c': { data: [], view: [] },
  'M_adp_c': { data: [], view: [] },
  'M_gtp_c': { data: [], view: [] },
  'M_ctp_c': { data: [], view: [] },
  'M_utp_c': { data: [], view: [] },
};

let atpChart: Chart;
const atpData = {
  'ATP_trsc': { data: [], view: [] },
  'ATP_translat': { data: [], view: [] },
  'ATP_mRNAdeg': { data: [], view: [] },
  'ATP_ptndeg': { data: [], view: [] },
  'ATP_DNArep': { data: [], view: [] },
  'ATP_transloc': { data: [], view: [] },
};

const commonOptions: any = {
  animation: false,
};

function multilineChart(id: string, data: { [key: string]: { view: number[] } }, yScale = "linear") {
  return new Chart(
    document.getElementById(id) as HTMLCanvasElement,
    {
      type: 'line',
      data: {
        labels: times.view,
        datasets: Object.keys(data).map(id => ({
          label: id,
          data: data[id].view,
          pointRadius: 0
        }))
      },
      options: {
        ...commonOptions,
        scales: {
          y: {
            type: yScale,
            beginAtZero: true,
            title: {
              display: true,
              text: 'Concentration'
            }
          },
          x: {
            type: "linear",
            title: {
              display: true,
              text: 'Time'
            }
          }
        }
      }
    }
  );
}

export function initUI() {
  metaboliteChart = multilineChart('metabolites', metaboliteData);
  atpChart = multilineChart('atp', atpData, "logarithmic");

  volumeChart = new Chart(
    document.getElementById('volume') as HTMLCanvasElement,
    {
      type: 'line',
      data: {
        labels: times.view,
        datasets: [{
          label: 'Volume',
          data: volumeData.view,
          pointRadius: 0
        }]
      },
      options: {
        ...commonOptions,
        scales: {
          y: {
            type: "linear",
            title: {
              display: true,
              text: 'Volume'
            }
          },
          x: {
            type: "linear",
            title: {
              display: true,
              text: 'Time'
            }
          }
        }
      }
    }
  );
}

function decimate(data: number[], size: number): number[] {
  if (data.length <= size) return data;

  const decimated = [];

  const windowSize = Math.ceil(data.length / size);
  let startIndex = 0;

  while (decimated.length < size) {
    const endIndex = Math.min(startIndex + windowSize, data.length);
    if (endIndex <= startIndex) break;

    let sum = 0;

    for (let i = startIndex; i < endIndex; i++) {
      sum += data[i];
    }

    decimated.push(sum / (endIndex - startIndex));

    startIndex += windowSize;
  }

  return decimated;
}

function overwrite(data: number[], newData: number[]) {
  for (let i = 0; i < newData.length; i++) {
    if (i < data.length) {
      data[i] = newData[i];
    } else {
      data.push(newData[i]);
    }
  }

  if (newData.length < data.length) {
    data.splice(newData.length, data.length - newData.length);
  }
}

export function updateUI(t: number, data: Map<string, number>) {
  times.data.push(Math.round(t * 100) / 100);
  volumeData.data.push(data.get("CellV"));

  for (const id of Object.keys(metaboliteData)) {
    metaboliteData[id].data.push(partTomM(data.get(id), data));
  }

  for (const id of Object.keys(atpData)) {
    atpData[id].data.push(partTomM(data.get(id), data));
  }

  if (t - lastUpdateTime < 0.4) return false;

  overwrite(times.view, decimate(times.data, DECIMATION));
  overwrite(volumeData.view, decimate(volumeData.data, DECIMATION));

  for (const id of Object.keys(metaboliteData)) {
    overwrite(metaboliteData[id].view, decimate(metaboliteData[id].data, DECIMATION));
  }

  for (const id of Object.keys(atpData)) {
    overwrite(atpData[id].view, decimate(atpData[id].data, DECIMATION));
  }

  volumeChart.update("none");
  metaboliteChart.update("none");
  atpChart.update("none");

  lastUpdateTime = t;

  return true;
}