// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

const WIDTH = 320;
const HEIGHT = 180;
const ORBIT_FRAME_INTERVAL = 1000 / 15;
const SETTLE_DELAY = 180;

const canvas = document.querySelector("#scene");
const context = canvas.getContext("2d", {alpha: false});
const previewCanvas = document.createElement("canvas");
const previewContext = previewCanvas.getContext("2d", {alpha: false});
const status = document.querySelector("#status");
const modelInfo = document.querySelector("#model-info");
const pauseButton = document.querySelector("#pause");
const recordButton = document.querySelector("#record");
const fileInput = document.querySelector("#model-file");
const clipInput = document.querySelector("#clip");
const clipPosition = document.querySelector("#clip-position");
const targetInputs = ["#target-x", "#target-y", "#target-z"]
  .map((selector) => document.querySelector(selector));
const goToButton = document.querySelector("#go-to");
const modelCenterButton = document.querySelector("#model-center");

const worker = new Worker("render-worker.js", {type: "module"});
let ready = false;
let renderPending = false;
let paused = false;
let dirty = true;
let azimuth = 0.8;
let elevation = 0.42;
let distance = 2.7;
let dragging = false;
let previousPointer = null;
let previousTime = performance.now();
let modelCenter = [0, 0, 0];
let interactionUntil = 0;
let nextOrbitFrame = 0;
let needsFullFrame = true;
let recording = false;
let pausedBeforeRecording = false;

function markInteraction() {
  interactionUntil = performance.now() + SETTLE_DELAY;
  needsFullFrame = true;
  dirty = true;
}

function setCoordinateInputs(coordinates) {
  coordinates.forEach((value, index) => {
    targetInputs[index].value = Number(value.toPrecision(8));
  });
}

function showModel(name, model) {
  modelInfo.textContent = `${name} · ${model.cells} cells · ${model.surfaces} surfaces`;
  modelCenter = model.center;
  setCoordinateInputs(model.center);
}

worker.onmessage = ({data}) => {
  if (data.type === "ready") {
    ready = true;
    const mode = data.threaded ? `OpenMP · ${data.threads} workers` : "single thread";
    status.textContent = `Ready · ${mode}`;
    showModel("pin-cluster.mcnp", data.model);
    dirty = true;
  } else if (data.type === "loaded") {
    ready = true;
    status.textContent = "MCNP geometry loaded";
    showModel(data.name, data.model);
    dirty = true;
  } else if (data.type === "targeted") {
    status.textContent = `Target · ${data.target.map((value) => value.toPrecision(6)).join(", ")}`;
    dirty = true;
  } else if (data.type === "frame") {
    const image = new ImageData(
      new Uint8ClampedArray(data.pixels), data.width, data.height);
    if (data.preview) {
      previewCanvas.width = data.width;
      previewCanvas.height = data.height;
      previewContext.putImageData(image, 0, 0);
      context.imageSmoothingEnabled = true;
      context.drawImage(previewCanvas, 0, 0, WIDTH, HEIGHT);
    } else {
      context.putImageData(image, 0, 0);
    }
    const quality = data.preview ? "preview" : "full";
    document.querySelector("#render-time").textContent =
      `${data.renderMs.toFixed(1)} ms · ${quality}`;
    renderPending = false;
  } else if (data.type === "error") {
    status.textContent = `Error · ${data.message}`;
    renderPending = false;
    ready = true;
  }
};

const forceSingleThread = new URLSearchParams(location.search).has("single");
worker.postMessage({type: "init", threaded: crossOriginIsolated && !forceSingleThread});

function animate(now) {
  const delta = Math.min((now - previousTime) / 1000, 0.1);
  previousTime = now;
  if (!paused && !dragging) {
    azimuth += delta * 0.22;
    dirty = true;
  }
  const moving = dragging || !paused || now < interactionUntil;
  const preview = moving && !recording;
  if (!moving && needsFullFrame) dirty = true;
  const throttleReady = !preview || now >= nextOrbitFrame;
  if (ready && dirty && !renderPending && throttleReady) {
    dirty = false;
    renderPending = true;
    if (preview) {
      needsFullFrame = true;
      nextOrbitFrame = now + ORBIT_FRAME_INTERVAL;
    } else {
      needsFullFrame = false;
    }
    worker.postMessage({
      type: "render", azimuth, elevation, distance,
      clip: clipInput.checked,
      clipFraction: Number(clipPosition.value),
      preview,
    });
  }
  requestAnimationFrame(animate);
}

canvas.addEventListener("pointerdown", (event) => {
  dragging = true;
  previousPointer = event;
  canvas.setPointerCapture(event.pointerId);
  markInteraction();
});
canvas.addEventListener("pointermove", (event) => {
  if (!dragging) return;
  azimuth -= (event.clientX - previousPointer.clientX) * 0.012;
  elevation = Math.max(-1.35, Math.min(1.35,
    elevation + (event.clientY - previousPointer.clientY) * 0.012));
  previousPointer = event;
  markInteraction();
});
canvas.addEventListener("pointerup", () => {
  dragging = false;
  markInteraction();
});
canvas.addEventListener("wheel", (event) => {
  event.preventDefault();
  distance = Math.max(1.5, Math.min(7, distance * Math.exp(event.deltaY * 0.001)));
  markInteraction();
}, {passive: false});

pauseButton.addEventListener("click", () => {
  paused = !paused;
  pauseButton.textContent = paused ? "Orbit" : "Pause";
  needsFullFrame = true;
  dirty = true;
});
clipInput.addEventListener("change", markInteraction);
clipPosition.addEventListener("input", markInteraction);

function goToCoordinates(coordinates) {
  if (!coordinates.every(Number.isFinite)) {
    status.textContent = "Error · X, Y, and Z must be valid numbers";
    return;
  }
  worker.postMessage({
    type: "target", x: coordinates[0], y: coordinates[1], z: coordinates[2],
  });
}

goToButton.addEventListener("click", () => {
  goToCoordinates(targetInputs.map((input) => Number(input.value)));
});
targetInputs.forEach((input) => input.addEventListener("keydown", (event) => {
  if (event.key === "Enter") goToButton.click();
}));
modelCenterButton.addEventListener("click", () => {
  setCoordinateInputs(modelCenter);
  goToCoordinates(modelCenter);
});

fileInput.addEventListener("change", async () => {
  const file = fileInput.files[0];
  if (!file) return;
  ready = false;
  status.textContent = `Parsing ${file.name}…`;
  const buffer = await file.arrayBuffer();
  worker.postMessage({type: "load", name: file.name, buffer}, [buffer]);
});

recordButton.addEventListener("click", () => {
  const stream = canvas.captureStream(30);
  const preferred = "video/webm;codecs=vp9";
  const options = MediaRecorder.isTypeSupported(preferred) ? {mimeType: preferred} : {};
  const recorder = new MediaRecorder(stream, options);
  const chunks = [];
  pausedBeforeRecording = paused;
  recording = true;
  paused = false;
  pauseButton.textContent = "Pause";
  dirty = true;
  recorder.ondataavailable = ({data}) => { if (data.size) chunks.push(data); };
  recorder.onstop = () => {
    const blob = new Blob(chunks, {type: recorder.mimeType || "video/webm"});
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "libalea-mcnp-orbit.webm";
    link.click();
    setTimeout(() => URL.revokeObjectURL(link.href), 1000);
    recordButton.disabled = false;
    recordButton.textContent = "Record 10s";
    recording = false;
    paused = pausedBeforeRecording;
    pauseButton.textContent = paused ? "Orbit" : "Pause";
    needsFullFrame = true;
    dirty = true;
  };
  recordButton.disabled = true;
  recordButton.textContent = "Recording…";
  recorder.start();
  setTimeout(() => recorder.stop(), 10_000);
});

requestAnimationFrame(animate);
