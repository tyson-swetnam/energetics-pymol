/* Energetics — two-pane gallery with switchable Mol* / 3Dmol.js viewers. */
"use strict";

// diagnostics: count runtime errors (read in tests) without spamming the console
window.__errs = 0;
window.addEventListener("error", () => { window.__errs++; });
window.addEventListener("unhandledrejection", (e) => { window.__errs++; e.preventDefault(); });

// category accent, ordered low -> high energy density (matches the CSS energy gradient)
const ENERGY_ORDER = ["Sugar monomers", "Hemicellulose", "Cellulose", "Lignin",
                      "Fats", "Extractives", "Fossil carbon", "Combustion"];
const CAT_COLOR = {
  "Sugar monomers": "#8aa06a", "Hemicellulose": "#9ba06a", "Cellulose": "#b9a96a",
  "Lignin": "#d98a3a", "Fats": "#e0712d", "Extractives": "#f4a93b",
  "Fossil carbon": "#7a5238", "Combustion": "#c2402f",
};

const $ = (id) => document.getElementById(id);
const esc = (s) => s.replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const absUrl = (p) => new URL(p, location.href).href;

const state = {
  data: null, byId: {}, current: null,
  engine: "molstar", bg: "dark", spin: false,
  molstar: null, gl: null, host: null,
  molstarQueue: Promise.resolve(), // serialize Mol* load/clear so they never race
};

// ------------------------------------------------------------------ load + build
async function init() {
  state.data = await (await fetch("molecules.json")).json();
  state.data.molecules.forEach((m) => (state.byId[m.id] = m));
  $("subtitle").textContent = state.data.meta.subtitle;
  $("metaStat").textContent = `${state.data.meta.count} molecules · ${state.data.meta.paper.split("—")[0].trim()}`;
  buildSidebar();
  wireControls();

  state.host = document.createElement("div");
  state.host.className = "engine-host";
  $("viewport").appendChild(state.host);

  const start = location.hash.slice(1);
  if (start && state.byId[start]) select(start);
}

function buildSidebar(filter = "") {
  const nav = $("nav");
  nav.innerHTML = "";
  const q = filter.trim().toLowerCase();
  for (const cat of ENERGY_ORDER) {
    let items = state.data.molecules.filter((m) => m.category === cat);
    if (q) items = items.filter((m) =>
      m.name.toLowerCase().includes(q) || m.formula.toLowerCase().includes(q) ||
      m.id.includes(q) || cat.toLowerCase().includes(q));
    if (!items.length) continue;

    const sec = document.createElement("div");
    sec.className = "cat";
    sec.innerHTML = `<div class="cat-head"><span class="cat-dot" style="background:${CAT_COLOR[cat]}"></span>${cat}<span class="cat-count">${items.length}</span></div>`;
    for (const m of items) {
      const b = document.createElement("button");
      b.className = "mol" + (state.current === m.id ? " active" : "");
      b.dataset.id = m.id;
      b.innerHTML =
        `<img class="mol-thumb" src="${m.img}" alt="" loading="lazy">` +
        `<span class="mol-meta"><span class="mol-nm">${esc(prettyName(m.name))}</span>` +
        `<span class="mol-fm">${esc(m.formula)}</span></span>` +
        `<span class="mol-energy" style="background:${CAT_COLOR[cat]}"></span>`;
      b.addEventListener("click", () => select(m.id));
      sec.appendChild(b);
    }
    nav.appendChild(sec);
  }
  if (!nav.children.length)
    nav.innerHTML = `<p style="padding:16px;color:var(--paper-faint);font-size:13px">No molecules match “${esc(filter)}”.</p>`;
}

// drop only a trailing "(...)" qualifier for a tidier list label (keep inline locants)
function prettyName(n) { return n.replace(/\s*\([^)]*\)\s*$/, ""); }

// ------------------------------------------------------------------ selection
async function select(id) {
  const m = state.byId[id];
  if (!m) return;
  state.current = id;
  history.replaceState(null, "", "#" + id);

  document.querySelectorAll(".mol").forEach((el) =>
    el.classList.toggle("active", el.dataset.id === id));

  $("molName").textContent = m.name;
  $("molFormula").textContent = m.formula;
  $("molTag").textContent = m.category;
  $("ctlPdb").href = m.pdb;
  $("ctlPdb").setAttribute("download", id + ".pdb");
  $("viewportEmpty").style.display = "none";

  renderReadout(m);
  await showMolecule(m);
}

function renderReadout(m) {
  $("readout").hidden = false;
  const meta = state.data.meta;
  const fill = $("gaugeFill"), ref = $("gaugeRef");
  if (m.hhvKind !== "fuel") {
    fill.style.width = "0%";
    ref.style.display = "none";
    $("hhvNum").textContent = m.hhvKind === "oxidant" ? "oxidant" : "0 MJ/kg";
    $("gaugeCap").textContent = m.hhvKind === "oxidant"
      ? "Molecular oxygen — the oxidant that releases biomass energy on combustion."
      : "Combustion product — fully oxidized, no remaining fuel value.";
  } else {
    fill.style.width = Math.min(100, (m.hhv / meta.gaugeMax) * 100) + "%";
    ref.style.display = "block";
    ref.style.left = (meta.wholeTreeHHV / meta.gaugeMax) * 100 + "%";
    $("hhvNum").textContent = m.hhv.toFixed(1) + " MJ/kg";
    $("gaugeCap").textContent =
      `Approx. higher heating value. Whole-tree Pinus mean ≈ ${meta.wholeTreeHHV} MJ/kg (paper, Table EC).`;
  }
  $("facts").innerHTML =
    `<dt>Residue</dt><dd class="mono">${esc(m.resname)}</dd>` +
    `<dt>Class</dt><dd>${esc(m.category)}</dd>` +
    (m.size !== "na" ? `<dt>Fragment</dt><dd>${esc(m.size)}</dd>` : "") +
    `<dt>SMILES</dt><dd class="mono">${esc(m.smiles)}</dd>`;
  $("molNote").textContent = m.note;
}

// ------------------------------------------------------------------ viewers
async function showMolecule(m) {
  if (state.engine === "molstar") {
    if (!state.molstar) await mountMolstar();
    await loadMolstar(m);
  } else {
    if (!state.gl) mount3dmol();
    await load3dmol(m);
  }
}

async function setEngine(engine) {
  if (engine === state.engine) return;
  teardown();
  state.engine = engine;
  $("engMolstar").setAttribute("aria-pressed", String(engine === "molstar"));
  $("eng3dmol").setAttribute("aria-pressed", String(engine === "threedmol"));
  $("engineHint").textContent = engine === "molstar"
    ? "Mol* · drag to rotate · scroll to zoom" : "3Dmol · drag to rotate · scroll to zoom";
  if (state.current) await showMolecule(state.byId[state.current]);
}

function teardown() {
  try { state.molstar?.plugin?.dispose(); } catch (e) {}
  state.molstar = null;
  try { state.gl?.clear(); } catch (e) {}
  state.gl = null;
  if (state.host) state.host.innerHTML = "";
}

/* ---- Mol* ---- */
async function mountMolstar() {
  state.molstar = await molstar.Viewer.create(state.host, {
    layoutIsExpanded: false, layoutShowControls: false, layoutShowRemoteState: false,
    layoutShowSequence: false, layoutShowLog: false, layoutShowLeftPanel: false,
    viewportShowExpand: false, viewportShowControls: false, collapseLeftPanel: true,
    pdbProvider: "rcsb", emdbProvider: "rcsb",
  });
  applyMolstarBg();
}
// serialize every Mol* clear/load so a new selection never interrupts an in-flight one
function loadMolstar(m) {
  state.molstarQueue = state.molstarQueue.then(() => _loadMolstar(m)).catch((e) => console.warn("Mol* load:", e));
  return state.molstarQueue;
}
async function _loadMolstar(m) {
  const v = state.molstar;
  if (!v) return;
  applyMolstarSpin(false);          // stop spinning before tearing down the scene
  await v.plugin.clear();
  if (m.id !== state.current) return; // a newer selection won; skip this stale load
  await v.loadStructureFromUrl(absUrl(m.pdb), "pdb", false);
  applyMolstarBg();
  applyMolstarSpin();
}
function applyMolstarBg() {
  const c = state.molstar?.plugin?.canvas3d;
  if (!c) return;
  try { c.setProps({ renderer: { backgroundColor: state.bg === "light" ? 0xfbfaf7 : 0x0d0b08 } }); }
  catch (e) { /* canvas not ready */ }
}
function applyMolstarSpin(force) {
  const c = state.molstar?.plugin?.canvas3d;
  if (!c) return;
  const on = force === undefined ? state.spin : force;
  try { c.setProps({ trackball: { animate: on ? { name: "spin", params: { speed: 1 } } : { name: "off", params: {} } } }); }
  catch (e) { /* canvas not ready */ }
}

/* ---- 3Dmol.js ---- */
function mount3dmol() {
  state.gl = $3Dmol.createViewer(state.host, { backgroundColor: state.bg === "light" ? "#fbfaf7" : "#0d0b08" });
}
async function load3dmol(m) {
  const v = state.gl;
  v.clear();
  const data = await (await fetch(absUrl(m.pdb))).text();
  v.addModel(data, "pdb");
  v.setStyle({}, { stick: { radius: 0.13, colorscheme: "Jmol" }, sphere: { scale: 0.22, colorscheme: "Jmol" } });
  v.setBackgroundColor(state.bg === "light" ? "#fbfaf7" : "#0d0b08");
  v.zoomTo();
  v.render();
  v.spin(state.spin ? "y" : false);
}

// ------------------------------------------------------------------ controls
function wireControls() {
  $("search").addEventListener("input", (e) => buildSidebar(e.target.value));
  $("engMolstar").addEventListener("click", () => setEngine("molstar"));
  $("eng3dmol").addEventListener("click", () => setEngine("threedmol"));

  $("ctlSpin").addEventListener("click", () => {
    state.spin = !state.spin;
    $("ctlSpin").setAttribute("aria-pressed", String(state.spin));
    if (state.engine === "molstar") applyMolstarSpin();
    else state.gl?.spin(state.spin ? "y" : false);
  });

  $("ctlReset").addEventListener("click", () => {
    if (state.engine === "molstar") state.molstar?.plugin?.canvas3d?.requestCameraReset();
    else { state.gl?.zoomTo(); state.gl?.render(); }
  });

  $("ctlBg").addEventListener("click", () => {
    state.bg = state.bg === "dark" ? "light" : "dark";
    document.body.dataset.viewerBg = state.bg;
    $("ctlBg").textContent = state.bg === "dark" ? "Light background" : "Dark background";
    if (state.engine === "molstar") applyMolstarBg();
    else { state.gl?.setBackgroundColor(state.bg === "light" ? "#fbfaf7" : "#0d0b08"); state.gl?.render(); }
  });

  window.addEventListener("resize", () => { if (state.engine === "threedmol" && state.gl) { state.gl.resize(); state.gl.render(); } });
  window.addEventListener("hashchange", () => {
    const id = location.hash.slice(1);
    if (id && state.byId[id] && id !== state.current) select(id);
  });
}

init().catch((e) => {
  console.error(e);
  document.getElementById("viewportEmpty").innerHTML =
    `<p>Could not load the molecule data.<br><small>${esc(String(e))}</small></p>`;
});
