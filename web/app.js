/* Browser-side execution of the same Python statistics as the CLI. */
"use strict";
const CDN = "https://cdn.jsdelivr.net/pyodide/v0.28.2/full/";
const byId = id => document.getElementById(id);
let python = null;
let exportedCsv = null;

function status(message, kind) {
  const el = byId("runtime");
  el.textContent = message;
  el.className = "status" + (kind ? " " + kind : "");
}
function fmt(value) {
  if (typeof value !== "number" || !Number.isFinite(value)) return "—";
  if (value === 0) return "0";
  if (Math.abs(value) < 0.0001) return value.toExponential(3);
  return Number(value.toPrecision(6)).toString();
}
function counts() {
  const result = ["aa","ab","bb"].map(id => {
    const input = byId(id);
    if (!input.value.trim() || !input.validity.valid) {
      throw Error("Enter whole-number counts between 0 and 1,000,000.");
    }
    const n = Number(input.value);
    if (!Number.isSafeInteger(n) || n < 0 || n > 1000000) {
      throw Error("Counts must be non-negative whole numbers up to 1,000,000.");
    }
    return n;
  });
  if (result.every(n => n === 0)) throw Error("Enter at least one observed genotype.");
  return result;
}
function appendCell(row, value, tag) {
  const cell = document.createElement(tag || "td");
  cell.textContent = String(value);
  row.append(cell);
}
function showAnalysis(data) {
  const allele = data.allele_frequencies;
  const chi = data.chi_square_test;
  const exact = data.exact_test;
  const f = data.inbreeding_coefficient;
  const monomorphic = allele.p === 0 || allele.q === 0;
  byId("empty").classList.add("hidden");
  byId("results").classList.remove("hidden");
  [["p",allele.p],["q",allele.q],["exact",exact.p_value],
    ["chi-p",chi.p_value],["chi",chi.chi_square],["f",f.F]]
    .forEach(pair => { byId(pair[0]).textContent = fmt(pair[1]); });
  const interpretation = monomorphic ?
    "Monomorphic locus: HWE departure is not informative." :
    exact.p_value <= 0.05 ?
    "Exact test: evidence of departure from HWE at α = 0.05." :
    "Exact test: no statistically significant departure from HWE at α = 0.05.";
  byId("interpretation").textContent =
    interpretation + " Sample size: n = " + allele.n + ".";
  byId("genotype-rows").replaceChildren();
  ["AA","Aa","aa"].forEach(g => {
    const tr = document.createElement("tr");
    appendCell(tr, g, "th");
    appendCell(tr, chi.observed[g]);
    appendCell(tr, fmt(chi.expected[g]));
    byId("genotype-rows").append(tr);
  });
  const warnings = [];
  if (monomorphic) warnings.push("Monomorphic data cannot meaningfully assess equilibrium.");
  warnings.push(...(data.diagnosis.warnings || []));
  const box = byId("warnings");
  box.replaceChildren();
  warnings.forEach(w => {
    const p = document.createElement("p");
    p.className = "warning";
    p.textContent = w;
    box.append(p);
  });
}
function runSingle(event) {
  event.preventDefault();
  byId("single-error").textContent = "";
  if (!python) return;
  try {
    const [a,b,c] = counts();
    const result = python.runPython(
      "json.dumps(hwe.full_hwe_analysis(" + a + "," + b + "," + c + "))");
    showAnalysis(JSON.parse(result));
  } catch (error) {
    byId("single-error").textContent = error.message;
  }
}
async function runBatch() {
  byId("batch-status").textContent = "";
  exportedCsv = null;
  byId("download").disabled = true;
  byId("batch-preview").replaceChildren();
  const file = byId("upload").files[0];
  if (!file) { byId("batch-status").textContent = "Select a CSV file."; return; }
  if (file.size > 2000000) {
    byId("batch-status").textContent = "The maximum supported CSV file size is 2 MB.";
    return;
  }
  byId("process").disabled = true;
  byId("batch-status").textContent = "Analyzing genotype counts locally…";
  try {
    python.FS.writeFile("uploaded.csv", await file.text());
    const result = python.runPython(
      'json.dumps(hwe.process_batch("uploaded.csv", "output.csv"))');
    const rows = JSON.parse(result);
    exportedCsv = python.FS.readFile("output.csv", {encoding:"utf8"});
    byId("download").disabled = false;
    const significant = rows.filter(r => r.bonferroni_significant).length;
    byId("batch-status").textContent = rows.length + " loci analyzed; " +
      significant + " Bonferroni-significant (α=0.05).";
    if (rows.length) {
      const table = document.createElement("table");
      const head = document.createElement("thead");
      const headings = document.createElement("tr");
      ["Locus","Chi-square p","F","Bonferroni"].forEach(
        h => appendCell(headings,h,"th"));
      head.append(headings); table.append(head);
      const body = document.createElement("tbody");
      rows.slice(0,10).forEach((r,i) => {
        const tr = document.createElement("tr");
        appendCell(tr, r.locus || r.id || String(i+1));
        appendCell(tr, fmt(Number(r.p_value)));
        appendCell(tr, fmt(Number(r.F)));
        appendCell(tr, r.bonferroni_significant ? "Yes" : "No");
        body.append(tr);
      });
      table.append(body);
      byId("batch-preview").append(table);
      if (rows.length > 10) byId("batch-status").textContent += " Previewing first 10.";
    }
  } catch (error) {
    byId("batch-status").textContent = "CSV processing failed: " + error.message;
  } finally {
    byId("process").disabled = !python;
  }
}
function downloadResult() {
  if (exportedCsv === null) return;
  const url = URL.createObjectURL(
    new Blob([exportedCsv], {type:"text/csv;charset=utf-8"}));
  const a = document.createElement("a");
  a.href = url;
  a.download = "hwe-results.csv";
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1500);
}
async function initialize() {
  byId("single-form").addEventListener("submit", runSingle);
  byId("example").addEventListener("click", () => {
    byId("aa").value = "40";
    byId("ab").value = "20";
    byId("bb").value = "40";
    if (python) runSingle(new Event("submit"));
  });
  byId("process").addEventListener("click", runBatch);
  byId("download").addEventListener("click", downloadResult);
  try {
    if (typeof loadPyodide !== "function") {
      throw Error("Pyodide script could not be loaded.");
    }
    python = await loadPyodide({indexURL: CDN});
    let response = await fetch("./hwe.py");
    // In the source checkout, the Python file is one directory above web/.
    if (!response.ok && location.pathname.endsWith("/web/")) {
      response = await fetch("../hwe.py");
    }
    if (!response.ok) throw Error("Could not fetch hwe.py (HTTP " + response.status + ").");
    python.FS.writeFile("hwe.py", await response.text());
    python.runPython("import json, hwe");
    status("Python statistical engine ready — analysis runs locally on this device.", "ready");
    byId("analyze").disabled = false;
    byId("process").disabled = false;
    runSingle(new Event("submit"));
  } catch (error) {
    python = null;
    status("Could not start local Python: " + error.message +
      " Check network access and reload.", "error");
  }
}
document.addEventListener("DOMContentLoaded", initialize);
