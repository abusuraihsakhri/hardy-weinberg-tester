# Hardy Weinberg Tester

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Hardy-Weinberg Equilibrium Tester — Real Implementation

Tests genotype frequencies for HWE using chi-square and exact tests,
calculates inbreeding coefficients, provides deviation diagnostics,
and supports multiple testing correction.

Zero external dependencies — Python stdlib only.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`allele_frequencies()`**: Calculate allele frequencies from genotype counts.

Args:
    hom_ref: Count of homozygous reference (AA) genotypes.
    het: Count of heterozygous (Aa) genotypes.
    hom_alt: Count of homozygous alternate (aa) genotypes.

Returns:
    Dict with p (ref freq), q (alt freq), n (sample size).
- **`expected_genotype_frequencies()`**: Calculate expected Hardy-Weinberg genotype frequencies.

p² (homozygous reference / AA)
2pq (heterozygous / Aa)
q² (homozygous alternate / aa)
- **`expected_genotype_counts()`**: Calculate expected genotype counts under HWE.
- **`chi_square_test()`**: Chi-square goodness-of-fit test for Hardy-Weinberg Equilibrium.

χ² = Σ((observed - expected)² / expected)

Returns:
    Dict with chi_square, p_value, df, in_hwe, expected counts, etc.
- **`exact_test()`**: Exact test for HWE (analogous to Fisher's exact test).

For small samples, computes the probability of observing genotype
counts as or more extreme than observed under HWE.

Uses the multinomial probability approach:
P(genotype | HWE) = n! / (nAA! × nAa! × naa!) × (p²)^nAA × (2pq)^nAa × (q²)^naa

We enumerate possible configurations and sum probabilities.

---

## 📐 Mathematical Formulation & Logic

```text
  calculates inbreeding coefficients, provides deviation diagnostics,
  """Calculate allele frequencies from genotype counts.
  """Calculate expected Hardy-Weinberg genotype frequencies.
  """Calculate expected genotype counts under HWE."""
  """Calculate inbreeding coefficient F = 1 - (Ho / He).
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --hom-ref <value> --het <value> --hom-alt <value> --p-values <value>
```

### Parameter Reference
- `--hom-ref`: Specifies input measurement or parameter value.
- `--het`: Specifies input measurement or parameter value.
- `--hom-alt`: Specifies input measurement or parameter value.
- `--p-values`: Specifies input measurement or parameter value.
- `--alpha`: Specifies input measurement or parameter value.
- `--method`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `id` | Parameter / observation metric | Required |
| `value` | Parameter / observation metric | Required |
| `qty` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t hardy-weinberg-tester .
docker run -p 8000:8000 hardy-weinberg-tester
```
