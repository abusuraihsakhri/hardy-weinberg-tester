# Hardy Weinberg Tester

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

---

## What It Does

Hardy-Weinberg Equilibrium Tester for genotype frequency analysis.

Tests genotype frequencies for HWE using chi-square and exact tests,
calculates inbreeding coefficients, provides deviation diagnostics,
and supports multiple testing correction.

Core HWE calculations use Python stdlib only (no external dependencies).

---

## Key Capabilities

### Analytical Functions

- **`allele_frequencies()`**: Calculate allele frequencies from genotype counts.
- **`expected_genotype_frequencies()`**: Calculate expected Hardy-Weinberg genotype frequencies (p², 2pq, q²).
- **`expected_genotype_counts()`**: Calculate expected genotype counts under HWE.
- **`chi_square_test()`**: Chi-square goodness-of-fit test for HWE.
- **`exact_test()`**: Exact test for small samples (Fisher's exact test equivalent).
- **`inbreeding_coefficient()`**: Calculate F = 1 - (Ho / He).
- **`diagnose_deviation()`**: Diagnose potential reasons for HWE deviation.
- **`bonferroni_correction()`**: Bonferroni multiple testing correction.
- **`benjamini_hochberg_correction()`**: Benjamini-Hochberg FDR correction.
- **`full_hwe_analysis()`**: Complete HWE analysis combining all methods.
- **`process_batch()`**: Batch process genotype data from CSV.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/hardy-weinberg-tester.git
cd hardy-weinberg-tester

# Install with web dependencies (optional)
pip install -e ".[web]"

# Install with test dependencies
pip install -e ".[test]"
```

### Environment Variables

The HMAC-SHA256 audit trail requires a secret key:

```bash
# Linux/macOS
export AUDIT_SECRET_KEY="your-secure-random-key"

# Windows
set AUDIT_SECRET_KEY=your-secure-random-key
```

**Important:** Never hardcode the audit secret key in source code. Use environment variables or a secrets manager.

---

## CLI Usage

### Chi-square HWE Test
```bash
python cli.py test --hom-ref 25 --het 50 --hom-alt 25
```

### Exact Test (Small Samples)
```bash
python cli.py exact --hom-ref 2 --het 1 --hom-alt 2
```

### Inbreeding Coefficient
```bash
python cli.py f --hom-ref 40 --het 20 --hom-alt 40
```

### Deviation Diagnosis
```bash
python cli.py diagnose --hom-ref 40 --het 20 --hom-alt 40
```

### Multiple Testing Correction
```bash
python cli.py correct --p-values 0.01,0.03,0.04 --method bonferroni
python cli.py correct --p-values 0.01,0.03,0.04 --method bh
```

### Full Analysis
```bash
python cli.py analyze --hom-ref 25 --het 50 --hom-alt 25
```

### Batch Processing from CSV
```bash
python cli.py batch -i input.csv -o results.csv
```

### Audit Commands (Enterprise Features)
```bash
python cli.py audit --task-id TASK-001 --primary-metric 12.0
python cli.py chat "What does this result mean?"
python cli.py verify-audit
```

---

## CSV Input Format

The batch processor accepts CSV files with genotype counts:

```csv
hom_ref,het,hom_alt
25,50,25
40,20,40
```

Alternative column names are also accepted: `AA`, `Aa`, `aa` or `homozygous_ref`, `heterozygous`, `homozygous_alt`.

---

## Security Features

- **Zero-PHI Outbound Interceptor:** AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
- **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation.

---

## Testing

Run the automated test suite:

```bash
pytest -v
```

All 38+ tests cover core HWE calculations, statistical tests, multiple testing correction, batch processing, PHI guard enforcement, and audit trail integrity.

---

## Container Deployment

```bash
docker build -t hardy-weinberg-tester .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret hardy-weinberg-tester
```

Or using Docker Compose:

```bash
docker-compose up
```

---

## Project Structure

```
hardy-weinberg-tester/
├── hwe.py                  # Core HWE calculations and statistical tests
├── cli.py                  # Command-line interface
├── test_hwe.py             # Core module tests
├── conftest.py             # Pytest configuration
├── agents/                 # Enterprise agent modules
│   ├── base.py             # PHI Guard, HMAC Audit Trail
│   ├── models.py           # Pydantic data models
│   ├── supervisor.py       # Multi-agent orchestrator
│   ├── workers.py          # Specialized worker agents
│   ├── api.py              # FastAPI REST endpoints
│   ├── metrics.py          # Prometheus metrics
│   ├── learning.py         # Bayesian calibration engine
│   ├── streamer.py         # WebSocket telemetry
│   └── llm_factory.py      # LLM provider abstraction
├── enrichment.py           # Domain enrichment feature suite
├── simulator.py            # High-throughput simulation
├── web/index.html          # Operations dashboard
├── tests/                  # Additional test suites
├── pyproject.toml          # Project metadata and dependencies
├── Dockerfile              # Container build
├── docker-compose.yml      # Container orchestration
├── openapi_spec.json       # OpenAPI 3.1 specification
└── benchmark_dataset.json  # Golden benchmark test cases
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
