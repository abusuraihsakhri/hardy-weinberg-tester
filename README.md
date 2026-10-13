# Hardy–Weinberg Tester

A Python toolkit for evaluating Hardy–Weinberg equilibrium (HWE) in diploid biallelic genotype counts. It provides a command-line interface, batch CSV analysis, a browser calculator, and an optional FastAPI service.

## Capabilities

- Reference and alternate allele frequencies, expected genotype counts and observed/expected heterozygosity
- Pearson chi-square test with one degree of freedom
- Two-sided conditional exact HWE test using numerically stable log probabilities
- Inbreeding coefficient `F = 1 − Ho / He` and descriptive deviation warnings
- Bonferroni and Benjamini–Hochberg multiple-testing corrections
- CSV batch processing, with exported chi-square statistics, F and Bonferroni results
- Browser calculations using the **same Python engine** via Pyodide; no backend needed
- Optional audit/supervisor CLI and API features, separate from the statistical methods

A significant HWE test does not identify a specific biological mechanism. A nonsignificant p-value does not demonstrate that equilibrium holds. The chi-square approximation can be unreliable with sparse expected counts; prefer the exact test in that setting. The tool assumes independent autosomal diploid observations and does not account for relatedness, population stratification, sex-linked inheritance or genotype uncertainty.

## Browser application

The browser application is served from `web/index.html`. For local use, run from the repository root:

```bash
python -m http.server 8000
```

Then open **http://localhost:8000/web/**. Enter AA, Aa and aa genotype counts to view allele frequencies, exact and chi-square p-values, expected counts, and F. Alternatively, choose a CSV with headers `hom_ref,het,hom_alt` or `AA,Aa,aa` and download results.

The first visit requires internet access to load Pyodide 0.28.2 from jsDelivr. The browser then fetches this repository's `hwe.py` and runs the statistics locally. Input counts and CSV data are not uploaded by the browser application. Requests to the CDN reveal ordinary request metadata to that provider. Very large batches may take longer to process.

GitHub Pages publishing is defined in `.github/workflows/pages.yml`. The repository's Pages setting must use **GitHub Actions** for publishing. Do not assume a public site is available until the deployment completes and its URL has been verified.

## Python installation

Requires Python 3.10 or newer.

```bash
git clone https://github.com/abusuraihsakhri/hardy-weinberg-tester.git
cd hardy-weinberg-tester
python -m pip install -e ".[test]"
```

Core mathematical functions in `hwe.py` use only the standard library. The installed CLI imports Pydantic through the optional supervisor code.

## CLI examples

```bash
hwe-tester test --hom-ref 25 --het 50 --hom-alt 25
hwe-tester exact --hom-ref 2 --het 1 --hom-alt 2
hwe-tester f --hom-ref 40 --het 20 --hom-alt 40
hwe-tester analyze --hom-ref 25 --het 50 --hom-alt 25
hwe-tester diagnose --hom-ref 40 --het 20 --hom-alt 40
hwe-tester correct --p-values 0.02,0.029,0.049 --method bh
hwe-tester batch -i sample.csv -o results.csv
```

The batch processor applies chi-square p-values and Bonferroni correction; its adjusted significance flags should not be used without inspecting sparse-count warnings and, when necessary, performing exact tests separately.

## Optional API and audit operations

```bash
python -m pip install -e ".[web]"
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_urlsafe(32))')"
uvicorn agents.api:app --host 127.0.0.1 --port 8000
```

Submit counts to `POST /api/hwe/analyze`, for example:

```bash
curl -sS -X POST http://127.0.0.1:8000/api/hwe/analyze \
  -H 'Content-Type: application/json' \
  -d '{"hom_ref":25,"het":50,"hom_alt":25}'
```

The same service exposes `/health`, `/metrics`, `/api/audit`, `/api/audit/logs` and `/api/chat`. The supervisor/chat provider in this project is a deterministic **mock**, not a validated clinical system or live language-model integration. The PHI guard uses pattern matching, which cannot guarantee de-identification, and the HMAC trail is stored **in process memory**, not in durable storage. Do not expose these service endpoints on an untrusted network without authentication, transport security, and appropriate safeguards. Do not submit patient-identifiable information.

For the container:

```bash
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_urlsafe(32))')"
docker compose up --build
```

The container runs the FastAPI service on port 8000. Docker Compose requires `AUDIT_SECRET_KEY` from the host environment.

## Testing

```bash
python -m pip install -e ".[web,test]"
python -m pytest -q
python -m compileall -q hwe.py cli.py agents tests
node --check web/app.js
python -m pip check
```

GitHub Actions runs the Python tests across Python 3.10–3.13. Regression tests independently check exact p-values, Benjamini–Hochberg step-up decisions, audit tampering detection and API input validation.

## Technology, compatibility and license

Python standard library, Pydantic v2, optional FastAPI/Uvicorn, and Pyodide 0.28.2 for browser Python. The browser app requires a current WebAssembly-enabled browser (Chrome, Firefox, Edge or Safari) and network access to the Pyodide CDN on first load. There is no JavaScript bundler or Node application dependency.

Licensed under the [MIT License](LICENSE).
