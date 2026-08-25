# Hardy-Weinberg Equilibrium Tester

Test genotype frequencies for Hardy-Weinberg Equilibrium using chi-square and exact tests, with deviation diagnostics and multiple testing correction.

## Features

- **Chi-Square Test**: χ² = Σ((observed - expected)² / expected), df=1
- **Exact Test**: Fisher's exact test equivalent for small samples
- **Inbreeding Coefficient**: F = 1 - (Ho/He)
- **Deviation Diagnoses**: Inbreeding, population stratification, selection, genotyping error
- **Multiple Testing**: Bonferroni and Benjamini-Hochberg FDR correction
- **Allele Frequencies**: p (ref), q (alt), expected genotype frequencies

## Installation

Zero dependencies — Python 3.8+ stdlib only.

```bash
cd hardy-weinberg-tester
```

## Usage

### CLI

```bash
# Chi-square HWE test
python cli.py test --hom-ref 25 --het 50 --hom-alt 25

# Exact test (small samples)
python cli.py exact --hom-ref 3 --het 2 --hom-alt 3

# Diagnose deviation
python cli.py diagnose --hom-ref 40 --het 20 --hom-alt 40

# Inbreeding coefficient
python cli.py f --hom-ref 25 --het 50 --hom-alt 25

# Multiple testing correction
python cli.py correct --p-values "0.01,0.03,0.04,0.5" --method bonferroni

# Full analysis
python cli.py analyze --hom-ref 25 --het 50 --hom-alt 25

# Batch process from CSV
python cli.py batch -i genotypes.csv -o results.csv
```

### Python API

```python
from hwe import (
    chi_square_test,
    exact_test,
    inbreeding_coefficient,
    diagnose_deviation,
    bonferroni_correction,
)

# Chi-square test
result = chi_square_test(25, 50, 25)
print(f"χ²={result['chi_square']}, p={result['p_value']}, HWE={result['in_hwe']}")

# Exact test
result = exact_test(3, 2, 3)
print(f"Exact p={result['p_value']}")

# Inbreeding coefficient
f = inbreeding_coefficient(40, 20, 40)
print(f"F={f['F']}, {f['interpretation']}")

# Diagnose
diag = diagnose_deviation(40, 20, 40)
for reason in diag['possible_reasons']:
    print(f"- {reason['reason']}: {reason['description']}")

# Bonferroni correction
corr = bonferroni_correction([0.01, 0.03, 0.04])
print(f"Adjusted α={corr['adjusted_alpha']}")
```

## Formulas

### Expected Genotype Frequencies
- p² (AA), 2pq (Aa), q² (aa) where p + q = 1

### Chi-Square Test
χ² = Σ((observed - expected)² / expected), df = 1

### Inbreeding Coefficient
F = 1 - (Ho / He)
- F > 0: excess homozygosity
- F < 0: excess heterozygosity
- F = 0: HWE

### Bonferroni Correction
Adjusted α = α / number_of_tests

## Common Reasons for HWE Deviation

| Reason | Direction | Mechanism |
|--------|-----------|-----------|
| Inbreeding | ↑ Homozygotes | Mating among relatives |
| Population stratification | ↑ Homozygotes | Mixed subpopulations (Wahlund effect) |
| Selection | ↑ or ↓ | Differential survival/reproduction |
| Genotyping error | ↑ Heterozygotes | Allele dropout, miscalling |
| Small sample size | Unreliable | Low statistical power |

## Tests

```bash
python -m pytest test_hwe.py -v
```

## License

MIT
