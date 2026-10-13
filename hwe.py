#!/usr/bin/env python3
"""
Hardy-Weinberg Equilibrium Tester — Real Implementation

Tests genotype frequencies for HWE using chi-square and exact tests,
calculates inbreeding coefficients, provides deviation diagnostics,
and supports multiple testing correction.

Zero external dependencies — Python stdlib only.
"""

import csv
import math
import sys
from typing import Dict, Any, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Validation Helpers
# ---------------------------------------------------------------------------

def _validate_genotype_counts(hom_ref: int, het: int, hom_alt: int) -> None:
    """Validate genotype counts are non-negative integers.

    Raises:
        TypeError: If any count is not an integer.
        ValueError: If any count is negative.
    """
    for name, val in [('hom_ref', hom_ref), ('het', het), ('hom_alt', hom_alt)]:
        if not isinstance(val, int):
            raise TypeError(f"{name} must be an integer, got {type(val).__name__}")
        if val < 0:
            raise ValueError(f"{name} must be non-negative, got {val}")


# ---------------------------------------------------------------------------
# Core HWE Calculations
# ---------------------------------------------------------------------------

def allele_frequencies(hom_ref: int, het: int, hom_alt: int) -> Dict[str, float]:
    """Calculate allele frequencies from genotype counts.

    Args:
        hom_ref: Count of homozygous reference (AA) genotypes.
        het: Count of heterozygous (Aa) genotypes.
        hom_alt: Count of homozygous alternate (aa) genotypes.

    Returns:
        Dict with p (ref freq), q (alt freq), n (sample size).

    Raises:
        TypeError: If any count is not an integer.
        ValueError: If any count is negative.
    """
    _validate_genotype_counts(hom_ref, het, hom_alt)
    n = hom_ref + het + hom_alt
    if n == 0:
        return {'p': 0.0, 'q': 0.0, 'n': 0}

    total_alleles = 2 * n
    ref_alleles = 2 * hom_ref + het
    alt_alleles = 2 * hom_alt + het

    p = ref_alleles / total_alleles  # Reference allele frequency
    q = alt_alleles / total_alleles  # Alternate allele frequency

    return {'p': p, 'q': q, 'n': n}


def expected_genotype_frequencies(p: float, q: float) -> Dict[str, float]:
    """Calculate expected Hardy-Weinberg genotype frequencies.

    p² (homozygous reference / AA)
    2pq (heterozygous / Aa)
    q² (homozygous alternate / aa)
    """
    return {
        'p_squared': p * p,      # Expected AA frequency
        'two_pq': 2.0 * p * q,   # Expected Aa frequency
        'q_squared': q * q,       # Expected aa frequency
    }


def expected_genotype_counts(hom_ref: int, het: int, hom_alt: int) -> Dict[str, float]:
    """Calculate expected genotype counts under HWE."""
    freqs = allele_frequencies(hom_ref, het, hom_alt)
    n = freqs['n']
    p, q = freqs['p'], freqs['q']
    hw = expected_genotype_frequencies(p, q)

    return {
        'expected_AA': hw['p_squared'] * n,
        'expected_Aa': hw['two_pq'] * n,
        'expected_aa': hw['q_squared'] * n,
        'observed_AA': hom_ref,
        'observed_Aa': het,
        'observed_aa': hom_alt,
        'n': n,
        'p': p,
        'q': q,
    }


# ---------------------------------------------------------------------------
# Chi-Square Test
# ---------------------------------------------------------------------------

def chi_square_test(hom_ref: int, het: int, hom_alt: int) -> Dict[str, Any]:
    """Chi-square goodness-of-fit test for Hardy-Weinberg Equilibrium.

    χ² = Σ((observed - expected)² / expected)

    Returns:
        Dict with chi_square, p_value, df, in_hwe, expected counts, etc.

    Raises:
        TypeError: If any count is not an integer.
        ValueError: If any count is negative.
    """
    _validate_genotype_counts(hom_ref, het, hom_alt)
    n = hom_ref + het + hom_alt
    if n == 0:
        return {'chi_square': 0.0, 'p_value': 1.0, 'df': 1, 'in_hwe': True,
                'error': 'No genotypes observed'}

    freqs = allele_frequencies(hom_ref, het, hom_alt)
    p, q = freqs['p'], freqs['q']

    exp = expected_genotype_frequencies(p, q)
    exp_AA = exp['p_squared'] * n
    exp_Aa = exp['two_pq'] * n
    exp_aa = exp['q_squared'] * n

    # Chi-square calculation
    chi_sq = 0.0
    for obs, expected in [(hom_ref, exp_AA), (het, exp_Aa), (hom_alt, exp_aa)]:
        if expected > 0:
            chi_sq += (obs - expected) ** 2 / expected

    df = 1  # 3 categories - 1 (constraint: allele freq) - 1 (df for estimation)
    p_value = _chi_square_p_value(chi_sq, df)

    return {
        'chi_square': round(chi_sq, 4),
        'p_value': p_value,
        'df': df,
        'in_hwe': p_value > 0.05,
        'alpha': 0.05,
        'observed': {'AA': hom_ref, 'Aa': het, 'aa': hom_alt},
        'expected': {'AA': round(exp_AA, 2), 'Aa': round(exp_Aa, 2), 'aa': round(exp_aa, 2)},
        'allele_freq': {'p': round(p, 6), 'q': round(q, 6)},
        'sample_size': n,
    }


def _chi_square_p_value(chi_sq: float, df: int) -> float:
    """Approximate p-value for chi-square distribution (stdlib only)."""
    if chi_sq <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(chi_sq / 2.0))
    elif df == 2:
        return math.exp(-chi_sq / 2.0)
    else:
        # Wilson-Hilferty approximation for general df
        a = df / 2.0
        z = chi_sq / 2.0
        if z < a + 1:
            return 1.0 - _gammainc_series(a, z)
        else:
            return _gammainc_cf(a, z)


def _gammainc_series(a: float, x: float) -> float:
    """Series expansion for regularized lower incomplete gamma function."""
    if x == 0:
        return 0.0
    ap = a
    total = 1.0 / a
    term = 1.0 / a
    for _ in range(300):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * 1e-15:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gammainc_cf(a: float, x: float) -> float:
    """Continued fraction for regularized upper incomplete gamma function."""
    b = x + 1.0 - a
    c = 1e30
    d = 1.0 / b
    h = d
    for i in range(1, 300):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-30:
            d = 1e-30
        c = b + an / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


# ---------------------------------------------------------------------------
# Exact Test (Fisher's Exact Test Equivalent)
# ---------------------------------------------------------------------------

def exact_test(hom_ref: int, het: int, hom_alt: int) -> Dict[str, Any]:
    """Two-sided conditional exact HWE test (probability-ordering definition).

    Condition on the observed reference/alternate allele counts. The
    probability of a genotype configuration is proportional to
    n! * 2**n_het / (n_ref_hom! * n_het! * n_alt_hom!).
    Log weights prevent factorial overflow and underflow for large samples.
    """
    _validate_genotype_counts(hom_ref, het, hom_alt)
    n = hom_ref + het + hom_alt
    if n == 0:
        return {'p_value': 1.0, 'method': 'exact', 'error': 'No genotypes observed'}

    n_ref = 2 * hom_ref + het
    n_alt = 2 * hom_alt + het
    rare = min(n_ref, n_alt)
    common = max(n_ref, n_alt)

    # The conditional sample space has fixed allele counts and heterozygote
    # parity. All listed states have identical numbers of each allele.
    def log_weight(heterozygotes: int) -> float:
        hom_rare = (rare - heterozygotes) // 2
        hom_common = (common - heterozygotes) // 2
        return (math.lgamma(n + 1) - math.lgamma(hom_rare + 1)
                - math.lgamma(heterozygotes + 1)
                - math.lgamma(hom_common + 1)
                + heterozygotes * math.log(2.0))

    observed_log_weight = log_weight(het)
    log_weights = [
        log_weight(h) for h in range(rare % 2, rare + 1, 2)
    ]
    max_log_weight = max(log_weights)
    scaled_weights = [math.exp(log_w - max_log_weight) for log_w in log_weights]
    total_weight = math.fsum(scaled_weights)
    # Relative log-space tolerance handles equivalent states without
    # incorrectly admitting every state when absolute probabilities are tiny.
    tolerance = 1e-10
    tail_weight = math.fsum(
        weight for log_w, weight in zip(log_weights, scaled_weights)
        if log_w <= observed_log_weight + tolerance
    )
    p_value = min(1.0, max(0.0, tail_weight / total_weight))
    freqs = allele_frequencies(hom_ref, het, hom_alt)
    result = {
        'p_value': p_value,
        'method': 'exact',
        'in_hwe': p_value > 0.05,
        'alpha': 0.05,
        'observed': {'AA': hom_ref, 'Aa': het, 'aa': hom_alt},
        'allele_freq': {'p': round(freqs['p'], 6), 'q': round(freqs['q'], 6)},
        'sample_size': n,
    }
    if rare == 0:
        result['warning'] = 'Monomorphic locus: no HWE deviation can be assessed'
    return result


# ---------------------------------------------------------------------------
# Inbreeding Coefficient
# ---------------------------------------------------------------------------

def inbreeding_coefficient(hom_ref: int, het: int, hom_alt: int) -> Dict[str, Any]:
    """Calculate inbreeding coefficient F = 1 - (Ho / He).

    F > 0: excess homozygosity (inbreeding, Wahlund effect)
    F < 0: excess heterozygosity (outbreeding, selection)
    F = 0: consistent with HWE

    Raises:
        TypeError: If any count is not an integer.
        ValueError: If any count is negative.
    """
    _validate_genotype_counts(hom_ref, het, hom_alt)
    n = hom_ref + het + hom_alt
    if n == 0:
        return {'F': 0.0, 'interpretation': 'No data', 'n': 0}

    freqs = allele_frequencies(hom_ref, het, hom_alt)
    p, q = freqs['p'], freqs['q']

    # Observed heterozygosity
    ho = het / n

    # Expected heterozygosity under HWE
    he = 2 * p * q

    if he == 0:
        return {'F': 0.0, 'Ho': ho, 'He': 0.0,
                'interpretation': 'Monomorphic locus — no heterozygosity expected',
                'n': n}

    f = 1.0 - (ho / he)

    # Interpretation
    if f > 0.1:
        interp = 'Excess homozygosity — possible inbreeding, population structure, or Wahlund effect'
        deviation = 'homozygote_excess'
    elif f < -0.1:
        interp = 'Excess heterozygosity — possible outbreeding, selection against homozygotes, or genotyping error'
        deviation = 'heterozygote_excess'
    else:
        interp = 'Consistent with Hardy-Weinberg Equilibrium'
        deviation = 'none'

    return {
        'F': round(f, 4),
        'Ho': round(ho, 4),
        'He': round(he, 4),
        'observed_heterozygosity': round(ho, 4),
        'expected_heterozygosity': round(he, 4),
        'interpretation': interp,
        'deviation_type': deviation,
        'n': n,
    }


# ---------------------------------------------------------------------------
# Deviation Diagnostics
# ---------------------------------------------------------------------------

def diagnose_deviation(hom_ref: int, het: int, hom_alt: int) -> Dict[str, Any]:
    """Diagnose potential reasons for HWE deviation.

    Checks:
        1. Excess homozygosity (inbreeding, population stratification)
        2. Excess heterozygosity (outbreeding, selection, genotyping error)
        3. Small sample size warning
        4. Rare allele warning

    Raises:
        TypeError: If any count is not an integer.
        ValueError: If any count is negative.
    """
    _validate_genotype_counts(hom_ref, het, hom_alt)
    n = hom_ref + het + hom_alt
    freqs = allele_frequencies(hom_ref, het, hom_alt)
    p, q = freqs['p'], freqs['q']

    f_result = inbreeding_coefficient(hom_ref, het, hom_alt)
    chi_result = chi_square_test(hom_ref, het, hom_alt)

    reasons = []
    warnings = []

    # Check for deviation
    if not chi_result['in_hwe']:
        if f_result['F'] > 0.1:
            reasons.append({
                'reason': 'Inbreeding',
                'description': 'Excess homozygosity suggests consanguinity or population substructure',
                'evidence': f"F = {f_result['F']:.4f}",
            })
            reasons.append({
                'reason': 'Population stratification (Wahlund effect)',
                'description': 'Mixing subpopulations with different allele frequencies creates apparent homozygote excess',
                'evidence': 'Consider testing subpopulations separately',
            })
        elif f_result['F'] < -0.1:
            reasons.append({
                'reason': 'Selection against homozygotes',
                'description': 'Heterozygote advantage (overdominance) can maintain excess heterozygosity',
                'evidence': f"F = {f_result['F']:.4f}",
            })
            reasons.append({
                'reason': 'Genotyping error',
                'description': 'Systematic miscalling can inflate heterozygote counts',
                'evidence': 'Check genotyping quality metrics',
            })

    # Small sample size warning
    if n < 25:
        warnings.append(f'Small sample size (n={n}): chi-square test may be unreliable. '
                        f'Use exact test instead.')

    # Rare allele warning
    if q < 0.01 or p < 0.01:
        warnings.append(f'Rare minor allele (MAF={min(p, q):.4f}): '
                        f'expected counts may be too small for reliable chi-square test.')

    # Check expected counts
    exp = expected_genotype_counts(hom_ref, het, hom_alt)
    if exp['expected_AA'] < 5 or exp['expected_Aa'] < 5 or exp['expected_aa'] < 5:
        warnings.append('Some expected genotype counts < 5: chi-square approximation may be invalid. '
                        'Consider using the exact test.')

    return {
        'in_hwe': chi_result['in_hwe'],
        'chi_square': chi_result['chi_square'],
        'p_value': chi_result['p_value'],
        'F': f_result['F'],
        'deviation_type': f_result['deviation_type'],
        'possible_reasons': reasons,
        'warnings': warnings,
        'recommendation': 'Use exact test' if n < 25 or min(p, q) < 0.01 else 'Chi-square test is appropriate',
    }


# ---------------------------------------------------------------------------
# Multiple Testing Correction
# ---------------------------------------------------------------------------

def bonferroni_correction(p_values: List[float],
                           alpha: float = 0.05) -> Dict[str, Any]:
    """Apply Bonferroni correction for multiple testing.

    Adjusted alpha = alpha / number_of_tests.
    A test is significant if p_value < adjusted_alpha.

    Raises:
        ValueError: If any p_value is not in [0, 1] or alpha is not in (0, 1].
    """
    if not 0 < alpha <= 1:
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    for i, p in enumerate(p_values):
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p_values[{i}] = {p} is outside valid range [0, 1]")
    m = len(p_values)
    if m == 0:
        return {'adjusted_alpha': alpha, 'results': [], 'num_tests': 0}

    adjusted_alpha = alpha / m

    results = []
    for i, p in enumerate(p_values):
        results.append({
            'test': i + 1,
            'p_value': p,
            'significant': p < adjusted_alpha,
            'adjusted_alpha': adjusted_alpha,
        })

    num_significant = sum(1 for r in results if r['significant'])

    return {
        'num_tests': m,
        'original_alpha': alpha,
        'adjusted_alpha': adjusted_alpha,
        'num_significant': num_significant,
        'results': results,
    }


def benjamini_hochberg_correction(p_values: List[float],
                                    alpha: float = 0.05) -> Dict[str, Any]:
    """Apply Benjamini-Hochberg FDR correction.

    Controls false discovery rate rather than family-wise error rate.
    Less conservative than Bonferroni.

    Raises:
        ValueError: If any p_value is not in [0, 1] or alpha is not in (0, 1].
    """
    if not 0 < alpha <= 1:
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    for i, p in enumerate(p_values):
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p_values[{i}] = {p} is outside valid range [0, 1]")
    m = len(p_values)
    if m == 0:
        return {'results': [], 'num_tests': 0}

    # Benjamini-Hochberg is a step-up test: once the largest passing rank
    # is found, all lower ranks must also be rejected, even if their own
    # rank-wise comparison failed.
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    ordered = [p for _, p in indexed]
    adjusted = [0.0] * m
    minimum = 1.0
    for i in range(m - 1, -1, -1):
        minimum = min(minimum, ordered[i] * m / (i + 1))
        adjusted[i] = minimum

    results = [None] * m
    for rank, ((original_idx, p), adjusted_p) in enumerate(
            zip(indexed, adjusted), start=1):
        results[original_idx] = {
            'test': original_idx + 1,
            'p_value': p,
            'rank': rank,
            'critical_value': rank / m * alpha,
            'adjusted_p_value': adjusted_p,
            'significant': adjusted_p <= alpha,
        }

    return {
        'num_tests': m,
        'alpha': alpha,
        'num_significant': sum(row['significant'] for row in results),
        'results': results,
    }

# ---------------------------------------------------------------------------
# Full HWE Analysis
# ---------------------------------------------------------------------------

def full_hwe_analysis(hom_ref: int, het: int, hom_alt: int) -> Dict[str, Any]:
    """Complete HWE analysis combining all methods and diagnostics.

    Raises:
        TypeError: If any count is not an integer.
        ValueError: If any count is negative.
    """
    _validate_genotype_counts(hom_ref, het, hom_alt)
    chi = chi_square_test(hom_ref, het, hom_alt)
    exact = exact_test(hom_ref, het, hom_alt)
    f_coeff = inbreeding_coefficient(hom_ref, het, hom_alt)
    diagnosis = diagnose_deviation(hom_ref, het, hom_alt)
    freqs = allele_frequencies(hom_ref, het, hom_alt)
    hw_freqs = expected_genotype_frequencies(freqs['p'], freqs['q'])

    return {
        'allele_frequencies': freqs,
        'expected_frequencies': hw_freqs,
        'chi_square_test': chi,
        'exact_test': exact,
        'inbreeding_coefficient': f_coeff,
        'diagnosis': diagnosis,
    }


# ---------------------------------------------------------------------------
# Batch Processing
# ---------------------------------------------------------------------------

def _parse_genotype_value(value: Any, field_name: str, row_num: int) -> int:
    """Parse a genotype count from CSV, with clear error messages."""
    if value is None or (isinstance(value, str) and value.strip() == ''):
        raise ValueError(
            f"Row {row_num}: Missing value for '{field_name}'. "
            f"Expected columns: hom_ref, het, hom_alt (or AA, Aa, aa)."
        )
    try:
        parsed = int(value)
    except (ValueError, TypeError) as e:
        raise ValueError(
            f"Row {row_num}: Invalid value '{value}' for '{field_name}'. "
            f"Must be a non-negative integer."
        ) from e
    if parsed < 0:
        raise ValueError(
            f"Row {row_num}: Negative value {parsed} for '{field_name}'. "
            f"Genotype counts must be non-negative."
        )
    return parsed


def process_batch(input_csv: str, output_csv: str) -> List[Dict[str, Any]]:
    """Process a CSV of genotype counts and test each for HWE.

    Expected CSV columns: hom_ref (or AA), het (or Aa), hom_alt (or aa).
    Each row is tested for HWE and results are written to output_csv.

    Args:
        input_csv: Path to input CSV with genotype counts.
        output_csv: Path to output CSV for results.

    Returns:
        List of result dicts, one per input row.

    Raises:
        FileNotFoundError: If input_csv does not exist.
        ValueError: If CSV is missing required columns or has invalid data.
    """
    import os
    if not os.path.isfile(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    with open(input_csv, mode='r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"CSV file has no headers: {input_csv}")
        rows = list(reader)

    # Detect column naming convention
    has_standard = {'hom_ref', 'het', 'hom_alt'}.issubset(set(fieldnames))
    has_aa = {'AA', 'Aa', 'aa'}.issubset(set(fieldnames))
    if not has_standard and not has_aa:
        raise ValueError(
            f"CSV must contain columns (hom_ref, het, hom_alt) or (AA, Aa, aa). "
            f"Found: {fieldnames}"
        )

    out_fields = fieldnames + ['p_value', 'chi_square', 'F', 'in_hwe', 'deviation_type']
    out_rows = []
    p_values = []

    for idx, r in enumerate(rows, start=1):
        try:
            if has_standard:
                hom_ref = _parse_genotype_value(r.get('hom_ref'), 'hom_ref', idx)
                het = _parse_genotype_value(r.get('het'), 'het', idx)
                hom_alt = _parse_genotype_value(r.get('hom_alt'), 'hom_alt', idx)
            else:
                hom_ref = _parse_genotype_value(r.get('AA'), 'AA', idx)
                het = _parse_genotype_value(r.get('Aa'), 'Aa', idx)
                hom_alt = _parse_genotype_value(r.get('aa'), 'aa', idx)
        except ValueError:
            raise

        chi = chi_square_test(hom_ref, het, hom_alt)
        f_coeff = inbreeding_coefficient(hom_ref, het, hom_alt)

        row_dict = dict(r)
        row_dict['p_value'] = round(chi['p_value'], 6)
        row_dict['chi_square'] = chi['chi_square']
        row_dict['F'] = f_coeff['F']
        row_dict['in_hwe'] = chi['in_hwe']
        row_dict['deviation_type'] = f_coeff['deviation_type']
        out_rows.append(row_dict)
        p_values.append(chi['p_value'])

    # Apply Bonferroni correction
    if p_values:
        correction = bonferroni_correction(p_values)
        for i, row in enumerate(out_rows):
            row['bonferroni_significant'] = correction['results'][i]['significant']
            row['bonferroni_alpha'] = correction['adjusted_alpha']

    with open(output_csv, mode='w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fields + ['bonferroni_significant', 'bonferroni_alpha'])
        writer.writeheader()
        writer.writerows(out_rows)

    return out_rows
