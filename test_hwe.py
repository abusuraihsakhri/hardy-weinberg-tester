#!/usr/bin/env python3
"""Tests for Hardy-Weinberg Equilibrium Tester — 20+ real functional tests."""
import pytest
from hwe import (
    allele_frequencies,
    expected_genotype_frequencies,
    expected_genotype_counts,
    chi_square_test,
    exact_test,
    inbreeding_coefficient,
    diagnose_deviation,
    bonferroni_correction,
    benjamini_hochberg_correction,
    full_hwe_analysis,
    process_batch,
)


# ---------------------------------------------------------------------------
# Allele Frequency Tests
# ---------------------------------------------------------------------------

class TestAlleleFrequencies:
    def test_equal_freq(self):
        result = allele_frequencies(25, 50, 25)
        assert abs(result['p'] - 0.5) < 1e-10
        assert abs(result['q'] - 0.5) < 1e-10
        assert result['n'] == 100

    def test_rare_alt(self):
        result = allele_frequencies(90, 9, 1)
        assert result['p'] > 0.9
        assert result['q'] < 0.1

    def test_all_homozygous_ref(self):
        result = allele_frequencies(50, 0, 0)
        assert result['p'] == 1.0
        assert result['q'] == 0.0

    def test_zero_samples(self):
        result = allele_frequencies(0, 0, 0)
        assert result['p'] == 0.0
        assert result['n'] == 0


# ---------------------------------------------------------------------------
# Expected Genotype Frequency Tests
# ---------------------------------------------------------------------------

class TestExpectedFrequencies:
    def test_frequencies_sum_to_one(self):
        hw = expected_genotype_frequencies(0.5, 0.5)
        total = hw['p_squared'] + hw['two_pq'] + hw['q_squared']
        assert abs(total - 1.0) < 1e-10

    def test_equal_alleles(self):
        hw = expected_genotype_frequencies(0.5, 0.5)
        assert abs(hw['p_squared'] - 0.25) < 1e-10
        assert abs(hw['two_pq'] - 0.5) < 1e-10
        assert abs(hw['q_squared'] - 0.25) < 1e-10

    def test_rare_allele(self):
        hw = expected_genotype_frequencies(0.99, 0.01)
        assert hw['q_squared'] < 0.001  # Very rare homozygous alt

    def test_expected_counts(self):
        result = expected_genotype_counts(25, 50, 25)
        assert abs(result['expected_AA'] - 25.0) < 0.1
        assert abs(result['expected_Aa'] - 50.0) < 0.1
        assert abs(result['expected_aa'] - 25.0) < 0.1


# ---------------------------------------------------------------------------
# Chi-Square Test Tests
# ---------------------------------------------------------------------------

class TestChiSquareTest:
    def test_perfect_hwe(self):
        result = chi_square_test(25, 50, 25)
        assert result['chi_square'] < 0.01
        assert result['in_hwe'] is True

    def test_deviation(self):
        # Excess homozygosity
        result = chi_square_test(40, 20, 40)
        assert result['chi_square'] > 10.0
        assert result['in_hwe'] is False

    def test_mild_deviation(self):
        result = chi_square_test(30, 40, 30)
        # This is actually close to HWE for p=0.5
        assert isinstance(result['chi_square'], float)

    def test_zero_samples(self):
        result = chi_square_test(0, 0, 0)
        assert 'error' in result

    def test_p_value_range(self):
        result = chi_square_test(25, 50, 25)
        assert 0.0 <= result['p_value'] <= 1.0

    def test_df_is_one(self):
        result = chi_square_test(25, 50, 25)
        assert result['df'] == 1


# ---------------------------------------------------------------------------
# Exact Test Tests
# ---------------------------------------------------------------------------

class TestExactTest:
    def test_perfect_hwe(self):
        result = exact_test(25, 50, 25)
        assert result['p_value'] > 0.05
        assert result['in_hwe'] is True

    def test_small_sample(self):
        result = exact_test(2, 1, 2)
        assert 'p_value' in result
        assert result['method'] == 'exact'

    def test_extreme_deviation(self):
        result = exact_test(10, 0, 10)
        assert result['p_value'] < 0.05

    def test_zero_samples(self):
        result = exact_test(0, 0, 0)
        assert 'error' in result


# ---------------------------------------------------------------------------
# Inbreeding Coefficient Tests
# ---------------------------------------------------------------------------

class TestInbreedingCoefficient:
    def test_no_inbreeding(self):
        result = inbreeding_coefficient(25, 50, 25)
        assert abs(result['F']) < 0.1

    def test_inbreeding(self):
        result = inbreeding_coefficient(40, 20, 40)
        assert result['F'] > 0.1
        assert 'homozygosity' in result['interpretation'].lower() or 'inbreeding' in result['interpretation'].lower()

    def test_excess_het(self):
        result = inbreeding_coefficient(10, 80, 10)
        assert result['F'] < -0.1

    def test_monomorphic(self):
        result = inbreeding_coefficient(50, 0, 0)
        assert result['F'] == 0.0


# ---------------------------------------------------------------------------
# Deviation Diagnosis Tests
# ---------------------------------------------------------------------------

class TestDiagnosis:
    def test_hwe_no_reasons(self):
        result = diagnose_deviation(25, 50, 25)
        assert result['in_hwe'] is True
        assert len(result['possible_reasons']) == 0

    def test_inbreeding_diagnosis(self):
        result = diagnose_deviation(40, 20, 40)
        assert result['in_hwe'] is False
        assert len(result['possible_reasons']) > 0

    def test_small_sample_warning(self):
        result = diagnose_deviation(3, 2, 3)
        assert any('small' in w.lower() for w in result['warnings'])

    def test_rare_allele_warning(self):
        # (99, 1, 0): q = 1/200 = 0.005 < 0.01 -> rare allele warning
        result = diagnose_deviation(99, 1, 0)
        assert any('rare' in w.lower() or 'maf' in w.lower() for w in result['warnings'])


# ---------------------------------------------------------------------------
# Bonferroni Correction Tests
# ---------------------------------------------------------------------------

class TestBonferroni:
    def test_single_test(self):
        result = bonferroni_correction([0.03])
        assert result['adjusted_alpha'] == 0.05
        assert result['results'][0]['significant'] is True

    def test_multiple_tests(self):
        result = bonferroni_correction([0.01, 0.03, 0.04], alpha=0.05)
        assert result['adjusted_alpha'] == pytest.approx(0.05 / 3)
        assert result['results'][0]['significant'] is True  # 0.01 < 0.0167
        assert result['results'][1]['significant'] is False  # 0.03 > 0.0167

    def test_empty(self):
        result = bonferroni_correction([])
        assert result['num_tests'] == 0


# ---------------------------------------------------------------------------
# Benjamini-Hochberg Tests
# ---------------------------------------------------------------------------

class TestBenjaminiHochberg:
    def test_basic(self):
        result = benjamini_hochberg_correction([0.01, 0.04, 0.03], alpha=0.05)
        assert result['num_tests'] == 3
        assert 'num_significant' in result

    def test_empty(self):
        result = benjamini_hochberg_correction([])
        assert result['num_tests'] == 0


# ---------------------------------------------------------------------------
# Full Analysis Tests
# ---------------------------------------------------------------------------

class TestFullAnalysis:
    def test_full_analysis(self):
        result = full_hwe_analysis(25, 50, 25)
        assert 'chi_square_test' in result
        assert 'exact_test' in result
        assert 'inbreeding_coefficient' in result
        assert 'diagnosis' in result
        assert 'allele_frequencies' in result


# ---------------------------------------------------------------------------
# Batch Processing Tests
# ---------------------------------------------------------------------------

class TestBatchProcessing:
    def test_batch_csv(self, tmp_path):
        csv_in = tmp_path / 'in.csv'
        csv_out = tmp_path / 'out.csv'
        csv_in.write_text(
            'hom_ref,het,hom_alt\n'
            '25,50,25\n'
            '40,20,40\n',
            encoding='utf-8'
        )
        results = process_batch(str(csv_in), str(csv_out))
        assert len(results) == 2
        assert csv_out.exists()
        assert 'p_value' in results[0]
        assert 'in_hwe' in results[0]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
