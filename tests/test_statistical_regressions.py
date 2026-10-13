"""Regression checks for conditional exact inference, BH and audit signatures."""
import math

import pytest

from hwe import exact_test, benjamini_hochberg_correction
from agents.base import AuditTrail


def reference_exact(aa, het, bb):
    """Independent integer-weight enumeration conditioned on allele totals."""
    n = aa + het + bb
    reference_alleles = 2 * aa + het
    probabilities = []
    observed_weight = None
    for ref_hom in range(n + 1):
        heterozygotes = reference_alleles - 2 * ref_hom
        alt_hom = n - ref_hom - heterozygotes
        if heterozygotes < 0 or alt_hom < 0:
            continue
        weight = math.factorial(n) * 2 ** heterozygotes // (
            math.factorial(ref_hom) *
            math.factorial(heterozygotes) *
            math.factorial(alt_hom)
        )
        probabilities.append(weight)
        if (ref_hom, heterozygotes, alt_hom) == (aa, het, bb):
            observed_weight = weight
    return sum(w for w in probabilities if w <= observed_weight) / sum(probabilities)


@pytest.mark.parametrize("aa,het,bb", [
    (2, 1, 2), (3, 2, 3), (5, 10, 5), (10, 0, 10),
    (1, 9, 0), (0, 0, 8), (7, 1, 0),
])
def test_exact_matches_integer_reference(aa, het, bb):
    assert exact_test(aa, het, bb)["p_value"] == pytest.approx(
        reference_exact(aa, het, bb), abs=1e-10
    )


def test_exact_large_sample_does_not_underflow_to_one():
    observed = exact_test(2500, 0, 2500)["p_value"]
    assert 0 <= observed < 1e-50


def test_bh_step_up_includes_earlier_nonpassing_rank():
    result = benjamini_hochberg_correction([0.02, 0.029, 0.049])
    assert result["num_significant"] == 3
    assert all(r["significant"] for r in result["results"])
    assert [r["adjusted_p_value"] for r in result["results"]] == pytest.approx(
        [0.03, 0.0435, 0.049]
    )


def test_bh_preserves_input_order_and_handles_ties():
    result = benjamini_hochberg_correction([0.6, 0.01, 0.01, 0.1])
    assert [r["test"] for r in result["results"]] == [1, 2, 3, 4]
    assert [r["adjusted_p_value"] for r in result["results"]] == pytest.approx(
        [0.6, 0.02, 0.02, 0.13333333333333333]
    )


def test_audit_recomputes_entry_signature():
    audit = AuditTrail(secret_key="unit-test-secret")
    audit.log("test", "operator", "CREATED", {"status": "ok"})
    audit.log("test", "operator", "CHECKED", {"status": "ok"})
    assert audit.verify_integrity()
    snapshot = audit.get_trail()
    snapshot[0]["actor"] = "other"
    assert audit.verify_integrity()  # snapshots cannot edit the live ledger
    audit.logs[0]["actor"] = "other"
    assert not audit.verify_integrity()


def test_audit_detects_last_entry_tampering():
    audit = AuditTrail(secret_key="unit-test-secret")
    audit.log("test", "operator", "CREATED", {"status": "ok"})
    audit.logs[0]["payload_hash"] = "forged"
    assert not audit.verify_integrity()
