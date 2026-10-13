"""API smoke tests for HWE count validation and analysis."""
from fastapi.testclient import TestClient
from agents.api import app

client = TestClient(app)


def test_hwe_endpoint_returns_real_statistics():
    response = client.post("/api/hwe/analyze", json={
        "hom_ref": 25, "het": 50, "hom_alt": 25
    })
    assert response.status_code == 200
    data = response.json()
    assert data["allele_frequencies"]["p"] == 0.5
    assert data["chi_square_test"]["p_value"] == 1.0
    assert data["exact_test"]["p_value"] > 0.05


def test_hwe_endpoint_rejects_negative_or_noninteger_counts():
    for invalid in [-1, 0.5, "8", True]:
        response = client.post("/api/hwe/analyze", json={
            "hom_ref": invalid, "het": 50, "hom_alt": 25
        })
        assert response.status_code == 422
