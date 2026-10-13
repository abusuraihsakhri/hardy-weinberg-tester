"""
FastAPI REST API Server for Hardy Weinberg Tester.
"""
from typing import Dict, Any, List
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field, StrictInt
from hwe import full_hwe_analysis
from .base import AuditLogger, PHIGuard
from .models import SystemTaskPayload, ConsensusDossier
from .supervisor import SystemSupervisor

supervisor = SystemSupervisor(model_provider="mock")

app = FastAPI(
    title="Hardy Weinberg Tester API",
    description="Enterprise Distributed Component Platform (Clinical & Biomedical AI)",
    version="3.0.0-ENTERPRISE",
)


class ChatRequest(BaseModel):
    query: str


class HWECounts(BaseModel):
    hom_ref: StrictInt = Field(ge=0)
    het: StrictInt = Field(ge=0)
    hom_alt: StrictInt = Field(ge=0)


@app.get("/health")
def health():
    return {"status": "HEALTHY", "service": "hardy-weinberg-tester", "domain": "Clinical & Biomedical AI", "standard": "CAP / CLSI / ISO Standards", "version": "3.0.0-ENTERPRISE"}


@app.get("/metrics")
def metrics():
    return {
        "dossiers_processed_total": len(supervisor.dossier_registry),
        "audit_blocks_total": len(AuditLogger.get_trail()),
        "system_status": "NOMINAL_OPTIMAL"
    }


@app.post("/api/audit")
def api_audit(payload: SystemTaskPayload):
    dossier = supervisor.process_task(payload)
    return dossier.to_dict()


@app.post("/api/chat")
def api_chat(req: ChatRequest):
    try:
        ans = supervisor.query_supervisory_chat(req.query)
        return {"response": ans}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/audit/logs")
def api_audit_logs():
    return {"audit_trail": AuditLogger.get_trail(), "verified": AuditLogger.verify_integrity()}


@app.post("/api/hwe/analyze")
def api_hwe_analyze(counts: HWECounts):
    """Analyze biallelic diploid genotype counts, without storing inputs."""
    return full_hwe_analysis(counts.hom_ref, counts.het, counts.hom_alt)
