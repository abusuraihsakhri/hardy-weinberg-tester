#!/usr/bin/env python3
"""
CLI for Hardy-Weinberg Equilibrium Tester.

Subcommands:
  test      — Chi-square HWE test
  exact     — Exact test for small samples
  diagnose  — Diagnose HWE deviation reasons
  f         — Inbreeding coefficient
  correct   — Multiple testing correction
  analyze   — Full HWE analysis
  batch     — Batch process from CSV
  audit     — Dispatch task payload across workers
  chat      — Air-gapped supervisory conversational assistant
  verify-audit — Retrieve and verify HMAC audit trail
"""
import argparse
import json
import sys

from hwe import (
    allele_frequencies,
    expected_genotype_frequencies,
    chi_square_test,
    exact_test,
    inbreeding_coefficient,
    diagnose_deviation,
    bonferroni_correction,
    benjamini_hochberg_correction,
    full_hwe_analysis,
    process_batch,
)

from agents.supervisor import SystemSupervisor
from agents.base import AuditLogger
from agents.models import SystemTaskPayload


def cmd_test(args):
    """Chi-square HWE test."""
    result = chi_square_test(args.hom_ref, args.het, args.hom_alt)
    print(json.dumps(result, indent=2))


def cmd_exact(args):
    """Exact test for HWE."""
    result = exact_test(args.hom_ref, args.het, args.hom_alt)
    print(json.dumps(result, indent=2))


def cmd_diagnose(args):
    """Diagnose HWE deviation."""
    result = diagnose_deviation(args.hom_ref, args.het, args.hom_alt)
    print(json.dumps(result, indent=2))


def cmd_f(args):
    """Inbreeding coefficient."""
    result = inbreeding_coefficient(args.hom_ref, args.het, args.hom_alt)
    print(json.dumps(result, indent=2))


def cmd_correct(args):
    """Multiple testing correction."""
    p_values = [float(x) for x in args.p_values.split(',')]
    if args.method == 'bonferroni':
        result = bonferroni_correction(p_values, args.alpha)
    else:
        result = benjamini_hochberg_correction(p_values, args.alpha)
    print(json.dumps(result, indent=2))


def cmd_analyze(args):
    """Full HWE analysis."""
    result = full_hwe_analysis(args.hom_ref, args.het, args.hom_alt)
    print(json.dumps(result, indent=2))


def cmd_batch(args):
    """Batch process from CSV."""
    results = process_batch(args.input, args.output)
    print(f"Processed {len(results)} loci -> {args.output}")


def cmd_audit(args):
    """Dispatch task payload across workers."""
    supervisor = SystemSupervisor(model_provider="mock")
    payload = SystemTaskPayload(
        task_id=args.task_id,
        target_identifier=args.target_id,
        primary_metric=args.primary_metric,
        secondary_metric=args.secondary_metric,
        status_descriptor=args.status_descriptor,
        is_critical_flag=args.is_critical,
    )
    dossier = supervisor.process_task(payload)
    print(json.dumps(dossier.to_dict(), indent=2, default=str))
    return 0


def cmd_chat(args):
    """Air-gapped supervisory conversational assistant."""
    supervisor = SystemSupervisor(model_provider="mock")
    response = supervisor.query_supervisory_chat(" ".join(args.query))
    print(json.dumps({"response": response}, indent=2))
    return 0


def cmd_verify_audit(args):
    """Retrieve and verify HMAC audit trail."""
    trail = AuditLogger.get_trail()
    verified = AuditLogger.verify_integrity()
    print(json.dumps({"audit_trail": trail, "verified": verified, "block_count": len(trail)}, indent=2, default=str))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='hwe-tester',
        description='Hardy-Weinberg Equilibrium Tester — Chi-square, exact test, inbreeding, diagnostics'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Common arguments for genotype counts
    def add_genotype_args(p):
        p.add_argument('--hom-ref', type=int, required=True, help='Homozygous reference (AA) count')
        p.add_argument('--het', type=int, required=True, help='Heterozygous (Aa) count')
        p.add_argument('--hom-alt', type=int, required=True, help='Homozygous alternate (aa) count')

    # test
    p_test = subparsers.add_parser('test', help='Chi-square HWE test')
    add_genotype_args(p_test)
    p_test.set_defaults(func=cmd_test)

    # exact
    p_exact = subparsers.add_parser('exact', help='Exact test (small samples)')
    add_genotype_args(p_exact)
    p_exact.set_defaults(func=cmd_exact)

    # diagnose
    p_diag = subparsers.add_parser('diagnose', help='Diagnose HWE deviation')
    add_genotype_args(p_diag)
    p_diag.set_defaults(func=cmd_diagnose)

    # f
    p_f = subparsers.add_parser('f', help='Inbreeding coefficient')
    add_genotype_args(p_f)
    p_f.set_defaults(func=cmd_f)

    # correct
    p_corr = subparsers.add_parser('correct', help='Multiple testing correction')
    p_corr.add_argument('--p-values', required=True, help='Comma-separated p-values')
    p_corr.add_argument('--alpha', type=float, default=0.05, help='Significance level')
    p_corr.add_argument('--method', default='bonferroni', choices=['bonferroni', 'bh'],
                        help='Correction method')
    p_corr.set_defaults(func=cmd_correct)

    # analyze
    p_ana = subparsers.add_parser('analyze', help='Full HWE analysis')
    add_genotype_args(p_ana)
    p_ana.set_defaults(func=cmd_analyze)

    # batch
    p_batch = subparsers.add_parser('batch', help='Batch process from CSV')
    p_batch.add_argument('-i', '--input', required=True, help='Input CSV')
    p_batch.add_argument('-o', '--output', default='results.csv', help='Output CSV')
    p_batch.set_defaults(func=cmd_batch)

    # audit
    p_audit = subparsers.add_parser('audit', help='Dispatch task payload across workers')
    p_audit.add_argument('--task-id', required=True, help='Task identifier')
    p_audit.add_argument('--target-id', default='SPECIMEN-001', help='Target identifier')
    p_audit.add_argument('--primary-metric', type=float, default=10.0, help='Primary metric value')
    p_audit.add_argument('--secondary-metric', type=float, default=5.0, help='Secondary metric value')
    p_audit.add_argument('--status-descriptor', default='NOMINAL', help='Status descriptor')
    p_audit.add_argument('--is-critical', action='store_true', help='Critical flag')
    p_audit.set_defaults(func=cmd_audit)

    # chat
    p_chat = subparsers.add_parser('chat', help='Air-gapped supervisory conversational assistant')
    p_chat.add_argument('query', nargs='+', help='Query text')
    p_chat.set_defaults(func=cmd_chat)

    # verify-audit
    p_verify = subparsers.add_parser('verify-audit', help='Retrieve and verify HMAC audit trail')
    p_verify.set_defaults(func=cmd_verify_audit)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    main()
