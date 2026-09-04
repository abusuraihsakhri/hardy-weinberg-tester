"""Pytest configuration — sets test environment variables before imports."""
import os

# Set a deterministic test key for HMAC audit operations during testing.
# This key is used only for test isolation and is not a production secret.
os.environ.setdefault("AUDIT_SECRET_KEY", "test-only-hardy-weinberg-audit-key-do-not-use-in-production")
