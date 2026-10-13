"""End-to-end smoke test of the browser Python runtime and CSV workflow.

Run: python -m pip install playwright && python -m playwright install chromium
     python tests/browser_smoke.py
Requires access to the pinned Pyodide CDN.
"""
import subprocess
import sys
import time
from pathlib import Path
from urllib.request import urlopen

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent.parent
PORT = 8765
URL = f"http://127.0.0.1:{PORT}/web/"


def main():
    server = subprocess.Popen(
        [sys.executable, "-m", "http.server", str(PORT), "--bind", "127.0.0.1"],
        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    try:
        for attempt in range(40):
            try:
                with urlopen(URL, timeout=1) as response:
                    assert response.status == 200
                break
            except OSError:
                time.sleep(0.25)
        else:
            raise AssertionError("Local server did not start")
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            try:
                page = browser.new_page(accept_downloads=True)
                errors = []
                page.on("pageerror", lambda exception: errors.append(str(exception)))
                page.goto(URL, wait_until="domcontentloaded")
                page.wait_for_selector("#runtime.ready", timeout=120000)
                assert page.locator("#p").inner_text() == "0.5"
                assert page.locator("#q").inner_text() == "0.5"
                assert float(page.locator("#exact").inner_text()) > 0.05

                page.locator("#example").click()
                assert page.locator("#aa").input_value() == "40"
                assert float(page.locator("#exact").inner_text()) < 0.05
                assert "departure" in page.locator("#interpretation").inner_text()

                page.locator("#upload").set_input_files({
                    "name": "genotypes.csv",
                    "mimeType": "text/csv",
                    "buffer": b"hom_ref,het,hom_alt\n25,50,25\n40,20,40\n",
                })
                page.locator("#process").click()
                page.get_by_text("2 loci analyzed", exact=False).wait_for(timeout=30000)
                with page.expect_download(timeout=30000) as download_event:
                    page.locator("#download").click()
                text = Path(download_event.value.path()).read_text(encoding="utf-8")
                assert "bonferroni_significant" in text
                assert "25,50,25" in text
                assert not errors, errors
                print("Browser HWE calculation and CSV export passed")
            finally:
                browser.close()
    finally:
        server.terminate()
        try:
            server.wait(timeout=5)
        except subprocess.TimeoutExpired:
            server.kill()


if __name__ == "__main__":
    main()
