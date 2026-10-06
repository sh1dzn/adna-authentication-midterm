#!/usr/bin/env python3
"""Retrieve run metadata and FASTQ files from the ENA (task 2).

Everything this project analyses is fetched by this script, so that the
repository holds the *provenance* of the data rather than the data itself. Each
download is verified against the MD5 the archive publishes, because a silently
truncated FASTQ is far more damaging in an ancient DNA project than a failed
one: short, damaged reads look abnormal by design, so corruption does not
announce itself downstream.

The archive is a live service and occasionally returns 5xx or drops a
connection mid-transfer. We therefore retry with exponential backoff, resume
partial transfers with an HTTP Range request, and treat an MD5 mismatch as a
reason to re-fetch once before giving up.

Usage
-----
    python scripts/fetch_ena.py metadata PRJEB29990 > results/ena_metadata.tsv
    python scripts/fetch_ena.py download ERR3457840 --out data/raw
    python scripts/fetch_ena.py provenance config/samples.tsv --out results/provenance.tsv
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# The metadata we need in order to reason about platform (task 3), library
# protocol (task 3) and archaeological provenance (task 2).
FIELDS = [
    "study_accession", "sample_accession", "run_accession", "experiment_accession",
    "sample_alias", "sample_title", "scientific_name", "tax_id",
    "instrument_platform", "instrument_model",
    "library_layout", "library_strategy", "library_source", "library_selection",
    "library_name", "library_construction_protocol",
    "read_count", "base_count", "fastq_bytes", "fastq_ftp", "fastq_md5",
    "collection_date", "country", "location", "first_public",
]

RETRIES = 5
BACKOFF = 3.0
CHUNK = 1 << 20


def _urlopen(url: str, headers: dict | None = None, timeout: int = 120):
    req = urllib.request.Request(url, headers=headers or {})
    return urllib.request.urlopen(req, timeout=timeout)


def _with_retries(fn, what: str):
    """Run fn(), retrying transient archive failures with exponential backoff."""
    last: Exception | None = None
    for attempt in range(1, RETRIES + 1):
        try:
            return fn()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
            # A 404 means the accession is wrong; retrying will not help.
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 404:
                raise
            last = exc
            wait = BACKOFF * attempt
            print(f"  attempt {attempt}/{RETRIES} failed for {what} ({exc}); "
                  f"retrying in {wait:.0f}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"gave up on {what} after {RETRIES} attempts: {last}")


def filereport(accession: str, fields: list[str] | None = None) -> list[dict]:
    """Query the ENA portal for every run belonging to an accession."""
    cols = fields or FIELDS
    url = (f"{ENA_PORTAL}?accession={accession}&result=read_run"
           f"&fields={','.join(cols)}&format=tsv")

    def go():
        with _urlopen(url) as resp:
            return resp.read().decode()

    text = _with_retries(go, f"metadata {accession}").strip()
    if not text:
        return []
    lines = text.splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def download_one(ftp_path: str, expected_md5: str, out_dir: Path) -> Path:
    """Download a single FASTQ, resuming and verifying it."""
    url = f"https://{ftp_path}"
    dest = out_dir / Path(ftp_path).name

    if dest.exists() and expected_md5:
        if md5sum(dest) == expected_md5:
            print(f"  have  {dest.name} (md5 ok)", file=sys.stderr)
            return dest
        print(f"  {dest.name} exists but md5 differs; re-fetching", file=sys.stderr)
        dest.unlink()

    for attempt in (1, 2):
        part = dest.with_suffix(dest.suffix + ".part")
        offset = part.stat().st_size if part.exists() else 0

        def go():
            headers = {"Range": f"bytes={offset}-"} if offset else {}
            mode = "ab" if offset else "wb"
            with _urlopen(url, headers=headers, timeout=300) as resp, open(part, mode) as fh:
                while chunk := resp.read(CHUNK):
                    fh.write(chunk)

        if offset:
            print(f"  resume {dest.name} from {offset / 1e6:.0f} MB", file=sys.stderr)
        else:
            print(f"  get   {dest.name}", file=sys.stderr)
        _with_retries(go, dest.name)

        part.rename(dest)
        if not expected_md5:
            print(f"  warn  {dest.name}: archive published no md5, cannot verify",
                  file=sys.stderr)
            return dest
        if md5sum(dest) == expected_md5:
            return dest
        print(f"  md5 mismatch on {dest.name} (attempt {attempt})", file=sys.stderr)
        dest.unlink()

    raise RuntimeError(f"{dest.name}: md5 never matched {expected_md5}")


def cmd_metadata(accessions: list[str]) -> None:
    rows = [r for acc in accessions for r in filereport(acc)]
    print("\t".join(FIELDS))
    for row in rows:
        print("\t".join(row.get(c, "") for c in FIELDS))


def cmd_download(accessions: list[str], out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for acc in accessions:
        rows = filereport(acc)
        if not rows:
            print(f"{acc}: no runs found", file=sys.stderr)
            continue
        for row in rows:
            ftps = [f for f in row["fastq_ftp"].split(";") if f]
            md5s = row["fastq_md5"].split(";")
            for ftp, md5 in zip(ftps, md5s + [""] * len(ftps)):
                download_one(ftp, md5, out)


def cmd_provenance(samples_tsv: str, out_path: str | None) -> None:
    """Audit how completely the archive documents each sample's origin.

    Task 2 asks us to assess provenance completeness, not merely to download.
    Ancient DNA papers carry archaeological context - site, date, burial,
    radiocarbon calibration - in the manuscript and its supplement, while the
    archive record frequently leaves the structured fields empty. Quantifying
    that gap is the finding, so we report which fields are actually populated.
    """
    wanted = ["collection_date", "country", "location", "sample_title",
              "library_construction_protocol", "instrument_model", "scientific_name"]

    runs: list[tuple[str, str]] = []
    with open(samples_tsv) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            runs.append((parts[0], parts[1].split(";")[0]))

    lines = ["sample\trun\t" + "\t".join(wanted) + "\tfields_populated\tcompleteness"]
    for sample, run in runs:
        rows = filereport(run)
        row = rows[0] if rows else {}
        vals = [(row.get(f) or "").strip() for f in wanted]
        n = sum(1 for v in vals if v)
        lines.append(
            f"{sample}\t{run}\t"
            + "\t".join(v if v else "(empty)" for v in vals)
            + f"\t{n}/{len(wanted)}\t{n / len(wanted):.2f}"
        )

    text = "\n".join(lines)
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        Path(out_path).write_text(text + "\n")
        print(f"wrote {out_path}", file=sys.stderr)
    print(text)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["metadata", "download", "provenance"])
    p.add_argument("accessions", nargs="+",
                   help="accessions, or the sample sheet path in provenance mode")
    p.add_argument("--out", default=None)
    a = p.parse_args()

    if a.mode == "metadata":
        cmd_metadata(a.accessions)
    elif a.mode == "download":
        cmd_download(a.accessions, a.out or "data/raw")
    else:
        cmd_provenance(a.accessions[0], a.out)


if __name__ == "__main__":
    main()
