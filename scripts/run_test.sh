#!/usr/bin/env bash
# One command to verify the pipeline from raw reads on the small committed test set.
#   bash scripts/run_test.sh
# Runs in a scratch copy so it never touches your real results/. Expect ~2 minutes.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
bash "$REPO/scripts/fetch_refs.sh" >/dev/null
W="$(mktemp -d)"; trap 'echo "scratch dir: $W"' EXIT
mkdir -p "$W/data/raw" "$W/data/ref" "$W/config"
cp -R "$REPO/scripts" "$REPO/workflow" "$W/"
cp "$REPO/config/config.yaml" "$W/config/"
cp "$REPO/config/samples_test.tsv" "$W/config/samples.tsv"
cp "$REPO/data/ref/ypestis_co92.fa" "$W/data/ref/"
for f in "$REPO"/data/test/*.fastq.gz; do
  b=$(basename "$f" | sed 's/_test_\([12]\)/_test_\1/'); cp "$f" "$W/data/raw/$b"
done
cd "$W"
snakemake -s workflow/Snakefile --cores "${CORES:-4}" \
  results/summary/authentication_table.tsv results/figures/damage_profiles.png \
  --config offline=true
echo; cat results/summary/authentication_report.md | sed -n '1,8p'
python - <<'PY'
import json
rows = {r["sample"]: r["verdict"] for r in json.load(open("results/summary/summary.json"))}
expect = {"TRP002.A": "ANCIENT", "LBG002.A": "UDG_TREATED_RESIDUAL_SIGNAL", "MOD1953": "MODERN_AS_EXPECTED"}
bad = {k: (rows.get(k), v) for k, v in expect.items() if rows.get(k) != v}
print("\nTEST", "PASSED" if not bad else f"FAILED {bad}")
raise SystemExit(1 if bad else 0)
PY
