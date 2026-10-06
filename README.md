# The Genomic Time Machine: Authenticating Ancient DNA

Introduction to Bioinformatics, midterm project 23 (AITU).

**Deliverable:** an ancient-DNA authentication pipeline that reports, per sample, the
post-mortem damage profile, contamination estimate and endogenous content, and states
whether the library can be called genuinely ancient.

## Headline result

Seven real libraries, four expected behaviours. The damage estimator fires where it
should and stays silent where it should (`results/summary/authentication_report.md`):

| Sample | Library | Verdict | 5' C>T excess | Median fragment |
|---|---|---|---|---|
| TRP002.A | ancient *Y. pestis*, partial UDG | **ANCIENT** | 0.192 | 42 bp |
| LBG002.A | ancient *Y. pestis*, full UDG | UDG_TREATED_RESIDUAL_SIGNAL | 0.001 | 56 bp |
| MOD1953 | **modern** *Y. pestis* (1953 isolate) | MODERN_AS_EXPECTED | 0.000 | 305 bp |
| LVC005.A | ancient, protocol undeclared | WEAK_ANCIENT_SIGNAL | 0.005 (3' G>A 0.018) | 53 bp |
| LVC001.C, Motala1 | ancient, protocol undeclared | NO_DAMAGE_DETECTED | 0.003-0.008 | 43-54 bp |
| LP31b | historic *Y. pestis* | UNRELIABLE (control moved) | 0.015 | 49 bp |

See the report for what each of these does and does not mean.

## Reproduce

```bash
# 1. environment (locked, osx-arm64) - or use envs/environment.yml on other platforms
micromamba create -p ./.env --file envs/environment.lock.txt
export PATH="$PWD/.env/bin:$PATH"

# 2. ONE COMMAND: verify the pipeline from raw reads on the committed test set (~20 s)
bash scripts/run_test.sh          # ends with: TEST PASSED

# 3. full analysis (downloads ~1 GB from ENA, md5-verified; ~15 min on a laptop)
bash scripts/fetch_refs.sh
snakemake -s workflow/Snakefile --cores 6
```

`scripts/run_test.sh` works in a scratch directory, so it never touches `results/`.
A `Dockerfile` is provided but was **not built** (no Docker daemon was available during
development); the locked conda environment is the tested route.

## Tasks -> where they live

| # | Task | Course topic | Where |
|---|------|--------------|-------|
| 1 | Post-mortem damage biology | Intro | report section 2 |
| 2 | Retrieve data + audit provenance | Databases | `scripts/fetch_ena.py`, `results/summary/provenance.tsv` |
| 3 | Library prep: single- vs double-stranded | Platforms | report section 3 |
| 4 | Short-read QC adapted to aDNA | NGS | `rule trim` in `workflow/Snakefile`, `config/config.yaml` |
| 5 | Why long reads do not apply | TGS | `scripts/fragment_lengths.py`, `results/figures/fragment_lengths.png` |
| 6 | Permissive alignment + reference bias | Alignment | `rule align_*`, `scripts/reference_bias.py` |
| 7 | Damage-based authentication | Authentication | `scripts/damage_profile.py`, `scripts/summarize.py` |
| 8 | Contamination and endogenous content | Authentication | `scripts/authenticate.py` |

## Layout

```
config/      config.yaml (every non-default parameter justified), samples.tsv, samples_test.tsv
envs/        environment.yml, environment.lock.txt (explicit, osx-arm64)
scripts/     fetch_ena.py (retry/resume/md5), damage_profile.py, fragment_lengths.py,
             authenticate.py, reference_bias.py, summarize.py, plots.py, run_test.sh
workflow/    Snakefile
data/        raw/ ref/ (fetched, not committed); test/ (8.5 MB, committed)
results/     pipeline output (not committed)
report/      report.pdf and sources
notebooks/   analysis.ipynb (executed; reads results/, recomputes nothing)
```

## Known limitations (also in the report)

- **UDG-treated libraries cannot be authenticated by damage.** LBG002.A is reported as
  untestable, not as modern.
- The damage threshold (0.02) and `d_ref` (0.30) are design choices, not fitted values.
- The `no damage detected` verdict on the undeclared-protocol samples is ambiguous by
  construction: ENA does not state their library protocol.
- Peak memory was measured with `/usr/bin/time -l`; Snakemake's own memory benchmark
  reports `NA` on macOS.
- Contamination is estimated from damage partitioning and competitive depth, not from
  an X-chromosome or mtDNA haplotype test, because no sample here is a high-coverage
  human genome.

## Contributions

_TODO (must match the git history): who led which component._

## AI assistance disclosure

_TODO: see report appendix. Claude Code was used for scaffolding, pipeline code and drafting._
