# glyparse corpus validation and benchmark

This workflow parses every sequence in
`data/benchmark/glycan_sequences`, benchmarks each parser on a deterministic
evenly spaced sample, and preserves one result row per source row.

## Validation contract

1. IUPAC-condensed source strings are tested for parse coverage with glyparse,
   GlycanFormatConverter, glycowork, and glypy, but no correctness claim is made
   for that source format.
2. Every other successful glyparse result is serialized with `as.character()`.
   When the same GlyTouCan accession exists in the IUPAC-condensed corpus, the
   reference is independently normalized with `parse_iupac_condensed()` and the
   two canonical glyrepr serializations are compared.
3. GlycanFormatConverter, glycowork, and glypy outputs are retained raw, then
   normalized through current `parse_iupac_condensed()`. Converter availability,
   conversion success, normalization success, and semantic agreement are
   separate fields.
4. Canonical glyrepr serialization is the primary semantic key. Accession
   references and GlycanFormatConverter are compared strictly. Because
   glycowork and glypy do not encode the reducing-end anomer, their comparison
   keys remove the terminal reducing-end token from both the glyparse and tool
   canonical strings (for example, `Gal(b1-3)GalNAc(a1-` and
   `Gal(b1-3)GalNAc(?1-` are both compared as `Gal(b1-3)GalNAc`). Full
   normalized strings and explicit comparison keys are retained for audit.
   Internal anomers/linkages, residue identities, substituents, topology,
   reducing-end alditol state, floating components, floating parent domains,
   and floating substituents remain part of the semantic comparison.
5. A difference is assigned to glyparse only when the parsed accession-matched
   reference or a supported external consensus demonstrates that glyparse is
   the outlier. External information loss and unsupported syntax remain external
   evidence limitations.

## Run

From the repository root:

```sh
Rscript src/benchmark/glyparse.R
```

The stages can be resumed independently:

```sh
Rscript src/benchmark/glyparse.R --stage glyparse
Rscript src/benchmark/glyparse.R --stage external
Rscript src/benchmark/glyparse.R --stage adjudicate
Rscript src/benchmark/glyparse.R --stage report
```

Environment overrides:

- `GLYPARSE_SOURCE_DIR`
- `GLYREPR_SOURCE_DIR`
- `GLYCAN_FORMAT_CONVERTER_JAR`
- `GLYPARSE_GFC_REUSE_CACHE`
- `GLYPARSE_BENCHMARK_ROOT`

Generated cache files and final audit artifacts are written to
`results/data/glyparse_validation`.

Final artifacts include:

- `format_summary.csv` and `parser_benchmark.csv`
- `external_tool_summary.csv`
- `sequence_coverage.csv.gz`, with one row per tested source sequence and
  direct parseability, evidence-scope, and failure-reason fields for glyparse,
  GlycanFormatConverter, glypy, and glycowork
- `sequence_coverage_summary.csv`, with format-by-tool coverage counts
- `parse_failures.csv.gz`, with one detailed reason per failed source row
- `semantic_disagreements.csv.gz`, with one component-level reason per
  normalized glyparse/comparator pair
- `validation_rows.csv.gz`, the complete joined evidence ledger
- `report.md`, `artifact.json`, and the self-contained `report.html`

Sequence coverage uses four statuses: `parsed` and `failed` are direct tests of
the source sequence; `fallback_only` means the available result came from an
accession-matched WURCS sequence rather than the source format; and `not_tested`
means the tool was not run for that format. The corresponding `can_parse`
field is `TRUE` or `FALSE` only for direct tests and blank otherwise. For glypy
IUPAC-condensed coverage, the dangling reducing-end token is removed before
calling its simple-IUPAC parser because glypy does not encode reducing-end
anomers; this normalization is recorded in the row's evidence-scope field.

The adjudication stage can reuse exact-row diagnostic evidence from a previous
audit when `GLYPARSE_PRIOR_DIAGNOSTIC_DIR` points to it. If that directory is
unavailable or the accession plus full source sequence does not match, the
current internal parser stage is replayed instead.

Python converter caches are reused only when every cached row records the same
installed tool version as the active conda environment. A glycowork or glypy
upgrade therefore forces fresh raw conversion before normalization.

To package the portable report after running the `report` stage, use the
Data Analytics report builder:

```sh
node /path/to/data-analytics/skills/build-report/scripts/deliver_portable_artifact.mjs \
  --input results/data/glyparse_validation/artifact.json \
  --output results/data/glyparse_validation/report.html
```
