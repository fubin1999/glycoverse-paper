# glyparse corpus validation and benchmark

This workflow parses every sequence in
`data/benchmark/glycan_sequences`, benchmarks each parser on a deterministic
evenly spaced sample, and preserves one result row per source row.

## Validation contract

1. `parse_iupac_condensed()` is measured for parse coverage and speed only.
   No correctness claim is made for that source format.
2. Every other successful glyparse result is serialized with `as.character()`.
   When the same GlyTouCan accession exists in the IUPAC-condensed corpus, the
   reference is independently normalized with `parse_iupac_condensed()` and the
   two canonical glyrepr serializations are compared.
3. GlycanFormatConverter, glycowork, and glypy outputs are retained raw, then
   normalized through current `parse_iupac_condensed()`. Converter availability,
   conversion success, normalization success, and semantic agreement are
   separate fields.
4. Canonical glyrepr serialization is the primary semantic key. It removes
   non-semantic branch and node-order differences while preserving residue
   identities, substituents, linkages, reducing-end anomer/alditol state,
   floating components, floating parent domains, and floating substituents.
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
`data/benchmark/glyparse_validation`.

Final artifacts include:

- `format_summary.csv` and `parser_benchmark.csv`
- `external_tool_summary.csv`
- `parse_failures.csv.gz`, with one detailed reason per failed source row
- `semantic_disagreements.csv.gz`, with one component-level reason per
  normalized glyparse/comparator pair
- `validation_rows.csv.gz`, the complete joined evidence ledger
- `report.md`, `artifact.json`, and the self-contained `report.html`

The adjudication stage can reuse exact-row diagnostic evidence from a previous
audit when `GLYPARSE_PRIOR_DIAGNOSTIC_DIR` points to it. If that directory is
unavailable or the accession plus full source sequence does not match, the
current internal parser stage is replayed instead.

To package the portable report after running the `report` stage, use the
Data Analytics report builder:

```sh
node /path/to/data-analytics/skills/build-report/scripts/deliver_portable_artifact.mjs \
  --input data/benchmark/glyparse_validation/artifact.json \
  --output data/benchmark/glyparse_validation/report.html
```
