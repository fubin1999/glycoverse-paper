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
```

Environment overrides:

- `GLYPARSE_SOURCE_DIR`
- `GLYREPR_SOURCE_DIR`
- `GLYCAN_FORMAT_CONVERTER_JAR`
- `GLYPARSE_GFC_REUSE_CACHE`
- `GLYPARSE_BENCHMARK_ROOT`

Generated cache files and final audit artifacts are written to
`data/benchmark/glyparse_validation`.
