# glyparse corpus validation and benchmark

## Technical summary

The benchmark evaluated **240,267** source sequences across seven formats. glyparse returned canonical glycans for **214,887** rows (**89.44%**). The remaining **25,380** failures have row-level stage, class, and reason evidence.

Among **104,619** accession-matched non-condensed rows, **104,527** were semantically comparable: **104,522** matched the independently normalized IUPAC-condensed reference and **5** differed; **92** failed in glyparse before comparison. Direct external evidence supports **2** glyparse outlier format-row(s), covering accession(s): G12345BK.

## Scope and definitions

IUPAC-condensed input is assessed only for parse coverage and performance, as requested. All other formats are checked against an accession-matched IUPAC-condensed corpus where available and against GlycanFormatConverter 2.10.3, glycowork 1.8.0, and glypy 1.0.17.

A semantic comparison requires both sides to normalize successfully through the current `parse_iupac_condensed()` and then have identical glyrepr canonical serialization. This removes branch and node-order differences while preserving residue identity, substituents, topology, linkages, reducing-end anomer and alditol state, floating components, and candidate-parent domains.

## Methodology

Each source sequence is parsed independently so one invalid or generic row cannot contaminate another row's result. Accession reference matching is exact. External converter availability, conversion success, current-model normalization, and semantic agreement are separate fields. Accessions converted through an accession-matched WURCS fallback are labeled as fallback and are never used to attribute a source-format error to glyparse.

Performance uses deterministic, evenly spaced samples of 1,000 successfully parsed concrete glycans, three vectorized repeats per parser. This isolates steady-state vector throughput from the scalar row-level diagnostic run.

## Findings

The normalized accession reference agrees with glyparse for **99.995%** of comparable non-condensed rows. The five reference differences represent two accessions and concern anomeric-carbon specificity on mannose residues. Direct external conversion supports glyparse being the outlier for 2 format-row(s); the other reference differences remain unresolved or support the reference/source distinction described in the evidence table.

GlycanFormatConverter direct coverage is format-dependent. It directly handles GlycoCT and WURCS in this run, but only 27 IUPAC-extended rows converted successfully; GLYCAM-IUPAC, GWB, and compact results from this tool use accession-matched WURCS fallback and are not direct source-parser evidence.

## Limitations and robustness

Canonical comparison is intentionally limited to chemistry representable by the current glyrepr model. An external output that cannot be normalized is reported as unavailable rather than as a semantic disagreement. Agreement between two direct converters without an accession reference is retained as an unresolved external-consensus difference, not attributed to glyparse, because shared information loss cannot be excluded. The external libraries differ in supported residue vocabularies and reducing-end conventions, so raw string equality is never used.

Detailed failure classifications reuse exact-row diagnostics only when accession and full source sequence match the prior diagnostic corpus; otherwise the current internal parser stage is replayed. Every final failure row retains the diagnostic source.

## Recommended next steps

Review the supported outlier accession(s) first, then the unresolved five-row reference ledger. For parser coverage, prioritize the largest glyrepr-restriction families separately from glyparse implementation limitations. Re-run the workflow after any parser or glyrepr semantic change; the source revisions and corpus hashes in the cache metadata make changes auditable.
