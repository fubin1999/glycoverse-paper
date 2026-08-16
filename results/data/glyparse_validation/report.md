# glyparse corpus validation and benchmark

## Technical summary

The benchmark evaluated **240,267** source sequences across seven formats. glyparse returned canonical glycans for **220,320** rows (**91.70%**). The remaining **19,947** failures have row-level stage, class, and reason evidence.

IUPAC-compact, IUPAC-condensed, and IUPAC-extended each achieved **100% glyparse parse coverage**. The remaining failure ledger contains **0** row(s) whose detailed reason cites the former mixed generic/concrete residue restriction.

Among **104,619** accession-matched non-condensed rows, **104,617** were semantically comparable: **104,617** matched the independently normalized IUPAC-condensed reference and **0** differed; **2** failed in glyparse before comparison. Direct external evidence supports **0** glyparse outlier format-row(s), covering accession(s): none.

A separate sequence-level coverage ledger records whether each of the four tested tools parsed each exact source sequence. Direct failures, accession-fallback-only results, and untested formats remain distinct.

## Scope and definitions

IUPAC-condensed input is assessed with all four tools for parse coverage and with glyparse for performance, but is excluded from correctness adjudication as requested. All other formats are checked against an accession-matched IUPAC-condensed corpus where available and against GlycanFormatConverter 2.10.3, glycowork 1.10.0 (git dev@d132d2fda005), and glypy 1.0.17.

A semantic comparison requires both sides to normalize successfully through the current `parse_iupac_condensed()`. Chemically fixed donor positions are then filled with `fill_anomer_pos()`, and equality is tested with an isomorphism-invariant labeled-graph fingerprint, so implicit donor positions and alternate serializations of symmetric branches compare fairly. Accession references and GlycanFormatConverter retain the reducing-end anomer. Because glycowork and glypy do not encode it, only that graph-level attribute is ignored for those comparisons. This makes `Gal(b1-3)GalNAc(a1-` and `Gal(b1-3)GalNAc(?1-` equivalent while retaining internal anomers/linkages, residue identity, substituents, topology, reducing-end alditol state, floating components, and candidate-parent domains.

## Methodology

Each source sequence is parsed independently so one invalid or generic row cannot contaminate another row's result. Accession reference matching is exact. External converter availability, conversion success, current-model normalization, comparison contract, display key, and semantic agreement are separate fields. Labeled-graph fingerprints encode the rooted residue/linkage graph, substituents, alditol state, floating components, and candidate-parent relations. Direct cross-tool consensus uses the same graph comparison while ignoring only the terminal reducing-end anomer, so missing root-anomer information cannot create a false outlier. Accessions converted through an accession-matched WURCS fallback are labeled as fallback and are never used to attribute a source-format error to glyparse.

Performance uses deterministic, evenly spaced samples of 1,000 successfully parsed concrete glycans, three vectorized repeats per parser. This isolates steady-state vector throughput from the scalar row-level diagnostic run.

## Direct sequence coverage keeps fallback evidence separate

The row-level coverage report contains one record for every source sequence. `parsed` and `failed` are direct tests of that exact source string. `fallback_only` means the available conversion used an accession-matched WURCS sequence instead, and `not_tested` means the tool was not run for that source format. Boolean `can_parse` values are therefore populated only for direct tests. For glypy IUPAC-condensed input, the terminal open reducing-end token is removed before its simple-IUPAC parser is called, matching the agreed convention that glypy does not encode reducing-end anomers; the row retains this evidence scope.

- **glyparse:** parsed 220,320 of 240,267 directly tested rows (91.70%); 0 fallback-only and 0 not tested.
- **GlycanFormatConverter:** parsed 101,938 of 158,068 directly tested rows (64.49%); 82,199 fallback-only and 0 not tested.
- **glypy:** parsed 165,658 of 208,848 directly tested rows (79.32%); 31,419 fallback-only and 0 not tested.
- **glycowork:** parsed 230,344 of 240,267 directly tested rows (95.87%); 0 fallback-only and 0 not tested.

## Findings

The normalized accession reference agrees with glyparse for **100%** of the 104,617 comparable non-condensed rows. The previously reported G12345BK differences were false positives: GlycoCT, WURCS, and compact assigned an unknown linkage to different but automorphic terminal Man nodes, while GWB explicitly wrote Man's chemically fixed donor position where the reference omitted it. Chemistry-aware labeled-graph isomorphism now recognizes all four serializations as the same glycan.

GlycanFormatConverter direct coverage is format-dependent. It directly attempts GlycoCT, WURCS, IUPAC-extended, and IUPAC-condensed in this run; 19,483 of 19,621 IUPAC-condensed rows converted, while only 27 IUPAC-extended rows converted. GLYCAM-IUPAC, GWB, and compact results from this tool use accession-matched WURCS fallback and are not direct source-parser evidence.

## Limitations and robustness

Graph comparison is intentionally limited to chemistry representable by the current glyrepr model. An external output that cannot be normalized is reported as unavailable rather than as a semantic disagreement. Agreement between two direct converters without an accession reference is retained as an unresolved external-consensus difference, not attributed to glyparse, because shared information loss cannot be excluded. The reducing-end-neutral contract for glycowork and glypy ignores only the graph-level terminal root anomer; it does not forgive internal anomer/linkage, alditol, topology, residue, substituent, or floating-component differences. Raw string equality is only a fast path; every unequal serialization that could share the same residue/linkage inventory is adjudicated by labeled-graph isomorphism.

Detailed failure classifications reuse exact-row diagnostics only when accession and full source sequence match the prior diagnostic corpus; otherwise the current internal parser stage is replayed. Every final failure row retains the diagnostic source.

## Recommended next steps

No accession-reference semantic outlier remains for manual review. For parser coverage, prioritize the largest glyrepr-restriction families separately from glyparse implementation limitations. Re-run the workflow after any parser or glyrepr semantic change; the source revisions and corpus hashes in the cache metadata make changes auditable.
