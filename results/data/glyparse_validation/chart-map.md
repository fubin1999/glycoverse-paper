# Chart map

| Segment | Question | Family/type | Fields | Supported claim | Palette |
|---|---|---|---|---|---|
| Parse coverage | How often does each parser return a glycan? | Comparison / bar | format_label, parse_success_rate, rows | Coverage differs materially by format | Single blue root; direct labels |
| External coverage | How much external output is semantically comparable? | Comparison / grouped bar | format_label, normalized_rate, tool, direct_rows, fallback_rows | Tool coverage and format dialect support differ | Three approved roots plus labels |
| Throughput | What is steady-state parser throughput? | Ranking / bar | format_label, median_rows_per_second, sample_rows, repeats | GLYCAM normalization is slower than the other parsers | Single blue root; exact tooltips |

Sequence-level coverage is presented as a table because exact lookup across four mutually exclusive states is more important than visual shape.
