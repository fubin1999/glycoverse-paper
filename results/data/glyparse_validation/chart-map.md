# Chart map

| Segment | Question | Family/type | Fields | Supported claim | Palette |
|---|---|---|---|---|---|
| Parse coverage | How often does each parser return a glycan? | Comparison / bar | format_label, parse_success_rate, rows | Coverage differs materially by format | Single blue root; direct labels |
| Sequence coverage | What share of directly tested sequences can each tool parse in each format? | Comparison / grouped bar | format_label, direct_parse_rate, tool, directly_tested, parsed | Direct coverage varies by tool and format; missing bars are untested | Four approved roots plus legend |
| External coverage | How much external output is semantically comparable? | Comparison / grouped bar | format_label, normalized_rate, tool, direct_rows, fallback_rows | Tool coverage and format dialect support differ | Three approved roots plus labels |
| Throughput | What is steady-state parser throughput? | Ranking / bar | format_label, median_rows_per_second, sample_rows, repeats | Throughput varies by parser; exact medians and repeat ranges are retained | Single blue root; exact tooltips |

The sequence-level coverage table remains alongside the chart for exact lookup across parsed, failed, fallback-only, and not-tested states.
