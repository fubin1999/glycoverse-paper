#!/usr/bin/env python3

import argparse
import csv
import gzip
import importlib.metadata
import os
import re
import sys
import time
import traceback
import warnings


FORMATS = {
    "glycam_iupac": (
        "glycan_sequences_glycam_iupac.csv",
        "sequence_glycam_iupac",
    ),
    "glycoct": ("glycan_sequences_glycoct.csv", "sequence_glycoct"),
    "gwb": ("glycan_sequences_gwb.csv", "sequence_gwb"),
    "iupac_compact": (
        "glycan_sequences_iupac_compact.csv",
        "sequence_iupac_compact",
    ),
    "iupac_extended": (
        "glycan_sequences_iupac_extended.csv",
        "sequence_iupac_extended",
    ),
    "wurcs": ("glycan_sequences_wurcs.csv", "sequence_wurcs"),
}


def normalize_error(error):
    text = " ".join(str(error).replace("\r", " ").replace("\n", " ").split())
    return text or error.__class__.__name__


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as stream:
        yield from csv.DictReader(stream)


def read_wurcs_by_accession(corpus_dir):
    path = os.path.join(corpus_dir, FORMATS["wurcs"][0])
    return {
        row["glytoucan_ac"]: row[FORMATS["wurcs"][1]] for row in read_rows(path)
    }


def glycowork_converter(format_name):
    from glycowork.motif.processing import (
        canonicalize_iupac,
        glycam_to_iupac,
        glycoct_to_iupac,
        glycoworkbench_to_iupac,
        iupac_extended_to_condensed,
        wurcs_to_iupac,
    )

    def convert(sequence, accession, wurcs_by_accession):
        del accession, wurcs_by_accession
        if format_name == "glycam_iupac":
            raw = glycam_to_iupac(sequence)
            scope = "direct_source_converter"
        elif format_name == "glycoct":
            raw = glycoct_to_iupac(sequence.replace(" ", "\n"))
            scope = "direct_source_converter"
        elif format_name == "gwb":
            raw = glycoworkbench_to_iupac(sequence)
            scope = "direct_source_converter"
        elif format_name == "iupac_extended":
            raw = iupac_extended_to_condensed(sequence)
            scope = "direct_source_converter"
        elif format_name == "wurcs":
            raw = wurcs_to_iupac(sequence)
            scope = "direct_source_converter"
        elif format_name == "iupac_compact":
            raw = sequence
            scope = "direct_general_canonicalizer"
        else:
            raise ValueError(f"Unsupported format: {format_name}")
        return canonicalize_iupac(raw), scope

    return convert


def normalize_open_extended(sequence):
    sequence = re.sub(r"-\((?:1|2)(?:->|→)$", "", sequence)
    return sequence.rstrip("-")


def glypy_converter(format_name):
    from glypy.io import glycoct, gws, iupac, wurcs

    def serialize(structure):
        return iupac.dumps(structure, dialect="simple")

    def from_wurcs(sequence):
        return serialize(wurcs.loads(sequence))

    def convert(sequence, accession, wurcs_by_accession):
        if format_name == "glycoct":
            structure = glycoct.loads(sequence.replace(" ", "\n"))
            return serialize(structure), "direct_source_converter"
        if format_name == "gwb":
            structure, _metadata = gws.loads(sequence)
            return serialize(structure), "direct_source_converter"
        if format_name == "iupac_extended":
            normalized = normalize_open_extended(sequence)
            structure = iupac.loads(normalized, dialect="extended")
            return (
                serialize(structure),
                "direct_with_open_reducing_end_normalization",
            )
        if format_name == "wurcs":
            return from_wurcs(sequence), "direct_source_converter"
        if format_name in {"glycam_iupac", "iupac_compact"}:
            matched = wurcs_by_accession.get(accession)
            if not matched:
                raise LookupError("No accession-matched WURCS sequence is available")
            return from_wurcs(matched), "accession_matched_wurcs_fallback"
        raise ValueError(f"Unsupported format: {format_name}")

    return convert


def run_format(tool, format_name, corpus_dir, output_dir, wurcs_by_accession):
    file_name, sequence_column = FORMATS[format_name]
    input_path = os.path.join(corpus_dir, file_name)
    output_path = os.path.join(output_dir, f"{tool}-{format_name}.csv.gz")
    converter = (
        glycowork_converter(format_name)
        if tool == "glycowork"
        else glypy_converter(format_name)
    )
    version = importlib.metadata.version(tool)
    fields = [
        "tool",
        "tool_version",
        "format",
        "row_index",
        "glytoucan_ac",
        "evidence_scope",
        "conversion_status",
        "raw_iupac_condensed",
        "conversion_warning",
        "conversion_error",
    ]
    started = time.perf_counter()
    converted = 0
    failed = 0
    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row_index, row in enumerate(read_rows(input_path), start=1):
            captured_warnings = []
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    output, scope = converter(
                        row[sequence_column],
                        row["glytoucan_ac"],
                        wurcs_by_accession,
                    )
                    captured_warnings = [str(item.message) for item in caught]
                if not isinstance(output, str) or not output.strip():
                    raise ValueError("Converter returned an empty result")
                status = "converted"
                error = ""
                converted += 1
            except Exception as condition:
                output = ""
                scope = (
                    "accession_matched_wurcs_fallback"
                    if tool == "glypy"
                    and format_name in {"glycam_iupac", "iupac_compact"}
                    else "direct_source_converter"
                )
                status = "failed"
                error = normalize_error(condition)
                failed += 1
            writer.writerow(
                {
                    "tool": tool,
                    "tool_version": version,
                    "format": format_name,
                    "row_index": row_index,
                    "glytoucan_ac": row["glytoucan_ac"],
                    "evidence_scope": scope,
                    "conversion_status": status,
                    "raw_iupac_condensed": output,
                    "conversion_warning": " | ".join(captured_warnings),
                    "conversion_error": error,
                }
            )
            if row_index % 1000 == 0:
                print(
                    f"{tool} {format_name}: {row_index} rows complete",
                    file=sys.stderr,
                    flush=True,
                )
    elapsed = time.perf_counter() - started
    return {
        "tool": tool,
        "tool_version": version,
        "format": format_name,
        "rows": converted + failed,
        "converted": converted,
        "failed": failed,
        "elapsed_seconds": f"{elapsed:.6f}",
        "rows_per_second": f"{(converted + failed) / elapsed:.6f}",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", choices=("glycowork", "glypy"), required=True)
    parser.add_argument("--corpus-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    arguments = parser.parse_args()
    os.makedirs(arguments.output_dir, exist_ok=True)
    csv.field_size_limit(sys.maxsize)
    wurcs_by_accession = read_wurcs_by_accession(arguments.corpus_dir)
    summaries = []
    for format_name in FORMATS:
        summaries.append(
            run_format(
                arguments.tool,
                format_name,
                arguments.corpus_dir,
                arguments.output_dir,
                wurcs_by_accession,
            )
        )
    summary_path = os.path.join(
        arguments.output_dir,
        f"{arguments.tool}-run-summary.csv",
    )
    with open(summary_path, "w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
