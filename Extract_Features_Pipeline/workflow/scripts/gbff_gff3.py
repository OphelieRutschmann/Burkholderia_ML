#!/usr/bin/env python3
"""
Convert GenBank GBFF format to GFF3 format.

This script:
1. Parses GenBank GBFF annotation file
2. Extracts sequence features (genes, CDSs, etc.)
3. Converts to GFF3 format with proper attribute escaping
4. Outputs GFF3 file

Input: GBFF file (GenBank format)
Output: GFF3 file (GFF version 3 format)
"""

import argparse
from Bio import SeqIO


def gff_escape(value):
    """Escape special characters for GFF3 attribute field."""
    return value.replace(";", "%3B").replace("=", "%3D").replace(",", "%2C")


def convert_gbff_to_gff3(input_file, output_file):
    """Convert GenBank GBFF file to GFF3 format."""
    with open(output_file, "w") as out_fh:
        for rec in SeqIO.parse(input_file, "genbank"):
            for f in rec.features:
                start = int(f.location.start) + 1
                end = int(f.location.end)
                strand = "+" if f.location.strand != -1 else "-"

                ftype = f.type

                attrs = []

                for key in ["locus_tag", "gene", "product", "protein_id"]:
                    if key in f.qualifiers:
                        attrs.append(
                            f"{key}={gff_escape(f.qualifiers[key][0])}"
                        )

                attr_str = ";".join(attrs) if attrs else "."

                print(
                    rec.id,
                    "GenBank",
                    ftype,
                    start,
                    end,
                    ".",
                    strand,
                    ".",
                    attr_str,
                    sep="\t",
                    file=out_fh
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert GenBank GBFF format to GFF3 format"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input GBFF file (GenBank format)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output GFF3 file"
    )

    args = parser.parse_args()
    convert_gbff_to_gff3(args.input, args.output)