#!/usr/bin/env python3
"""
Convert snippy-core output to VCF format for annotation with SNPEff.

Input: core.tab (SNP table from snippy-core)
Output: core_snps.vcf (VCF format)
"""

import sys
from datetime import datetime

def core_tab_to_vcf(core_tab, ref_fasta, output_vcf):
    """Convert core.tab to VCF format."""

    # Get reference info from FASTA header
    ref_id = "unknown"
    with open(ref_fasta, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith('>'):
            ref_id = first_line[1:].split()[0]

    # Read core.tab and write VCF
    with open(output_vcf, 'w') as vcf:
        # Write VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        vcf.write("##source=snippy-core\n")
        vcf.write(f"##reference={ref_fasta}\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Parse core.tab and write SNPs
        # core.tab format: CHROM, POS, TYPE, REF, ALT, EVIDENCE, ...
        with open(core_tab, 'r') as f:
            next(f)  # Skip header
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue

                chrom = fields[0]
                pos = fields[1]
                snp_type = fields[2]
                ref = fields[3]
                alt = fields[4]

                # Skip non-SNP entries (like deletions)
                if snp_type != 'snp':
                    continue

                # Format: CHROM POS ID REF ALT QUAL FILTER INFO
                vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n"
                vcf.write(vcf_line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <core.tab> <ref.fasta> <output.vcf>")
        sys.exit(1)

    core_tab_to_vcf(sys.argv[1], sys.argv[2], sys.argv[3])
