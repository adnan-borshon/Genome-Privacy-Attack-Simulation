import gzip
import csv
import os

VCF_PATH = "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL_PATH = "integrated_call_samples_v3.20130502.ALL.panel"
OUT_CSV_PATH = "vcf_variants_preview.csv"
OUT_SAMPLES_CSV_PATH = "vcf_sample_metadata_summary.csv"

NUM_VARIANTS_TO_EXTRACT = 100
NUM_SAMPLES_TO_INCLUDE = 20  # First 20 samples to keep CSV easily viewable in Excel / Notepad

def parse_info_field(info_str):
    """Extract common fields like AF (Allele Frequency) and VT (Variant Type)."""
    fields = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            fields[k] = v
        else:
            fields[item] = True
    return fields

def simplify_genotype(gt_str):
    """
    VCF genotypes are phased like '0|0', '0|1', '1|0', '1|1'.
    Returns genotype string and allele dosage (0, 1, or 2).
    """
    gt = gt_str.split(":")[0]  # Just in case there are other sub-fields
    if gt in ("0|0", "0/0"):
        return f"{gt} (0)"
    elif gt in ("0|1", "1|0", "0/1", "1/0"):
        return f"{gt} (1)"
    elif gt in ("1|1", "1/1"):
        return f"{gt} (2)"
    return gt

def main():
    print(f"Opening VCF: {VCF_PATH} ...")
    
    # 1. Read panel metadata if available
    panel_data = {}
    if os.path.exists(PANEL_PATH):
        with open(PANEL_PATH, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                panel_data[row["sample"]] = row

    with gzip.open(VCF_PATH, "rt", encoding="utf-8") as vcf:
        header_cols = []
        for line in vcf:
            if line.startswith("#CHROM"):
                header_cols = line.strip().split("\t")
                break
        
        if not header_cols:
            print("Error: Could not find #CHROM header line.")
            return

        all_samples = header_cols[9:]
        selected_samples = all_samples[:NUM_SAMPLES_TO_INCLUDE]
        
        print(f"Total samples in dataset: {len(all_samples)}")
        print(f"Extracting first {NUM_VARIANTS_TO_EXTRACT} variants for first {len(selected_samples)} samples...")

        # Build output CSV header
        csv_headers = [
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", 
            "VARIANT_TYPE", "ALLELE_FREQ", "TOTAL_DEPTH"
        ] + selected_samples

        records = []
        for _ in range(NUM_VARIANTS_TO_EXTRACT):
            line = vcf.readline()
            if not line:
                break
            parts = line.strip().split("\t")
            chrom, pos, var_id, ref, alt, qual, filter_val, info, fmt = parts[:9]
            sample_gts = parts[9:9 + len(selected_samples)]
            
            info_dict = parse_info_field(info)
            vt = info_dict.get("VT", "N/A")
            af = info_dict.get("AF", "N/A")
            dp = info_dict.get("DP", "N/A")

            row = [
                chrom, pos, var_id, ref, alt, qual, filter_val,
                vt, af, dp
            ] + [simplify_genotype(gt) for gt in sample_gts]

            records.append(row)

    # 2. Write variants CSV
    with open(OUT_CSV_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(csv_headers)
        writer.writerows(records)

    print(f"Successfully saved {len(records)} variants to: {OUT_CSV_PATH}")

    # 3. Write metadata summary for the selected samples
    with open(OUT_SAMPLES_CSV_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Sample_ID", "Population", "Super_Population", "Gender"])
        for s in selected_samples:
            meta = panel_data.get(s, {})
            writer.writerow([
                s,
                meta.get("pop", "Unknown"),
                meta.get("super_pop", "Unknown"),
                meta.get("gender", "Unknown")
            ])
            
    print(f"Successfully saved sample demographic metadata to: {OUT_SAMPLES_CSV_PATH}")

if __name__ == "__main__":
    main()
