import os
import sys
import pandas as pd
import warnings

warnings.filterwarnings("ignore")


def merge_demographics():
    # Base directory (always points to the project root where this script lives)
    project_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(project_dir, "Genome dataset")

    panel_file = os.path.join(data_dir, "integrated_call_samples_v3.20130502.ALL.panel")
    ped_file = os.path.join(data_dir, "integrated_call_samples_v3.20250704.ALL.ped")
    output_file = os.path.join(project_dir, "merged_sample_demographics.csv")

    print("=" * 65)
    print("1000 Genomes Project: Demographic & Pedigree Merger")
    print("=" * 65)
    print(f"Project Directory : {project_dir}")
    print(f"Dataset Directory : {data_dir}")

    # Check existence of input files
    for name, path in [("Panel File", panel_file), ("Pedigree File", ped_file)]:
        if not os.path.exists(path):
            print(f"[!] Error: {name} not found at: {path}")
            sys.exit(1)
        print(f"[OK] Found {name}: {os.path.basename(path)}")

    # 1. Load Demographic Panel
    print("\n[1/3] Loading panel data (population labels)...")
    panel = pd.read_csv(
        panel_file,
        sep="\t",
        usecols=[0, 1, 2, 3],
        names=["Sample_ID", "Population", "Super_Population", "Gender"],
        header=0,
    )
    print(f"      Panel rows: {len(panel)}, columns: {list(panel.columns)}")

    # 2. Load Pedigree Data
    print("\n[2/3] Loading pedigree data (family structure & sex)...")
    ped = pd.read_csv(
        ped_file,
        sep="\t",
        usecols=[0, 1, 2, 3, 4],
        header=0,
    )
    ped.columns = ["Family_ID", "Sample_ID", "Father_ID", "Mother_ID", "Sex_Code"]
    # Sex_Code: 1 = male, 2 = female
    ped["Sex_Label"] = ped["Sex_Code"].map({1: "male", 2: "female"}).fillna("unknown")
    print(f"      Pedigree rows: {len(ped)}, columns: {list(ped.columns)}")

    # 3. Merge Panel and Pedigree on Sample_ID
    print("\n[3/3] Merging panel and pedigree on 'Sample_ID'...")
    merged_df = panel.merge(
        ped[["Sample_ID", "Family_ID", "Father_ID", "Mother_ID", "Sex_Label"]],
        on="Sample_ID",
        how="left",
    )
    print(f"      Merged rows: {len(merged_df)}")

    # Save to CSV in main project directory
    merged_df.to_csv(output_file, index=False)
    file_size_kb = os.path.getsize(output_file) / 1024
    print("\n" + "-" * 65)
    print(f"[OK] Saved merged file: {output_file}")
    print(f"    Size: {file_size_kb:.2f} KB | Total Records: {len(merged_df)}")
    print("-" * 65)

    # Summary Statistics
    print("\nSuper-Population Counts:")
    print(merged_df["Super_Population"].value_counts().to_string())

    print("\nGender Counts:")
    print(merged_df["Gender"].value_counts().to_string())

    print("\nPreview of first 5 rows:")
    print(merged_df.head().to_string(index=False))


if __name__ == "__main__":
    merge_demographics()
