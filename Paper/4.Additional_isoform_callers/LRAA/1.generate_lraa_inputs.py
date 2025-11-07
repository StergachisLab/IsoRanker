import os
import json
import argparse

def generate_lraa_inputs(input_folder, annot_gtf, genome_fasta, output_folder):
    # Find all .tagged.bam files
    bam_files = [f for f in os.listdir(input_folder) if f.endswith(".tagged.bam")]

    main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

    os.makedirs(output_folder, exist_ok=True)

    for bam_file in bam_files:
        sample_id = bam_file.replace(".tagged.bam", "")  # you can customize how you define sample_id
        full_bam_path = os.path.join(input_folder, bam_file)

        inputs = {
            "LRAA_wf.sample_id": sample_id,
            "LRAA_wf.inputBAM": full_bam_path,
            "LRAA_wf.annot_gtf": annot_gtf,
            "LRAA_wf.referenceGenome": genome_fasta,
            "LRAA_wf.main_chromosomes": main_chromosomes
        }

        output_json = os.path.join(output_folder, f"{sample_id}.lraa_inputs.json")
        with open(output_json, "w") as f:
            json.dump(inputs, f, indent=4)

        print(f"Created input JSON: {output_json}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate LRAA input JSONs from tagged BAM files")
    parser.add_argument("--input_folder", required=True, help="Folder with .tagged.bam files")
    parser.add_argument("--annot_gtf", required=True, help="Path to annotation GTF file")
    parser.add_argument("--genome_fasta", required=True, help="Path to reference genome fasta file")
    parser.add_argument("--output_folder", required=True, help="Folder to save generated JSONs")
    args = parser.parse_args()

    generate_lraa_inputs(args.input_folder, args.annot_gtf, args.genome_fasta, args.output_folder)
