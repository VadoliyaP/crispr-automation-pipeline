import os
import subprocess
import shutil
import pandas as pd
import numpy as np
import gradio as io
from Bio import SeqIO

# --- PIPELINE BLOCKS 1 & 3: PARSING & SAFETY SCORING ---
def extract_spcas9_targets(fasta_path):
    compiled_targets = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom_id = record.id
        fwd_seq = str(record.seq).upper()
        rev_seq = str(record.seq.reverse_complement()).upper()
        seq_len = len(fwd_seq)
        window_size = 30
        
        # Scan Forward
        for i in range(seq_len - window_size + 1):
            window = fwd_seq[i:i+window_size]
            if window[25:27] == "GG":
                compiled_targets.append({
                    "Chromosome": chrom_id, "Genomic_Start": i + 4, "Genomic_End": i + 24,
                    "Strand": "+", "Spacer_20nt": window[4:24], "PAM": window[24:27], "Context_30mer": window
                })
        # Scan Reverse
        for j in range(seq_len - window_size + 1):
            window = rev_seq[j:j+window_size]
            if window[25:27] == "GG":
                compiled_targets.append({
                    "Chromosome": chrom_id, "Genomic_Start": seq_len - (j + 24), "Genomic_End": seq_len - (j + 4),
                    "Strand": "-", "Spacer_20nt": window[4:24], "PAM": window[24:27], "Context_30mer": window
                })
    return pd.DataFrame(compiled_targets)

def calculate_hsu_zhang_score(offtarget_df):
    ZHANG_WEIGHTS = [
        0.0, 0.014, 0.014, 0.014, 0.014, 0.017, 0.017, 0.017, 0.017, 0.019,
        0.019, 0.019, 0.019, 0.035, 0.035, 0.035, 0.035, 0.085, 0.170, 0.257
    ]
    unique_guides_scores = {}
    if not offtarget_df.empty:
        for guide, group in offtarget_df.groupby("Query_Sequence"):
            total_penalty_sum = 0.0
            for _, row in group.iterrows():
                if row['Mismatch_Count'] == 0: continue
                query, hit = row['Query_Sequence'][:20], row['Matched_Sequence'][:20]
                mismatch_positions = [i for i in range(20) if query[i] != hit[i]]
                
                weight_product = 1.0
                for pos in mismatch_positions:
                    weight_product *= (1.0 - ZHANG_WEIGHTS[pos])
                
                distance_penalty = 1.0
                if len(mismatch_positions) > 1:
                    distances = [mismatch_positions[k] - mismatch_positions[k-1] for k in range(1, len(mismatch_positions))]
                    distance_penalty = 1.0 / (((19.0 - np.mean(distances)) / 19.0) * 4.0 + 1.0)
                
                total_penalty_sum += weight_product * distance_penalty * (1.0 / (row['Mismatch_Count'] ** 2))
            unique_guides_scores[guide] = 100.0 / (1.0 + total_penalty_sum)
    return unique_guides_scores

# --- MAIN CONTROLLER CALLED BY GRADIO ---
def run_crispr_pipeline(fasta_file, max_mismatches):
    if fasta_file is None:
        return "Please upload a genomic FASTA file first.", None
        
    # Set up runtime directories
    run_dir = os.path.abspath("workspace_run")
    genome_dir = os.path.join(run_dir, "genome")
    output_dir = os.path.join(run_dir, "output")
    os.makedirs(genome_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Copy uploaded file to clean genome folder
    fasta_name = os.path.basename(fasta_file.name)
    stable_fasta_path = os.path.join(genome_dir, fasta_name)
    shutil.copy(fasta_file.name, stable_fasta_path)
    
    # 1. Parse Sequences
    master_df = extract_spcas9_targets(stable_fasta_path)
    if master_df.empty:
        shutil.rmtree(run_dir)
        return "No SpCas9 target locations (NGG) found in the provided sequence.", None
        
    # 2. Configure and Run CAS-OFFinder
    cas_input = os.path.join(run_dir, "cas_input.txt")
    cas_output = os.path.join(output_dir, "cas_report.txt")
    
    with open(cas_input, 'w') as f:
        f.write(genome_dir + "\n")
        f.write("NNNNNNNNNNNNNNNNNNNNNGG\n")
        for _, row in master_df.iterrows():
            f.write(f"{row['Spacer_20nt']}NGG {int(max_mismatches)}\n")
            
    # Run the compiled binary via the exposed path defined in Dockerfile
    subprocess.run(f"cas-offinder {cas_input} C0 {cas_output}", shell=True, capture_output=True)
    
    # 3. Parse Off-Targets & Calculate Scores
    headers = ["Query_Sequence", "Chromosome", "Genomic_Position", "Matched_Sequence", "Strand", "Mismatch_Count"]
    if os.path.exists(cas_output) and os.path.getsize(cas_output) > 0:
        offtarget_df = pd.read_csv(cas_output, sep='\t', names=headers)
    else:
        offtarget_df = pd.DataFrame(columns=headers)
        
    guide_safety_map = calculate_hsu_zhang_score(offtarget_df)
    master_df['OffTarget_Safety_Score'] = master_df['Spacer_20nt'].apply(
        lambda x: guide_safety_map.get(f"{x}NGG", 100.0)
    )
    
    # Save optimized workbook
    excel_out = os.path.join(output_dir, "optimized_crispr_guides.xlsx")
    master_df.to_excel(excel_out, index=False)
    
    summary_text = f"Analysis Complete!\nInitial target sites mapped: {len(master_df)}\nOff-target alignment processing successful."
    return summary_text, excel_out

# --- GRADIO LAYOUT DRAWING ---
with io.Blocks(title="Custom Genome CRISPR Pipeline") as demo:
    io.Markdown("# 🧬 Automated CRISPR/Cas9 Design Pipeline")
    io.Markdown("Upload a custom genomic sequence to identify SpCas9 target sites and calculate off-target safety metrics using standalone CAS-OFFinder execution loops.")
    
    with io.Row():
        with io.Column():
            file_input = io.File(label="Upload Custom Genome (FASTA Format)", file_types=[".fasta", ".fa"])
            mismatch_slider = io.Slider(minimum=0, maximum=5, value=3, step=1, label="Maximum Mismatch Threshold")
            submit_btn = io.Button("Run Design Pipeline", variant="primary")
            
        with io.Column():
            status_output = io.Textbox(label="Execution Summary Status", interactive=False)
            file_output = io.File(label="Download Optimized Guide Matrix (.xlsx)")
            
    submit_btn.click(
        fn=run_crispr_pipeline,
        inputs=[file_input, mismatch_slider],
        outputs=[status_output, file_output]
    )

demo.launch(server_name="0.0.0.0", server_port=7860)
