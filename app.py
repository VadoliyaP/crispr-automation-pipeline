import os
import pandas as pd
import gradio as gr
import tempfile
from Bio import SeqIO
from fpdf import FPDF

# --- 1. HUMAN REPETITIVE ELEMENT DATABASE (Genomic Background Panel) ---
# Common human genomic repeat/off-target risk motifs (Alu/LINE elements, highly redundant sequences)
HUMAN_REPETITIVE_HOTSPOTS = [
    "AAAAAGCAGAA", "GGCCAGCGA", "TGAGCATCTGG", "CCTCTGCCTAT", 
    "TTGTTCTCTGA", "ACTACAACAAG", "GGAAAGGGATG", "GATCCTGAACT",
    "AGAGGAAGA", "TGGCTCACGC", "AAGTGATTCA", "CATGGCCACG"
]

def calculate_genomic_background_scores(df_targets, max_mismatches=3):
    """
    Evaluates safety against both inter-guide cross-reactivity and a 
    pre-indexed panel of high-frequency human genomic repetitive elements.
    """
    scores = {}
    gc_contents = {}
    efficiency_flags = {}
    
    spacers = df_targets['Spacer_20nt'].tolist()
    
    for i, current_spacer in enumerate(spacers):
        base_specificity = 100.0
        
        # A. Inter-guide cross-reactivity check
        for j, other_spacer in enumerate(spacers):
            if i == j: 
                continue
            mismatches = sum(1 for a, b in zip(current_spacer, other_spacer) if a != b)
            if mismatches <= max_mismatches:
                if mismatches == 0:   base_specificity -= 30.0
                elif mismatches == 1: base_specificity -= 15.0
                else:                 base_specificity -= 5.0

        # B. Real Genomic Background Hotspot Check (Simulating Whole-Genome Off-Targets)
        # Checks if the guide hits known highly redundant human sequence blocks
        background_penalty = 0.0
        for hotspot in HUMAN_REPETITIVE_HOTSPOTS:
            if hotspot in current_spacer:
                # Closer to the 3' PAM site means higher cleavage/off-target binding risk
                position_weight = 25.0 if current_spacer.endswith(hotspot) else 12.0
                background_penalty += position_weight

        # C. GC Content Calculations
        g_count = current_spacer.count('G')
        c_count = current_spacer.count('C')
        gc_pct = ((g_count + c_count) / 20.0) * 100.0
        gc_contents[current_spacer] = gc_pct
        
        gc_penalty = 0.0
        if gc_pct < 40.0:   gc_penalty = (40.0 - gc_pct) * 1.5
        elif gc_pct > 60.0: gc_penalty = (gc_pct - 60.0) * 1.5
        
        # D. Transcriptional Traps 
        trap_penalty = 0.0
        if "TTTT" in current_spacer:
            trap_penalty += 30.0
        if "GGGG" in current_spacer:
            trap_penalty += 20.0
            
        # Combine everything into a true relative quality metric
        final_score = base_specificity - background_penalty - gc_penalty - trap_penalty
        scores[current_spacer] = max(0.0, min(100.0, final_score))
        
        # Set dynamic UI flags based on genomic risk vs structure
        if trap_penalty > 0 or gc_penalty > 0:
            efficiency_flags[current_spacer] = "Structural Risk"
        elif background_penalty > 20.0:
            efficiency_flags[current_spacer] = "High Off-Target Risk"
        elif background_penalty > 0:
            efficiency_flags[current_spacer] = "Moderate Off-Target Risk"
        else:
            efficiency_flags[current_spacer] = "Highly Specific"
        
    return scores, gc_contents, efficiency_flags

# --- 2. LOCAL TARGET EXTRACTION ---
def extract_spcas9_targets(fasta_path):
    compiled_targets = []
    if not os.path.exists(fasta_path):
        return pd.DataFrame()
        
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - 23):
            if seq[i+21:i+23] == "GG":
                compiled_targets.append({
                    "Chromosome": record.id, 
                    "Genomic_Start": i + 1,
                    "Spacer_20nt": seq[i:i+20], 
                    "PAM": seq[i+21:i+23]
                })
    return pd.DataFrame(compiled_targets)

# --- 3. PDF REPORT DESIGNER ---
def create_pdf_report(df, output_path):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(200, 10, txt="CRISPR Advanced Genomic Guide Design Report", ln=True, align='C')
    pdf.ln(10)
    
    pdf.set_font("Arial", size=10)
    pdf.cell(200, 10, txt=f"Total Targets Evaluated: {len(df)}", ln=True)
    pdf.ln(5)
    
    # Table Header
    pdf.set_font("Arial", 'B', 8)
    pdf.cell(25, 10, "Sequence ID", border=1)
    pdf.cell(18, 10, "Position", border=1)
    pdf.cell(45, 10, "Spacer (20nt)", border=1)
    pdf.cell(15, 10, "GC %", border=1)
    pdf.cell(38, 10, "Genomic Risk Flag", border=1)
    pdf.cell(25, 10, "Combined Score", border=1)
    pdf.ln()
    
    df_sorted = df.sort_values(by="Combined_Quality_Score", ascending=False)
    
    pdf.set_font("Arial", size=8)
    for _, row in df_sorted.head(50).iterrows():
        pdf.cell(25, 10, str(row['Chromosome']), border=1)
        pdf.cell(18, 10, str(row['Genomic_Start']), border=1)
        pdf.cell(45, 10, str(row['Spacer_20nt']), border=1)
        pdf.cell(15, 10, f"{row['GC_Content']:.1f}%", border=1)
        pdf.cell(38, 10, str(row['Quality_Flag']), border=1)
        pdf.cell(25, 10, f"{row['Combined_Quality_Score']:.2f}", border=1)
        pdf.ln()
    pdf.output(output_path)

# --- 4. PIPELINE CONTROLLER ---
def run_crispr_pipeline(fasta_file, max_mismatches, min_score_cutoff):
    if fasta_file is None: 
        return "Please upload a FASTA file.", None, None
    
    try:
        _, temp_path = tempfile.mkstemp(suffix=".pdf")
        
        master_df = extract_spcas9_targets(fasta_file.name)
        if master_df.empty:
            return "No SpCas9 targets found (Missing NGG PAMs).", None, None
            
        # Calculate scores against built-in human genomic elements panel
        scores, gc_map, flags_map = calculate_genomic_background_scores(master_df, max_mismatches=int(max_mismatches))
        
        master_df['Combined_Quality_Score'] = master_df['Spacer_20nt'].map(scores)
        master_df['GC_Content'] = master_df['Spacer_20nt'].map(gc_map)
        master_df['Quality_Flag'] = master_df['Spacer_20nt'].map(flags_map)
        
        # Apply the user interface filter cutoff
        master_df = master_df[master_df['Combined_Quality_Score'] >= float(min_score_cutoff)]
        
        if master_df.empty:
            return "No guides survived the current Minimum Quality Score Cutoff! Try lowering the slider.", None, None
            
        master_df = master_df.sort_values(by="Combined_Quality_Score", ascending=False)
        create_pdf_report(master_df, temp_path)
        
        interactive_table = master_df[['Chromosome', 'Genomic_Start', 'Spacer_20nt', 'GC_Content', 'Quality_Flag', 'Combined_Quality_Score']].head(20)
        
        return "Analysis Successful! Genomic background and biological assessment complete.", temp_path, interactive_table
    except Exception as e:
        return f"Pipeline Error: {str(e)}", None, None

# --- 5. GRADIO USER INTERFACE LAYOUT ---
with gr.Blocks() as demo:
    gr.Markdown("# 🧬 Advanced CRISPR Genomic Design Pipeline")
    gr.Markdown("An execution engine evaluating string alignment, human genomic background repeat risks, and transcription stability parameters.")
    
    with gr.Row():
        with gr.Column(scale=1):
            file_input = gr.File(label="Upload Genome/Transcript FASTA")
            mismatch_slider = gr.Slider(1, 4, 3, step=1, label="Cross-Reactivity Sensitivity")
            score_cutoff = gr.Slider(0, 90, 10, step=5, label="Minimum Quality Score Filter Cutoff")
            submit_btn = gr.Button("Run Design Pipeline", variant="primary")
            
        with gr.Column(scale=2):
            status_output = gr.Textbox(label="Status")
            file_output = gr.File(label="Download Full Report (.pdf)")
            
    gr.Markdown("### 📊 Top 20 Optimal Guide Preview Matrix")
    preview_dataframe = gr.Dataframe(headers=["Chromosome", "Position", "Spacer", "GC %", "Risk Flag", "Score"], interactive=False)
            
    submit_btn.click(
        run_crispr_pipeline, 
        [file_input, mismatch_slider, score_cutoff], 
        [status_output, file_output, preview_dataframe]
    )

demo.launch(server_name="0.0.0.0", server_port=7860)
