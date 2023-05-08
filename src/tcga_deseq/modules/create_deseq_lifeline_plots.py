import os
import pandas as pd


def create_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs, threshold):
    lifeline_plots = []
    for project in PROJECTS:
        for gender in ['female', 'male', 'female_male']:
            for cutoff in cutoffs:
                cutoff = f'cutoff_{str(cutoff)}'
                deseq_path = os.path.join( OUTPUT_PATH, project, 'DESeq2/DESeq2_output', DRUG_str, gender, cutoff)
                for thresh in threshold:
                    thr = f'threshold_{str(thresh)}'
                    for count_type in ['norm', 'raw', 'nt']:
                        for in_de in ['IN', 'DE']:
                            try:
                                ENSGs = pd.read_table(os.path.join(deseq_path, f'DESeq2_heatmap_log2f{in_de}CREASE_{count_type}_counts.tsv')).index.tolist()
                            except Exception as e:
                                continue
                            for ENSG in ENSGs:
                                lifeline_plots.append(os.path.join(deseq_path, thr, f'DESeq2_log2f_{in_de}CREASE_{count_type}_{ENSG}_lifeline.pdf'))
    return lifeline_plots

def create_lifeline_plots_validation(lifeline_plots):
    return [i.replace('.pdf', '_UP_val.pdf') for i in lifeline_plots]

