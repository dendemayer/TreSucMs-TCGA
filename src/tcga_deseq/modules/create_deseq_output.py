import os


def create_deseq_output(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs):
    summary_tables = []
    for project in PROJECTS:
        for gender in ['female', 'male', 'female_male']:
            for cutoff in cutoffs:
                deseq_path = os.path.join(
                    OUTPUT_PATH, project, 'DESeq2/DESeq2_output', DRUG_str,
                    gender, cutoff)
                for count_type in ['norm', 'raw', 'nt', 'vsd']:
                    for in_de in ['IN', 'DE']:
                        summary_tables.append(
                            os.path.join(
                                deseq_path, "DESeq2_heatmap_log2f" +
                                f"{in_de}CREASE_{count_type}_counts.tsv",))
    return summary_tables

# summary_tables.append(
#     os.path.join(
#         OUTPUT_PATH, project,
#         'DESeq2/DESeq2_output', DRUG_str, gender,
#         cutoff, 'DESeq2_results_summary.tsv'))
