import os


def create_deseq_output(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs):
    summary_tables = []
    for project in PROJECTS:
        for gender in ['female', 'male', 'female_male']:
            for cutoff in cutoffs:
                cutoff = f'cutoff_{str(cutoff)}'
                summary_tables.append(
                    os.path.join(
                        OUTPUT_PATH, project,
                        'DESeq2/DESeq2_output', DRUG_str, gender,
                        cutoff, 'DESeq2_results_summary.tsv'))
    return summary_tables
