import os


def create_summary_table(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs):
    summary_tables = []
    for project in PROJECTS:
        for gender in ['female', 'male', 'female_male']:
            for cutoff in cutoffs:
                summary_tables.append(
                    os.path.join(
                        OUTPUT_PATH, project,
                        'DESeq2/DESeq2_input_table', DRUG_str, gender,
                        cutoff, 'summary_for_DESeq2.tsv'))
    return summary_tables
