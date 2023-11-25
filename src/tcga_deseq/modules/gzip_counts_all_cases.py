import os


def create_gz_counts(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs, count_types):
    gz_counts = []
    for project in PROJECTS:
        for gender in ['female', 'male', 'female_male']:
            for cutoff in cutoffs:
                deseq_path = os.path.join(
                    OUTPUT_PATH, project, 'DESeq2/DESeq2_output', DRUG_str,
                    gender, cutoff)
                for count_type in count_types:
                    gz_counts.append(
                        os.path.join(
                            deseq_path,
                            f"DESeq2_{count_type}s_all_cases.tsv.gz"))
    return gz_counts
