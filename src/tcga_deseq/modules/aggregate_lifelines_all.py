import os


def aggregate_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs,
                             threshold):
    """
    # iterating over all thresholds not necessary, they are summarized into one
    # dir like: threshold_0_threshold_5_threshold_10
    """
    threshold_str = '_'.join([f'threshold_{str(i)}' for i in threshold])
    aggregate_lifeline_plots_list = []
    # do not request any aggregation for which the data is missing, f.e. male
    # for CESC:
    for project in PROJECTS:
        for gender in ['female', 'female_male', 'male']:
            for cutoff in cutoffs:
                cutoff = f'cutoff_{str(cutoff)}'
                deseq_path = os.path.join(
                    OUTPUT_PATH, project, 'DESeq2/DESeq2_output', DRUG_str,
                    gender, cutoff)
                file_path = os.path.join(
                    deseq_path, threshold_str,
                    'DESeq2_lifelines_aggregated.tsv.gz')
                aggregate_lifeline_plots_list.append(file_path)

    return aggregate_lifeline_plots_list

# not needed, those names can be derived from the aggregated lifeline plots
# def evaluate_lifelines_all(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs,
#                              threshold):
#     threshold_str = '_'.join([f'threshold_{str(i)}' for i in threshold])
#     evaluate_lifeline_list = []
#     # do not request any aggregation for which the data is missing, f.e. male
#     # for CESC:
#     for project in PROJECTS:
#         for gender in ['female', 'female_male', 'male']:
#             for cutoff in cutoffs:
#                 cutoff = f'cutoff_{str(cutoff)}'
#                 deseq_path = os.path.join(
#                     OUTPUT_PATH, project, 'DESeq2/DESeq2_output', DRUG_str,
#                     gender, cutoff)
#                 file_path_tsv = os.path.join(
#                     deseq_path, threshold_str,
#                     f'DESeq2_lifelines_evaluated.tsv.gz')
#                 file_path_pdf = os.path.join(
#                     deseq_path, threshold_str,
#                     f'DESeq2_lifelines_evaluated.pdf')
#                 evaluate_lifeline_list.append(file_path_tsv)
#                 evaluate_lifeline_list.append(file_path_pdf)

#     return evaluate_lifeline_list
