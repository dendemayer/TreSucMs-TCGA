import os


def aggregate_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs,
                             threshold, pipeline):
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
                metilene_path = os.path.join(
                    OUTPUT_PATH, project, f'{pipeline}/{pipeline}_output',
                    DRUG_str, gender, cutoff)
                file_path = os.path.join(
                    metilene_path, threshold_str,
                    f'{pipeline}_lifelines_aggregated.tsv.gz')
                aggregate_lifeline_plots_list.append(file_path)

    return aggregate_lifeline_plots_list
