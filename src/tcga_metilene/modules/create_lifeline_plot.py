import re
import pandas as pd
import os


def create_lifeline_plots(metilene_plots, threshold):
    """
    for each metilene plot, we create the corresponding lifeline plot, which is
    based on the given region
    out of the region, the best p_value is assessed with a logrank_test
    at this position, the lifeline plot is created. on the basis of the beta
    value median, the cases are separated in groups UP and DOWN
    """
    # reduce the set to the lineplots naming, the boxplot is not needed,
    metilene_plots = [i for i in metilene_plots if re.search('_lineplot_median_beta_value_', i)]
    lifeline_plots = [i.replace('lineplot_median_beta_value', 'lifeline_plot') for i in metilene_plots]
    lifeline_t = []
    for t in threshold:
        t = f'threshold_{str(t)}'
        for path in lifeline_plots:
            lifeline_t.append(os.path.join(os.path.split(path)[0], t, os.path.split(path)[1]))
    return lifeline_t
