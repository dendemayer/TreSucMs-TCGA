import os


def validation_plots(lifeline_plots):
    final_list = []
    for i in lifeline_plots:
        DOWN_val = os.path.join(os.path.split(i)[0], f"{os.path.split(i)[1].split('plot')[0]}DOWN_val_plot{os.path.split(i)[1].split('plot')[1]}")
        final_list.append(DOWN_val)
        UP_val = os.path.join(os.path.split(i)[0], f"{os.path.split(i)[1].split('plot')[0]}UP_val_plot{os.path.split(i)[1].split('plot')[1]}")
        final_list.append(UP_val)
    return final_list
