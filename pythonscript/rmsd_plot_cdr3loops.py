# -*- coding: utf-8 -*-
"""
Original code:
Name: plot_benchmark_success_rates.py
Function: Script to plot T1, T5, T10, T20, T50, T100 for every case, colored by Capri criteria.
Date: 16-05-2023 14:52
Author: Yannick Aarts
"""

"""
Code was altered to make the success plots for T1, T2, T3, T4, T5 for only CDR3 loops using i-RMSD
Name: rmsd_plot_cdr3loops.py
Date: 06-01-2025

"""

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib import ticker
from matplotlib.patches import Circle, Rectangle
from matplotlib import rc
import matplotlib.patches as mpatches

import matplotlib.pyplot as plt
import numpy as np
import json

import pandas as pd

import plotly.express as px
import plotly.graph_objs as go

def main():
    lrmsd_f = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\rmsd_cdr3_B_good.tsv"
    success_plot_f = r"C:\Users\edsil\Documents\coevolutionproject\plots\af3_success_plot_cdr3B"
    success_rate_plot_f = r"C:\Users\edsil\Documents\coevolutionproject\plots\success_rate_plot_cdr3B"
    base_name = "AlphaFold3 CDR3 B"
    
    # Check if the correct number of arguments is provided no more or less
    if len(argv) != 7:
        print("Usage: python3 plot_benchmark_success_rates.py lrmsd.txt success_plot success_rate_plot base_name")

    lrmsd_dict = parse_results(lrmsd_f)
    rmsd_plot_name = success_plot_f
    success_plot_name = success_rate_plot_f

    labels = [1, 2, 3, 4, 5]

    color_matrix_all = eval_model_qualities(lrmsd_dict, labels, rmsd_plot_name)

    with open('data.txt', 'w') as f:
        for model in color_matrix_all:
            f.write(model + "\n")
            f.write(str(color_matrix_all[model]) + "\n")

    color_matrix_top = calc_max_values([1, 2, 3, 4, 5], color_matrix_all)

    color_matrix = make_rmsd_plot_all_at_one(color_matrix_top, labels, rmsd_plot_name, base_name)

    # Update colors for the plot (if required)
    color_matrix[color_matrix == 'blue'] = 'green'
    color_matrix[color_matrix == 'lightblue'] = 'lightgreen'
    color_matrix[color_matrix == 'darkblue'] = 'darkgreen'

    colors = ['darkgreen', 'green', 'lightgreen', 'lightgrey']

    color_counts, percenteages = calculate_color_percentages(color_matrix, labels)
    plot_success_rate(percenteages, success_plot_name, base_name)

def parse_results(datapath):
    """
    Parses a results file and returns a dictionary containing the data.

    Args:
        datapath: The path to the results file.

    Returns:
        A dictionary where the keys are model names and the values are lists of RMSD values.
    """
    data_dict = {}
    with open(datapath) as f:
        lines = f.readlines()
        for line in lines:
            splits = line.split()
            try:
                # Skip lines where the first value is not a valid model name
                model_name = splits[0].strip()
                
                # Skip non-numeric data in the rest of the line (headers, etc.)
                rmsd_values = [float(i) for i in splits[1:] if i.replace('.', '', 1).isdigit()]
                data_dict[model_name] = rmsd_values
            except ValueError:
                # Skip lines with any non-convertible data
                continue
    return data_dict

def calc_max_values(indices, data):
    """
    Calculates the maximum values for specified indices from a given data dictionary.

    Args:
        indices (list): A list of indices specifying the number of maximum values to consider.
        data (dict): A dictionary containing the data.

    Returns:
        dict: A dictionary where the keys are model names and the values are lists of maximum values.
    """
    result = {}
    for key, values in data.items():
        if key not in result:
            result[key] = []
        if len(values) > 0:
            for index in indices:
                result[key].append(max(values[:index]))
    return result

def eval_model_qualities(l_dict, labels, rmsd_plot_name="rmsd_plot"):
    """
    Evaluates model quality based on RMSD values.

    Args:
        l_dict (dict): A dictionary where the keys are model names and the values are RMSD values.
        labels (list): A list of labels representing model ranks.
        rmsd_plot_name (str): The name of the plot file.

    Returns:
        dict: A dictionary of model quality evaluations.
    """
    # Define thresholds and corresponding quality levels
    l_thresholds = [0, 1, 2, 4, float('inf')]  # Threshold values for color levels based on RMSD
    quality = [3, 2, 1, 0, -1]  # Model quality levels (3 = High, 2 = Medium, 1 = Acceptable, 0 = Incorrect)
    
    color_matrix_all = {}

    # Evaluate each model's RMSD values
    for pdb_id, rmsd_values in l_dict.items():
        print(f"Processing PDB ID: {pdb_id}")
        model_quality = []  # List to store quality scores for the current model
        
        for rmsd in rmsd_values:
            print(f"RMSD value: {rmsd}")
            
            # Assign quality based on RMSD value and thresholds
            if rmsd <= 1:
                model_quality.append(3)  # High quality
            elif rmsd <= 2:
                model_quality.append(2)  # Medium quality
            elif rmsd <= 4:
                model_quality.append(1)  # Acceptable quality
            else:
                model_quality.append(0)  # Incorrect quality (for RMSD > 4)
            
            print(f"Assigned quality: {model_quality[-1]}")

        # Add the model's quality scores to the dictionary
        color_matrix_all[pdb_id] = model_quality
        print(f"Final quality for {pdb_id}: {model_quality}\n")

    return color_matrix_all




def make_rmsd_plot_all_at_one(l_dict, labels, rmsd_plot_name="rmsd_plot", base_name="rmsd_plot"):
    # Define the categories
    unique_trav = ['7rm4', '8d5q', '8i5c']
    unique_trbv = ['7pbc', '7qpj', '8i5d']
    unique_mhc = ['7l1d', '7rrg', '8shi']
    unique_trav_trbv = ['7na5', '8wte', '8wul']

    # Remove duplicates between categories
    unique_trav = list(set(unique_trav) - set(unique_trav_trbv))  # Remove TRAV models in TRAV+TRBV
    unique_trbv = list(set(unique_trbv) - set(unique_trav_trbv))  # Remove TRBV models in TRAV+TRBV
    unique_trav_trbv = list(set(unique_trav_trbv))  # TRAV and TRBV models

    # Collect all models in the order you desire
    all_models = []
    all_models.extend(unique_trav)      # Add unique TRAV models
    all_models.extend(unique_trbv)      # Add unique TRBV models
    all_models.extend(unique_trav_trbv) # Add unique TRAV+TRBV models
    all_models.extend(unique_mhc)       # Add MHC models

    # Collect remaining models that are not in special categories
    medium_models = [model for model in l_dict if model not in all_models]
    all_models.extend(medium_models)   # Add medium models

    # Ensure medium models are added in the correct order
    # This assumes `medium_models` can be sorted in the way you want, 
    # if you have a specific order in mind, you can sort it here:
    # For example, you can sort medium models alphabetically or in any custom order
    all_models = unique_trav + unique_trbv + unique_trav_trbv + unique_mhc + sorted(medium_models)

    # Color mappings for each category based on RMSD value
    medium_colors = ['lightgrey', 'lightgreen', 'green', 'darkgreen']
    rigid_colors = ['lightgrey', 'lightblue', 'blue', 'darkblue']

    # Prepare x and y data for plotting
    x = np.arange(len(all_models))
    y = np.arange(len(labels))

    # Set square size and spacing
    square_size = 0.8
    spacing = 0.2

    # Create plot
    fig, ax = plt.subplots()

    # Model count
    model_num = len(all_models)
    colors = np.zeros((model_num, 6), dtype=object)  # Assuming there can be up to 6 RMSD values per model

    # Loop through each model and plot its RMSD values with the appropriate color
    for i, model in enumerate(all_models):
        if model in l_dict:
            rmsd_values = l_dict[model]

            for j, rmsd in enumerate(rmsd_values):
                # Default color to white (if no other condition is met)
                color = 'white'

                # Assign color based on the category the model belongs to
                if model in unique_trav_trbv:
                    color = medium_colors[rmsd]  # Clamp the RMSD index
                elif model in unique_trav:
                    color = medium_colors[rmsd]  # Clamp the RMSD index
                elif model in unique_trbv:
                    color = medium_colors[rmsd]  # Clamp the RMSD index
                elif model in unique_mhc:
                    color = medium_colors[rmsd]  # Clamp the RMSD index
                elif model in medium_models:
                    color = medium_colors[rmsd]  # Clamp the RMSD index

                if color:
                    colors[i, j] = color
                    if color == 'darkgreen':
                        color = '#003600'
                    if color == 'darkblue':
                        color = 'midnightblue'

                # Assign the color to the respective square
                rect = plt.Rectangle((i - square_size / 2, j - square_size / 2), square_size, square_size, color=color)
                ax.add_patch(rect)

    # Set plot margins and aspect ratio
    plt.margins(0.2)
    plt.gca().set_aspect('equal', adjustable='box')

    # Set x and y axis limits
    plt.xlim(-0.5, model_num - 0.5)
    plt.ylim(-0.5, len(labels) - 0.5)

    # Set x and y axis labels
    ax.set_xticks(np.arange(model_num))
    ax.set_xticklabels([model.upper() for model in all_models], rotation=90, fontsize=8, fontname='DejaVu Sans Mono')

    # Set y-axis tick labels
    ax.set_yticks(np.arange(len(labels)) - 0.2)
    ax.set_yticklabels(['Top ' + str(label) for label in labels], fontsize=8, fontname='Helvetica')
    ax.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False)

    # Set plot title
    ax.set_title(base_name, fontsize=10, fontweight='bold')

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Adjust the spacing between the plot and the legend
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.45, top=0.8)

    # Create custom legend for different RMSD color categories
    medium_patch_incorrect = mpatches.Patch(color='lightgrey', label='Incorrect')
    medium_patch_acceptable = mpatches.Patch(color='lightgreen', label='Acceptable')
    medium_patch_medium = mpatches.Patch(color='green', label='Medium')
    medium_patch_high = mpatches.Patch(color='darkgreen', label='High')

    # Place legend below the plot
    ax.legend(handles=[medium_patch_incorrect, medium_patch_acceptable, medium_patch_medium, medium_patch_high],
              loc='upper center', bbox_to_anchor=(0.5, -0.55),  # Move the legend further down
              ncol=3, fontsize=8, frameon=True, fancybox=True)

    # Save the plot
    plt.savefig(rmsd_plot_name+'.png', bbox_inches='tight', dpi=300)
    plt.savefig(rmsd_plot_name+'.pdf', bbox_inches='tight', dpi=300)

    return colors




def calculate_color_percentages(color_matrix, labels):
    percentages = {}
    color_counts = np.zeros((6, 4), dtype=float)
    for i, label in enumerate(labels):
        counts = {'darkgreen': 0, 'green': 0, 'lightgreen': 0, 'lightgrey': 0}

        for model_colors in color_matrix:
            counts[model_colors[i]] += 1
        total = sum(counts.values())
        percentages[label] = {color: (count / total) * 100 for color, count in counts.items()}
        color_counts[i] = list(percentages[label].values())
    return color_counts, percentages

def reorder_dict_keys(dictionary, keys_order):
    return {key: dictionary[key] for key in keys_order}

def plot_success_rate(data, success_plot_file, base_name):
    for key in data:
        data[key] = reorder_dict_keys(data[key], ['darkgreen', 'green', 'lightgreen', 'lightgrey'])

    df = pd.DataFrame(data).T.reset_index().rename(columns={'index': 'Top'})
    
    df = df.melt(id_vars=['Top'], var_name='color', value_name='percentage')
    df['Top'] = df['Top'].astype(str)

    desired_order = ['#003600', 'green', 'lightgreen', 'lightgrey']
    df['color'] = df['color'].replace('darkgreen', '#003600')
    df['color'] = pd.Categorical(df['color'], categories=desired_order, ordered=True)

    fig = px.bar(df, x="Top", y="percentage", color="color", title=str(base_name) + " Success Rate",
                 color_discrete_map={color: color for color in desired_order})

    reference_lines = [20, 40, 60, 80, 100]
    for line in reference_lines:
        fig.add_hline(y=line, line_dash="dash", line_color="black")

    fig.update_layout(bargap=0.1, font=dict(size=60),
                      xaxis=dict(tickmode='array', ticktext=['Top ' + str(val) for val in df['Top']]),
                      title_x=0.5, title=dict(y=0.98), legend_title_text='Model Quality')

    legend_labels = {'#003600': 'High', 'green': 'Medium', 'lightgreen': 'Acceptable', 'lightgrey': 'Incorrect'}
    
    for color, label in legend_labels.items():
        for trace in fig.data:
            if trace.marker.color == color:
                trace.name = label

    for trace in fig.data:
        label = trace.name
        trace.legendgroup = label
        trace.hovertemplate = trace.hovertemplate.replace(trace.name, label)

    for trace in fig.data:
        trace.hovertext = [f"{label}: {y:.2f}%" for label, y in zip(df['Top'], trace.y)]

    fig.write_html(success_plot_file + ".html")

if __name__ == "__main__":
    main()
