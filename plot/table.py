import matplotlib
import pandas as pd
from matplotlib import pyplot as plt
from plottable import ColumnDefinition, Table
from plottable.plots import bar
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
import numpy as np

def add_colorbars(fig, spec, n_cols, n_rows, colorbar_specs, score_transfunc):
    """Creates a group of colorbars that is centered below the table."""
    bar_width_fraction = 0.12
    padding_fraction = 0.04

    # Total width of all colorbars and the gaps between them
    total_width_fraction = len(colorbar_specs) * bar_width_fraction + (len(colorbar_specs) - 1) * padding_fraction

    # Calculate starting column to center the group
    start_triplet_fraction = (1 - total_width_fraction) / 2
    start_triplet_col = int(start_triplet_fraction * n_cols)

    # Calculate the width of each colorbar and spacing in columns
    bar_width_cols = int(bar_width_fraction * n_cols)
    padding_cols = int(padding_fraction * n_cols)

    for i, colorbar_spec in enumerate(colorbar_specs):
        cmap = colorbar_spec["cmap"]
        label = colorbar_spec["label"]

        norm = Normalize(vmin=0, vmax=1) # if score_transfunc == "" else Normalize(vmin=n_rows, vmax=1)
        scalar_mappable = ScalarMappable(norm=norm, cmap=cmap)

        # Compute the exact start and end columns for each colorbar
        start_col = (start_triplet_col + i * (bar_width_cols + padding_cols)) + 2
        end_col = start_col + bar_width_cols

        # Add the colorbar
        ax_colorbar = fig.add_subplot(spec[1, slice(start_col, end_col)])
        cbar = fig.colorbar(
            scalar_mappable,
            cax=ax_colorbar,
            orientation="horizontal",
        )
        # if score_transfunc == "_scaled":
        #     cbar.ax.invert_xaxis()

        # Add title above and numbers (ticks) below
        ax_colorbar.set_title(label, pad=8, fontsize=10)
        cbar.ax.xaxis.set_ticks_position("bottom")
        cbar.ax.xaxis.set_label_position("bottom")
        ax_colorbar.xaxis.label.set_size(8)


def make_discrete_cmap(cmap_name, num_colors):
    """Creates a discrete colormap with a specified number of colors."""
    base_cmap = plt.get_cmap(cmap_name)
    color_list = base_cmap(np.linspace(0, 1, num_colors))
    return mcolors.ListedColormap(color_list)


def white_yellow_green_cm():
    lut_size = 256
    spec = [
        (1.0, 1.0, 1.0),
        (0.90196078431372551, 0.96078431372549022, 0.81568627450980391),
        (0.30196078431372547, 0.5725490196078431, 0.12941176470588237),
    ]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("wYlGn", spec, lut_size)
    return cmap


def get_scalar_mapppable(col_data, norm_type=None):
    if norm_type == "minmax":
        vmin = col_data.min()
        vmax = col_data.max()
    if norm_type == "interquartile":
        # taken from plottable.cmap.normed_cmap
        num_stds = 2.5
        _median, _std = col_data.median(), col_data.std()
        vmin = _median - num_stds * _std
        vmax = _median + num_stds * _std
    else:
        vmin, vmax = 0, 1

    cmap = white_yellow_green_cm()
    norm = matplotlib.colors.Normalize(vmin, vmax)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    return m


def draw_table(pivot_path: str, ax: plt.Axes, cbar_specs, score_transfunc):
    """
    Adapted from:
    https://github.com/yoseflab/scib-metrics/blob/0.4.1/src/scib_metrics/benchmark/_core.py#L276-L364
    """
    df = pd.read_parquet(pivot_path)

    # extract cmaps for the respective metric categories
    batch_cmap = next(spec["cmap"] for spec in cbar_specs if spec["label"] == "Batch correction metrics")
    bio_cmap = next(spec["cmap"] for spec in cbar_specs if spec["label"] == "Bio conservation metrics")
    aggregate_cmap = next(spec["cmap"] for spec in cbar_specs if spec["label"] == "Aggregate scores")

    # extract labels from columns so we can construct queries later
    eval_keys = df.columns.get_level_values("eval_key").unique()
    metric_type_keys = df.columns.get_level_values("metric_type").unique()
    metric_keys = df.columns.get_level_values("metric").unique().drop(["mean_batch", "mean_bio", "mean_overall"])

    # flatten the MultiIndex so we can properly create a plottable
    df_for_table = df.reset_index()
    df_for_table.columns = [
        "_".join(filter(None, map(str, col))).strip() if isinstance(col, tuple) else col
        for col in df_for_table.columns.values
    ]

    # Add first column for the method names
    column_definitions = [
        ColumnDefinition(
            "Method", width=1.5, textprops={"ha": "left", "weight": "bold"}
        ),
    ]

    # define how circles and barplots look
    circle_props = {"ha": "center", "bbox": {"boxstyle": "circle", "pad": 0.25}}

    # we'll use this to resort the df columns later
    final_col_order = []

    # iterate over keys and generate columns for each
    for eval_key in eval_keys:
        # we'll iterate over this later for the aggregated scores
        mean_col_fstring_tuples = []

        for metric_type in metric_type_keys:
            # Skip 'aggregate_score' as it contains mean values
            if metric_type == 'aggregate_score':
                continue

            for metric in metric_keys:
                colname = f"{eval_key}_{metric_type}_{metric}"

                # we might not have all metrics for all eval_keys (esp during testing)
                if colname in df_for_table.columns:
                    col_def = ColumnDefinition(
                        colname,
                        title=metric.replace("_", "\n", 1),
                        textprops=circle_props,
                        width=1,
                        cmap=batch_cmap if metric_type == "batch_correction" else bio_cmap,
                        group=eval_key,
                        formatter="{:.2f}",
                    )
                    column_definitions.append(col_def)
                    final_col_order.append(colname)

            # Collect mean columns
            metric_substring = metric_type.split("_")[0]
            mean_colname = f"{eval_key}_aggregate_score_mean_{metric_substring}"
            if mean_colname in df_for_table.columns:
                mean_col_fstring_tuples.append((eval_key, 'aggregate_score', f"mean_{metric_substring}"))

        # Now, append 'mean_overall'
        mean_overall_colname = f"{eval_key}_aggregate_score_mean_overall"
        if mean_overall_colname in df_for_table.columns:
            mean_col_fstring_tuples.append((eval_key, 'aggregate_score', "mean_overall"))

        title_lookup = {
            "mean_batch": "mean\nbatch",
            "mean_bio": "mean\nbio",
            "mean_overall": f"mean\n{eval_key.replace('Metadata_', '')}",
        }

        border_lookup = {
            0: "left",
            len(mean_col_fstring_tuples) - 1: "right",
        }

        for i, (ek, mt, meancol) in enumerate(mean_col_fstring_tuples):
            colname = f"{ek}_{mt}_{meancol}"
            col_def = ColumnDefinition(
                colname,
                title=title_lookup.get(meancol, meancol).replace("_", "\n", 1),
                textprops=circle_props,
                width=1,
                cmap=aggregate_cmap,
                group=eval_key,
                formatter="{:.2f}",
                border=border_lookup.get(i),
            )
            column_definitions.append(col_def)
            final_col_order.append(colname)

    name_for_col_that_averages_across_eval_keys = "total_mean"
    non_aggregate_cols = [col for col in df_for_table.columns if "aggregate_score" not in col and col != "Method"]

    if score_transfunc == "":

        pass
    
    elif score_transfunc == "_scaled":
        
        df_for_table[non_aggregate_cols] = df_for_table[non_aggregate_cols].apply(lambda x: (x - x.min()) / (x.max() - x.min()))
    
    elif score_transfunc == "_ranked":
        
        # Identify non-aggregate columns and substitue scores with ranks, then scale to [0, 1]
        df_for_table[non_aggregate_cols] = df_for_table[non_aggregate_cols].rank(method="average", ascending=True)
        df_for_table[non_aggregate_cols] = df_for_table[non_aggregate_cols].apply(lambda x: (x - x.min()) / (x.max() - x.min()))

    # re-average per key and metric type
    for eval_key in eval_keys:

        metric_types = [mt for mt in metric_type_keys if mt != 'aggregate_score']
        
        for metric_type in metric_types:

            metric_cols = [col for col in non_aggregate_cols if col.startswith(f"{eval_key}_{metric_type}_")]

            if not metric_cols:
                continue  # Skip if no metrics found for this category
            
            # mean rank for metric_type (batch / bio) 
            mean_rank_col = f"{eval_key}_aggregate_score_mean_{metric_type.split('_')[0]}"
            df_for_table[mean_rank_col] = df_for_table[metric_cols].mean(axis=1)

        # mean rank for eval_key
        overall_metric_cols = [col for col in non_aggregate_cols if col.startswith(f"{eval_key}_")]
        mean_overall_col = f"{eval_key}_aggregate_score_mean_overall"
        if mean_overall_col in df_for_table.columns:
            df_for_table[mean_overall_col] = df_for_table[overall_metric_cols].mean(axis=1)

    # overall mean rank
    total_mean = df_for_table[[col for col in df_for_table.columns if "aggregate_score_mean_overall" in col]].mean(axis=1)
    df_for_table[name_for_col_that_averages_across_eval_keys] = total_mean


    col_def = ColumnDefinition(
        name_for_col_that_averages_across_eval_keys,
        title="Total",
        textprops=circle_props,
        width=1,
        cmap=aggregate_cmap,
        formatter="{:.2f}",
        border="both",
    )
    column_definitions.append(col_def)
    final_col_order.append(name_for_col_that_averages_across_eval_keys)

    # Reorder columns and rows
    df_for_table = df_for_table[["Method"] + final_col_order]
    df_for_table = df_for_table.sort_values(by=name_for_col_that_averages_across_eval_keys, ascending=False)

    # create the table
    plt.style.use("default")
    tab = Table(
        df_for_table,
        cell_kw={
            "linewidth": 0,
            "edgecolor": "k",
        },
        column_definitions=column_definitions,
        ax=ax,
        row_dividers=True,
        footer_divider=True,
        textprops={"fontsize": 10, "ha": "center"},
        row_divider_kw={"linewidth": 1, "linestyle": (0, (1, 5))},
        col_label_divider_kw={"linewidth": 1, "linestyle": "-"},
        column_border_kw={"linewidth": 1, "linestyle": "-"},
        index_col="Method",
    )
    # Adjust colnames for font color settings
    colnames = [col_def.name for col_def in column_definitions[1:]]  # Exclude "Method"
    tab.autoset_fontcolors(colnames=colnames)