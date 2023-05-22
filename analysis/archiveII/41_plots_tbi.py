#!/usr/bin/env python3

import argparse
import json
import statistics
import matplotlib.pyplot as plt
from scipy.stats import shapiro, kstest, ttest_ind, wilcoxon
import numpy as np
import ternary
import random
import math
import pandas as pd
import re 
import plotly.graph_objects as go
import scipy.stats as stat

TIMESTAMPS = ["0.04", "1.0", "10.0", "60.0", "600.0", "3600.0", "equilibrium", "eq"]
STRUCTURES = ["nat", "natc", "mfe", "eq", "nat_rev", "natc_rev", "mfe_rev", "eq_rev"]
MEASUREMENTS = ["TP", "FP", "FN", "TN", "TPR", "PPV", "MCC", "F1", "bpdis"]

def ternary_plot(input_file):
    with open(input_file) as f:
        db = json.load(f)

    points = []
    sums = []
    for header in db:
        x = db[header]["1.0"]["nat"]["bpdis"]
        y = db[header]["1.0"]["nat_rev"]["bpdis"]
        z = db[header][f"bpdis_nat_subopt"]  # TODO db[header]["eq"]["nat"]["MCC"]
        sum = x + y + z
        x /= sum
        y /= sum
        z /= sum
        sums.append(sum)
        points.append((x,y,z))
    max_sum = max(sums)
    sums_normalized = [x/max_sum for x in sums]
    # Scatter Plot
    figure, tax = ternary.figure(scale=1)
    # figure.set_size_inches(10, 10) TODO: triangle not equilateral if heatmap is shown and eight=width
    tax.scatter(points, s=20, c=sums_normalized, label="RNA", cmap='plasma', alpha=0.6, linewidths=0)
    # tax.legend()
    
    tax.set_title("Bp-distances to native structure", fontsize=16)
    tax.boundary(linewidth=2.0)
    tax.gridlines(multiple=0.1, color="black")
    tax.ticks(axis='lbr', linewidth=1, multiple=0.1, tick_formats="%.1f", offset=0.02, fontsize=8)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    offset=0.15
    fontsize=12
    tax.left_axis_label("Equilibrium", fontsize=fontsize, offset=offset)
    tax.right_axis_label("Reversed", fontsize=fontsize, offset=offset)
    tax.bottom_axis_label("Normal", fontsize=fontsize, offset=offset)
    tax.heatmap(data={x:y for x,y in zip(points, sums)}, cmap='plasma', cbarlabel="Total distance")

    tax.savefig("tRNA_ternary.png")

def plots(directory, rna_classes=["5s", "16s", "23s", "grpI", "grpII", "RNaseP", "srp", "telomerase", "tmRNA", "tRNA"], TS="1.0", MS="MCC"):
    d = {}
    for rna_class in rna_classes:
        with open(f"../{directory}/analysis_results/{rna_class}/{rna_class}.JSON") as f:
            db = json.load(f)
        d[rna_class] = {ts: {st: {ms: {"values": [], "mean": None, "stdev": None, "median": None} for ms in MEASUREMENTS} for st in STRUCTURES} for ts in TIMESTAMPS}  # stores statistics

        for header in db:
                nan_flag = False
                d_temp = {}
                d_temp[rna_class] = {ts: {st: {ms: None for ms in MEASUREMENTS} for st in STRUCTURES} for ts in TIMESTAMPS}
                for ts in TIMESTAMPS:
                    for st in STRUCTURES:
                        for ms in MEASUREMENTS:
                            try:
                                val = float(db[header][ts][st][ms])
                                if math.isnan(val):
                                    # val = 0  # TODO: !!!
                                    if ms == MS:
                                        # print("Removing seq:", header, ts, st, ms)
                                        nan_flag = True
                                d_temp[rna_class][ts][st][ms] = val
                            except KeyError:  # for reversed-st in "eq"
                                pass
                if not nan_flag:
                    for ts in TIMESTAMPS:
                        for st in STRUCTURES:
                            for ms in MEASUREMENTS:
                                d[rna_class][ts][st][ms]["values"].append(d_temp[rna_class][ts][st][ms])
                else:
                    print("Seq with nan in 'MS'", header)
            
        for ts in TIMESTAMPS:
            for st in STRUCTURES:
                for ms in MEASUREMENTS:
                    try:
                        d[rna_class][ts][st][ms]["mean"] = statistics.mean(d[rna_class][ts][st][ms]["values"])
                        d[rna_class][ts][st][ms]["stdev"] = statistics.stdev(d[rna_class][ts][st][ms]["values"])
                        d[rna_class][ts][st][ms]["median"] = statistics.median(d[rna_class][ts][st][ms]["values"])
                    except TypeError: # removed because includes nan
                        pass
                    except statistics.StatisticsError:  # for reversed-st in "eq"
                        pass

        # ts = "1.0"
        # ms = "MCC"
        # t1 = d[ts]["nat"][ms]["values"]
        # t2 = d[ts]["nat_rev"][ms]["values"]
        
        # # boxplot
        # data = [t1, t2]
        # plt.boxplot(data)
        # plt.title("Boxplot - MCC, 1.0s")
        # plt.xlabel("Transcription direction")
        # plt.ylabel("MCC")
        # plt.xticks([1, 2], ['Normal', 'Reversed'])
        # # show plot
        # plt.savefig("tRNA_boxplot.png")
        # plt.close()

        # # histogram
        # plt.hist(t1, bins=50, alpha=0.5, label="normal")
        # plt.hist(t2, bins=50, alpha=0.5, label="reversed")
        # plt.title("Histogram - MCC, 1.0s")
        # plt.xlabel(ms)
        # plt.ylabel("Frequency")
        # plt.legend()
        # plt.savefig("tRNA_hist.png")
        # plt.close()

        # mean over different timestamps
        plt.plot([d[rna_class][ts]["nat"][MS]["mean"] for ts in TIMESTAMPS[:-2]], label="normal")
        plt.plot([d[rna_class][ts]["nat_rev"][MS]["mean"] for ts in TIMESTAMPS[:-2]], label="reversed")
        plt.title(f"Time series - {MS} mean")
        plt.xlabel("Time [sec]")
        plt.xticks(np.arange(len(TIMESTAMPS[:-2])), TIMESTAMPS[:-2])
        plt.ylabel(f"{MS}")
        plt.legend()
        plt.savefig(f"../{directory}/plots/time_series/time_series_{MS}_{rna_class}.png")
        plt.close()

        # ternary_plot(f"../{directory}/analysis_results/{rna_class}/{rna_class}.JSON")


def combined_plots(directory, rna_classes=["5s", "16s", "23s", "grpI", "grpII", "RNaseP", "srp", "telomerase", "tmRNA", "tRNA"], TS="1.0", MS="MCC"):
    d = {}
    for rna_class in rna_classes:
        with open(f"../{directory}/analysis_results/{rna_class}/{rna_class}.JSON") as f:
            db = json.load(f)

        d[rna_class] = {ts: {st: {ms: {"values": [], "mean": None, "stdev": None, "median": None} for ms in MEASUREMENTS} for st in STRUCTURES} for ts in TIMESTAMPS}  # stores statistics

        for header in db:
            nan_flag = False
            d_temp = {}
            d_temp[rna_class] = {ts: {st: {ms: None for ms in MEASUREMENTS} for st in STRUCTURES} for ts in TIMESTAMPS}
            for ts in TIMESTAMPS:
                for st in STRUCTURES:
                    for ms in MEASUREMENTS:
                        try:
                            val = float(db[header][ts][st][ms])
                            if math.isnan(val):
                                val = 0  # TODO: !!!
                                if ms == MS:
                                    # print("Removing seq:", header, ts, st, ms)
                                    nan_flag = True
                            d_temp[rna_class][ts][st][ms] = val
                        except KeyError:  # for reversed-st in "eq"
                            pass
            if not nan_flag:
                for ts in TIMESTAMPS:
                    for st in STRUCTURES:
                        for ms in MEASUREMENTS:
                            d[rna_class][ts][st][ms]["values"].append(d_temp[rna_class][ts][st][ms])
            else:
                print("Seq with nan in 'MS'", header)
        
        for ts in TIMESTAMPS:
            for st in STRUCTURES:
                for ms in MEASUREMENTS:
                    try:
                        d[rna_class][ts][st][ms]["mean"] = statistics.mean(d[rna_class][ts][st][ms]["values"])
                        d[rna_class][ts][st][ms]["stdev"] = statistics.stdev(d[rna_class][ts][st][ms]["values"])
                        d[rna_class][ts][st][ms]["median"] = statistics.median(d[rna_class][ts][st][ms]["values"])
                    except TypeError: # removed because includes nan
                        pass
                    except statistics.StatisticsError:  # for reversed-st in "eq"
                        pass

    # create lists of plotted values
    data_nat = [d[rna_class][TS]["nat"][MS]["values"] for rna_class in rna_classes]
    data_nat_rev = [d[rna_class][TS]["nat_rev"][MS]["values"] for rna_class in rna_classes]
    data_diff = []
    for i in range(len(data_nat)):
        data_diff.append([x-y for x,y in zip(data_nat[i], data_nat_rev[i])])
    data_eq = [d[rna_class]["eq"]["nat"][MS]["values"] for rna_class in rna_classes]

    tick_names = rna_classes.copy()
    if "telomerase" in tick_names:
        tick_names[tick_names.index("telomerase")] = "telo."

    # boxplot (without equ), https://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    ticks = tick_names

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    plt.figure()

    bpl = plt.boxplot(data_nat, positions=np.array(range(len(data_nat)))*2-0.34, sym='', widths=0.55)
    bpr = plt.boxplot(data_nat_rev, positions=np.array(range(len(data_nat_rev)))*2+0.34, sym='', widths=0.55)
    set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
    set_box_color(bpr, '#2C7BB6')

    # draw temporary red and blue lines and use them to create a legend
    plt.plot([], c='#D7191C', label='normal')
    plt.plot([], c='#2C7BB6', label='reversed')
    plt.legend()

    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    # plt.xlim(-2, len(ticks)*2)
    # plt.ylim(0, 8)
    # plt.tight_layout()
    plt.title(f"Boxplot\n{TS}s")
    plt.ylabel(f"{MS}")
    plt.savefig(f'../{directory}/plots/boxplot/boxplot_{MS}_{TS}s.png')
    plt.close()
    
    # boxplot, https://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    ticks = tick_names

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    plt.figure(figsize=(8,4))

    bpl = plt.boxplot(data_nat, positions=np.array(range(len(data_nat)))*3-0.72, sym='', widths=0.6)
    bpm = plt.boxplot(data_nat_rev, positions=np.array(range(len(data_nat_rev)))*3+0.0, sym='', widths=0.6)
    bpr = plt.boxplot(data_eq, positions=np.array(range(len(data_eq)))*3+0.72, sym='', widths=0.6)
    set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
    set_box_color(bpm, '#2C7BB6')
    set_box_color(bpr, '#000')

    # draw temporary red and blue lines and use them to create a legend
    plt.plot([], c='#D7191C', label='normal')
    plt.plot([], c='#2C7BB6', label='reversed')
    plt.plot([], c='#000', label='equilibrium')
    plt.legend()

    plt.xticks(range(0, len(ticks) * 3, 3), ticks)
    plt.xlim(-3, len(ticks)*3)
    # plt.ylim(0, 8)
    # plt.tight_layout()
    plt.title(f"Boxplot\n{TS}s")
    plt.ylabel(f"{MS}")
    plt.savefig(f'../{directory}/plots/boxplot_with_equ/boxplot_with_equ_{MS}_{TS}s.png')
    plt.close()

    # difference_boxplot
    ticks = tick_names

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    plt.figure()

    plt.axhline(y = 0, color = 'black', linestyle = 'dashed', linewidth=1)  
    bpl = plt.boxplot(data_diff, positions=np.array(range(len(data_nat)))*2.0-0.0, sym='', widths=0.6, showmeans=True, meanprops=dict(marker='o', markeredgecolor='blue', markerfacecolor='blue', markersize=3))
    plt.scatter([], [], marker='o', c='blue', s=10, label='mean')
    plt.legend()
    # bpr = plt.boxplot(data_nat_rev, positions=np.array(range(len(data_nat_rev)))*2.0+0.4, sym='', widths=0.6)
    set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
    # set_box_color(bpr, '#2C7BB6')


    # draw temporary red and blue lines and use them to create a legend
    # plt.plot([], c='#D7191C', label='normal - reversed')
    # plt.plot([], c='#2C7BB6', label='nat_rev')
    # plt.legend()

    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-2, len(ticks)*2)
    # plt.ylim(0, 8)
    # plt.tight_layout()
    plt.title(f"Boxplot: Difference per RNA (normal - reversed)\n{TS}s")
    plt.ylabel(f"{MS} difference per RNA (normal - reversed)")
    plt.savefig(f'../{directory}/plots/difference_boxplot/difference_boxplot_{MS}_{TS}s.png')
    plt.close()


    # volcano, https://thecodingbiologist.com/posts/Making-volcano-plots-in-python-in-Google-Colab
    # two-sided wilcoxon test
    pvals = []
    mean_nat = []
    mean_nat_rev = []
    mean_ds = []
    for rna_class in rna_classes:
        print(rna_class)
        out = wilcoxon(d[rna_class][TS]["nat"][MS]["values"], d[rna_class][TS]["nat_rev"][MS]["values"])
        print(out)
        pvals.append(out[1])
        mean_d = d[rna_class][TS]["nat"][MS]["mean"] - d[rna_class][TS]["nat_rev"][MS]["mean"]
        mean_nat.append(d[rna_class][TS]["nat"][MS]["mean"])
        mean_nat_rev.append(d[rna_class][TS]["nat_rev"][MS]["mean"])
        print(mean_d)
        mean_ds.append(mean_d)

    foldchanges = list(np.log2(np.divide(mean_nat, mean_nat_rev)))
    print(f"{foldchanges=}")

    transformed_pvals = list(-1*np.log10(1*np.array(pvals)))  # here, no Bonferroni correction
    print(f"{transformed_pvals=}")

    plot_title = f"Volcano plot: normal vs. reversed direction, {TS}s" #@param {type:"string"}
    x_axis_title = f"log2 fold change ({MS}: normal/reversed)" #@param {type:"string"}
    y_axis_title = "-log10 pvalue" #@param {type:"string"}
    point_radius = 10 #@param {type:"slider", min:1, max:30, step:1}

    fig = go.Figure()

    fig.update_layout(
        title=plot_title,
        xaxis_title= x_axis_title,
        yaxis_title=y_axis_title,
        xaxis=dict(
            range=[-0.4, 0.4]
        ),
        yaxis=dict(
            range=[0, 25]
        )
        # paper_bgcolor= 'white',
        # plot_bgcolor='white',
    )

    colors = []

    for i in range(0, len(foldchanges)):

        if transformed_pvals[i] > 1.301029996:  # significance level: 0.05

            if foldchanges[i] > 0.0:
                colors.append('#db3232')
            elif foldchanges[i] < -0.0:
                colors.append('#3f65d4')
            else:
                colors.append('rgba(150,150,150,0.5)')
        else:
            colors.append('rgba(150,150,150,0.5)')

    # horizontal line representing significance level
    fig.add_shape(type='line',
                x0=-0.4,
                y0=1.301029996,
                x1=0.4,
                y1=1.301029996,
                line=dict(color='black', dash="dot"),
                xref='x',
                yref='y'
    )

    fig.add_trace(
        go.Scatter(
            x = foldchanges,
            y = transformed_pvals,
            mode = 'markers+text',
            text = rna_classes,
            textposition="top center",
            marker= {
                'color':colors,
                'size':point_radius,
            }
        )
    )

    print("Saving volcano...")
    fig.write_image(f'../{directory}/plots/volcano/volcano_{MS}_{TS}s_3.png', engine='kaleido')


def main():
    """Statistics on json-files.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Analysis of a database in dot-bracket format.')

    parser.add_argument("--input", metavar = '<str>',
            help = """Path of the input file""")
    parser.add_argument("--rna_classes", nargs='+', default = ["5s", "16s", "23s", "grpI", "grpII", "RNaseP", "srp", "telomerase", "tmRNA", "tRNA"], metavar = '<str>',
            help = """Space-separated list of RNA classes""")
    parser.add_argument("--ts", metavar = '<str>', default="1.0",
            help = """Time after transcription""")
    parser.add_argument("--directory",
            help = """Name of directory to analyze (e.g. archiveII_stef_o-prune_0.1).""")

    args = parser.parse_args()

    for ms in ["MCC", "bpdis"]:
        # plots(directory=args.directory, rna_classes=args.rna_classes, TS=args.ts, MS=ms)

        # combined_plots(rna_classes=args.rna_classes, TS=args.ts)
        for ts in ["0.04", "1.0", "10.0", "60.0", "600.0", "3600.0", "equilibrium"]:
            combined_plots(directory=args.directory, rna_classes=args.rna_classes, TS=ts, MS=ms)
 
if __name__ == '__main__':
    main()
