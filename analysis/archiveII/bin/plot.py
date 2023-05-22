#!/usr/bin/env python3

import ternary
import random
import json

def random_points(num_points=25, scale=40):
    points = []
    for i in range(num_points):
        x = random.randint(1, scale)
        y = random.randint(0, scale - x)
        z = scale - x - y
        points.append((x,y,z))
    return points

def ternary_plot(db):
    points = []
    sums = []
    for header in db:
        x = db[header][f"dist_nat_tend"]
        y = db[header][f"dist_nat_tend_rev"]
        z = db[header][f"dist_nat_subopt"]
        sum = x + y + z
        x /= sum
        y /= sum
        z /= sum
        sums.append(sum)
        points.append((x*100,y*100,z*100))
    max_sum = max(sums)
    sums_normalized = [x/max_sum for x in sums]
    # Scatter Plot
    figure, tax = ternary.figure(scale=100)
    # figure.set_size_inches(10, 10) TODO: triangle not equilateral if heatmap is shown and eight=width
    tax.scatter(points, s=30, c=sums_normalized, label="RNA", cmap='plasma', alpha=0.6, linewidths=0)
    tax.legend()

    tax.set_title("Distances to native structure", fontsize=20)
    tax.boundary(linewidth=2.0)
    tax.gridlines(multiple=10, color="blue")
    tax.ticks(axis='lbr', linewidth=1, multiple=10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    offset=0.1
    fontsize=16
    tax.left_axis_label("dist_nat_subopt", fontsize=fontsize, offset=offset)
    tax.right_axis_label("dist_nat_tend_rev", fontsize=fontsize, offset=offset)
    tax.bottom_axis_label("dist_nat_tend", fontsize=fontsize, offset=offset)
    tax.heatmap(data={x:y for x,y in zip(points, sums)}, cmap='plasma', cbarlabel="Total distance")

    tax.show()

with open(f"/home/jotue/drtransformer/analysis/slurm_out/srp_short/srp_short.json") as f:
    db = json.load(f)

ternary_plot(db)