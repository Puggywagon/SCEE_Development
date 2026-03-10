#!/usr/bin/python3

import argparse
import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import glob


################################################################################
################################################################################
################################################################################
def read_data(T, p):
    csvfile_list = glob.glob(f'../../Replica_Dipoles/Replica_*/{T}K/{p:.1f}Bar/hbanalysis.csv')
    df = pd.DataFrame()
    for csvfile in csvfile_list:
        tmp = pd.read_csv(csvfile)
        df = pd.concat([df, tmp])
    return df
################################################################################
################################################################################
################################################################################
parser = argparse.ArgumentParser()

parser.add_argument('--show', type=str, default='True', help='plot figures')

args = parser.parse_args()

show = (args.show == 'True')


fontsize = 10
ticksize = 10
legendsize = 6
################################################################################
################################################################################
################################################################################
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]

donors_min, donors_max = 0, 4           # adjust if needed
mu_max, mu_min = 1.65, 3.75       # reuse your limits

donors_bins = np.arange(donors_min - 0.5, donors_max + 1.5, 1.0)   # centred on integers
mu_bins = np.linspace(mu_max, mu_min, 20)                 # smoother in 

cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(p_list))]
################################################################################
T = 298
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
T = 400
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
T = 500
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
T = 600
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
T = 700
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
T = 800
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
T = 900
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
T = 1000
fig, ax = plt.subplots(figsize=(3.5, 3))

print(f"temperature = {T} K")

for p, colour in zip(p_list, color_list):
    df = read_data(T, p)

    # basic guard
    if df.empty or ("muL" not in df.columns) or ("donors" not in df.columns):
        continue

    # 2D histogram density: H has shape (len(donors_bins)-1, len(mu_bins)-1)
    H, xedges, yedges = np.histogram2d(
        df["muL"].to_numpy(),
        df["donors"].to_numpy(),
        bins=[donors_bins, mu_bins],
        density=True
    )

    # grid of bin centres for contouring
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent, indexing="ij")

    # choose a couple of contour levels relative to the peak for each pressure
    # this makes shapes comparable even if absolute density changes
    Hmax = np.nanmax(H)
    if not np.isfinite(Hmax) or Hmax == 0:
        continue

    levels = [0.2 * Hmax, 0.5 * Hmax]   # tweak if you want tighter or broader contours
    ax.contour(X, Y, H, levels=levels, colors=[colour], linewidths=1.0)

    # optional: label in legend as proxy line
    ax.plot([], [], color=colour, label=f"p={p:.0f} bar")

ax.set_ylim([donors_min - 0.5, donors_max + 0.5])
ax.set_xlim([mu_max, mu_min])

ax.set_ylabel(r"$n_{\mathrm{HB}}$", fontsize=fontsize)
ax.set_xlabel(r"$\mu$ / D", fontsize=fontsize)

ax.set_yticks(range(donors_min, donors_max + 1))
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)

ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", prop={"size": legendsize})
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig(f"../dipdist_donors_{T}.pdf")
fig.savefig(f"../dipdist_donors_{T}.png")

if show:
    plt.show()
else:
    plt.close(fig)
################################################################################
