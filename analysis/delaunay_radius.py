#!/usr/bin/env python3

import os
import re
import csv
import math
import itertools

import progressbar
import sys
import datetime

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

import numpy as np
import pandas as pd


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def closeFig(fig, file):
    if not file.startswith('output/'):
        file = 'output/' + file

    if not '.' in file:
        file = file + '.png'

    ensure_dir(file)

    fig.savefig(file)
    plt.close(fig)


def saveData(fig, file):
    lines = fig.gca().get_lines()

    X = lines[0].get_xdata()

    data = np.zeros(shape=(len(X), len(lines) + 1))

    data[:, 0] = X
    for i in range(len(lines)):
        data[:, i + 1] = lines[i].get_ydata()

    if not file.startswith('output/'):
        file = 'output/' + file

    if not '.' in file:
        file = file + '.csv'

    np.savetxt(file, data)


# setup pandas
pd.set_option('mode.chained_assignment', None)

# prepare regex strings
r_stats_file_name = re.compile(r'delaunay_(?P<d>\d)d_(?P<n>\d+)_(?P<k>\d+)_stats_(?P<rank>\d+)')
r_cells_file_name = re.compile(r'delaunay_(?P<d>\d)d_(?P<n>\d+)_(?P<k>\d+)_cells_(?P<rank>\d+)')

statsfile = 'delaunay_stats.pkl'
cellsfile = 'delaunay_cells.pkl'
header_stats = ['d',
                'n',
                'k',
                'chunks_per_dim',
                'chunk_size',
                'cells_per_chunk',
                'cells_per_dim',
                'cell_size',
                'cell_type',
                'max_radius',
                'chunk',
                'radius',
                'n_p_chunk',
                'n_p_tria',
                'n_simplices']

header_cells = ['d',
                'n',
                'k',
                'chunks_per_dim',
                'chunk_size',
                'cells_per_chunk',
                'cells_per_dim',
                'cell_size',
                'cell_type',
                'max_radius',
                'chunk',
                'n_p_cell']

dir = 'logs'

if os.path.exists(statsfile):
    gStats = pd.read_pickle(statsfile)
    statsfileTS = os.stat(statsfile).st_mtime

    print("loaded stats database from %s" % datetime.datetime.fromtimestamp(statsfileTS))
else:
    gStats = pd.DataFrame(columns=header_stats)
    statsfileTS = 0

    print("creating new stats database")

if os.path.exists(dir):
    files = [f for f in os.listdir(dir) if ('_stats_' in f) and (os.stat(os.path.join(dir, f)).st_mtime > statsfileTS)]
    bar = progressbar.ProgressBar()
    frames = []

    if len(files):
        print("%i new stats files found" % len(files))
        for file in bar(files):
            match = r_stats_file_name.match(file)
            if match:
                # new file format
                d = int(match.group('d'))
                n = int(match.group('n'))
                k = int(match.group('k'))
                rank = int(match.group('rank'))
                iteration = 0

            frames.append(pd.read_csv(os.path.join(dir, file), delim_whitespace=True, names=header_stats))

        gStats = pd.concat(frames, ignore_index=True)
        gStats.sort_values(['d', 'n', 'k'], inplace=True)
        gStats.to_pickle(statsfile)

if os.path.exists(cellsfile):
    gCells = pd.read_pickle(cellsfile)
    cellsfileTS = os.stat(cellsfile).st_mtime

    print("loaded cells database from %s" % datetime.datetime.fromtimestamp(cellsfileTS))
else:
    gCells = pd.DataFrame(columns=header_cells)
    cellsfileTS = 0

    print("creating new cells database")

if os.path.exists(dir):
    files = [f for f in os.listdir(dir) if ('_cells_' in f) and (os.stat(os.path.join(dir, f)).st_mtime > cellsfileTS)]
    bar = progressbar.ProgressBar()
    frames = []

    if len(files):
        print("%i new cells files found" % len(files))
        for file in bar(files):
            match = r_cells_file_name.match(file)
            if match:
                # new file format
                d = int(match.group('d'))
                n = int(match.group('n'))
                k = int(match.group('k'))
                rank = int(match.group('rank'))
                iteration = 0

            frames.append(pd.read_csv(os.path.join(dir, file), delim_whitespace=True, names=header_cells))

        gCells = pd.concat(frames, ignore_index=True)
        gCells.sort_values(['d', 'n', 'k'], inplace=True)
        gCells.to_pickle(cellsfile)

# done loading new data files

class CompareFunction:
    __slots__ = ['f', 'name']

    def __init__(self, name: str, f):
        self.f = f
        self.name = name


def make_lineplot(data: pd.DataFrame, file: str,
                  x_field: str, y_field: str, z_field: str,
                  x_label: str, y_label: str, z_label: str,
                  x_transfrom=None, y_transform=None, z_transform=None,
                  x_agg: str = None, y_agg: str = 'mean', z_agg: str = None,
                  x_tick: str = r'%i', y_tick: str = r'%i', z_tick: str = r'%i',
                  x_tick_transform=lambda x: x, y_tick_transform=lambda y: y, z_tick_transform=lambda z: z,
                  x_scale=None, y_scale=None,
                  descr: str = None,
                  cmps: list = None
                  ):
    fig, ax = plt.subplots()

    X = np.sort(np.unique(data[x_field]))
    Z = np.sort(np.unique(data[z_field]))

    colors = plt.cm.Set1(np.linspace(0, 1, len(Z)))
    markers = itertools.cycle(['x', 'D', 'o', '*', '+', '^'])
    linestyles = itertools.cycle(['--', '-.', ':'])

    for z, c, m in zip(Z, colors, markers):
        select = data.loc[data[z_field] == z]

        if y_transform != None:
            y_field = 'transform'
            select[y_field] = y_transform(select)

        grouped = getattr(select.groupby(x_field), y_agg)()
        ax.plot(grouped.index, grouped[y_field], color=c, marker=m, ls='-', label=z_tick % z_tick_transform(z))

    if cmps:
        for cmp, m, ls in zip(cmps, markers, linestyles):
            ax.plot(X, cmp.f(X), color='k', marker=m, ls=ls, label=cmp.name)

    ax.set_xlabel(x_label)
    if x_scale != None:
        ax.set_xscale(**x_scale)
    # ax.set_xticklabels([x_tick % x_tick_transform(x) for x in X[1::2]])

    ax.set_ylabel(y_label)
    if y_scale != None:
        ax.set_yscale(**y_scale)

    ax.legend(loc=0, title=r'%s - %s' % (z_label, y_agg))

    if descr != None:
        fig.text(0.12, 0.9, descr, va='bottom', ha='left')

    saveData(fig, file)
    closeFig(fig, file)


def make_boxplot(data: pd.DataFrame, file: str,
                 x_field: str, y_field: str, z_field: str,
                 x_label: str, y_label: str, z_label: str,
                 x_transfrom=None, y_transform=None, z_transform=None,
                 x_tick: str = r'%i', y_tick: str = r'%i', z_tick: str = r'%i',
                 x_tick_transform=lambda x: x, y_tick_transform=lambda y: y, z_tick_transform=lambda z: z,
                 x_scale=None, y_scale=None,
                 descr: str = None
                 ):
    fig, ax = plt.subplots()

    X = np.sort(np.unique(data[x_field]))
    Z = np.sort(np.unique(data[z_field]))

    colors = plt.cm.Set1(np.linspace(0, 1, len(Z)))

    width = 2 / (3 * len(X))
    spacing = width / 3
    offset = spacing

    legItems = []
    for z, c in zip(Z, colors):
        select = data.loc[data[z_field] == z]

        if y_transform != None:
            y_field = 'transform'
            select[y_field] = y_transform(select)

        grouped = select.groupby(x_field)
        Y = []
        for x, g in grouped:
            Y.append(g[y_field].values)

        ax.boxplot(Y,
                   positions=np.arange(offset, len(X) + offset, 1),
                   widths=width, whis=[5, 95],
                   boxprops={'color': c},
                   whiskerprops={'color': c},
                   capprops={'color': c},
                   flierprops={'markeredgecolor': tuple(c), 'marker': '.'},
                   medianprops={'color': c}
                   )
        legItems.append(mpatches.Patch(color=c, label=z_tick % z_tick_transform(z)))
        offset += width + spacing

    ax.set_xlabel(x_label)
    ax.set_xlim(0, len(X) + 1)
    ax.set_xticklabels([x_tick % x_tick_transform(x) for x in X[1::2]])
    ax.set_xticks(np.arange(2, len(X) + 1, 2))

    ax.set_ylabel(y_label)
    if y_scale != None:
        ax.set_yscale(**y_scale)

    ax.legend(handles=legItems, bbox_to_anchor=(1.1, 1.1), title=z_label)
    if descr != None:
        fig.text(0.12, 0.9, descr, va='bottom', ha='left')

    closeFig(fig, file)


D = np.sort(np.unique(gStats['d']))
N = np.sort(np.unique(gStats['n']))
K = np.sort(np.unique(gStats['k']))

for d in D:
    print("processing %i dimensions" % d)
    data = gStats.loc[gStats['d'] == d]

    #######################################################################################################################

    lRadius = lambda x: x['radius'] * x['cell_size']
    lTNorm = lambda x: x['n_p_tria'] / x['n_p_chunk']
    lTNorm2 = lambda x: (x['n_p_tria'] - x['n_p_chunk']) / (x['n_p_chunk'] ** (x['d'] - 1 / x['d']))
    lNCells = lambda x: (1.0 / x['cell_size']) ** (x['d'])
    lPpCell = lambda x: x['n_p_chunk'] / (1.0 / x['cell_size']) ** (x['d'])

    #######################################################################################################################

    cmp12 = CompareFunction(r'$\frac{c}{n^{\frac{1}{2}}}$', lambda x: (4 / x) ** (1 / 2))
    cmp13 = CompareFunction(r'$\frac{c}{n^{\frac{1}{3}}}$', lambda x: (5 / x) ** (1 / 3))

    cmps = {
        2: [cmp12],
        3: [cmp13]
    }

    #######################################################################################################################

    make_lineplot(data, 'output/%id/n_r_mean' % d,
                  x_field='n', y_field=None, z_field='k',
                  x_label=r'$n$', y_label=r'$r_\mathrm{halo}$', z_label=r'$k$',
                  y_transform=lRadius,
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 10},
                  descr='%iD' % d,
                  cmps=cmps[d]
                  )

    make_lineplot(data, 'output/%id/n_r_max' % d,
                  x_field='n', y_field=None, z_field='k',
                  x_label=r'$n$', y_label=r'$r_\mathrm{halo}$', z_label=r'$k$',
                  y_transform=lRadius,
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 10},
                  descr='%iD' % d,
                  cmps=cmps[d]
                  )

    make_boxplot(data, 'output/%id/n_r_box' % d,
                 x_field='n', y_field=None, z_field='k',
                 x_label=r'$n$', y_label=r'$r_\mathrm{halo}$', z_label=r'$k$',
                 y_transform=lRadius,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 10},
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_boxplot(data, 'output/%id/n_n_per_chunk_box' % d,
                 x_field='n', y_field='n_p_chunk', z_field='k',
                 x_label=r'$n$', y_label=r'$n_\mathrm{chunk}$', z_label=r'$k$',
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 2},
                 descr='%iD' % d
                 )

    make_boxplot(data, 'output/%id/n_cell_size_box' % d,
                 x_field='n', y_field='cell_size', z_field='k',
                 x_label=r'$n$', y_label=r'$r_\mathrm{cell}$', z_label=r'$k$',
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 10},
                 descr='%iD' % d
                 )

    make_boxplot(data, 'output/%id/n_cells_box' % d,
                 x_field='n', y_field=None, z_field='k',
                 x_label=r'$n$', y_label=r'$n_\mathrm{cells}$', z_label=r'$k$',
                 y_transform=lNCells,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 2},
                 descr='%iD' % d
                 )

    make_boxplot(data, 'output/%id/n_n_per_cell_box' % d,
                 x_field='n', y_field=None, z_field='k',
                 x_label=r'$n$', y_label=r'$n_\mathrm{per\ cell}$', z_label=r'$k$',
                 y_transform=lPpCell,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 2},
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/n_c_mean' % d,
                  x_field='n', y_field='radius', z_field='k',
                  x_label=r'$n$', y_label=r'$r_\mathrm{cells}$', z_label=r'$k$',
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/n_c_max' % d,
                  x_field='n', y_field='radius', z_field='k',
                  x_label=r'$n$', y_label=r'$r_\mathrm{cells}$', z_label=r'$k$',
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/n_c_box' % d,
                 x_field='n', y_field='radius', z_field='k',
                 x_label=r'$n$', y_label=r'$r_\mathrm{cells}$', z_label=r'$k$',
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/n_tNorm_mean' % d,
                  x_field='n', y_field=None, z_field='k',
                  x_label=r'$n$', y_label=r'$\frac{n_\mathrm{triangulated}}{\frac{n}{k}}$', z_label=r'$k$',
                  y_transform=lTNorm,
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/n_tNorm_max' % d,
                  x_field='n', y_field=None, z_field='k',
                  x_label=r'$n$', y_label=r'$\frac{n_\mathrm{triangulated}}{\frac{n}{k}}$', z_label=r'$k$',
                  y_transform=lTNorm,
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/n_tNorm_box' % d,
                 x_field='n', y_field=None, z_field='k',
                 x_label=r'$n$', y_label=r'$\frac{n_\mathrm{triangulated}}{\frac{n}{k}}$', z_label=r'$k$',
                 y_transform=lTNorm,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/n_tNorm2_mean' % d,
                  x_field='n', y_field=None, z_field='k',
                  x_label=r'$n$',
                  y_label=r'$\frac{n_\mathrm{triangulated} - n_\mathrm{chunk}}{n_\mathrm{chunk}^\frac{D-1}{D}}$',
                  z_label=r'$k$',
                  y_transform=lTNorm2,
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/n_tNorm2_max' % d,
                  x_field='n', y_field=None, z_field='k',
                  x_label=r'$n$',
                  y_label=r'$\frac{n_\mathrm{triangulated} - n_\mathrm{chunk}}{n_\mathrm{chunk}^\frac{D-1}{D}}$',
                  z_label=r'$k$',
                  y_transform=lTNorm2,
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/n_tNorm2_box' % d,
                 x_field='n', y_field=None, z_field='k',
                 x_label=r'$n$',
                 y_label=r'$\frac{n_\mathrm{triangulated} - n_\mathrm{chunk}}{n_\mathrm{chunk}^\frac{D-1}{D}}$',
                 z_label=r'$k$',
                 y_transform=lTNorm2,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/n_t_mean' % d,
                  x_field='n', y_field='n_p_tria', z_field='k',
                  x_label=r'$n$', y_label=r'$n_\mathrm{triangulated}$', z_label=r'$k$',
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/n_t_max' % d,
                  x_field='n', y_field='n_p_tria', z_field='k',
                  x_label=r'$n$', y_label=r'$n_\mathrm{triangulated}$', z_label=r'$k$',
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/n_t_box' % d,
                 x_field='n', y_field='n_p_tria', z_field='k',
                 x_label=r'$n$', y_label=r'$n_\mathrm{triangulated}$', z_label=r'$k$',
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 2},
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/k_r_mean' % d,
                  x_field='k', y_field=None, z_field='n',
                  x_label=r'$k$', y_label=r'$r_\mathrm{halo}$', z_label=r'$n$',
                  y_transform=lRadius,
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 10},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/k_r_max' % d,
                  x_field='k', y_field=None, z_field='n',
                  x_label=r'$k$', y_label=r'$r_\mathrm{halo}$', z_label=r'$n$',
                  y_transform=lRadius,
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 10},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/k_r_box' % d,
                 x_field='k', y_field=None, z_field='n',
                 x_label=r'$k$', y_label=r'$r_\mathrm{halo}$', z_label=r'$n$',
                 y_transform=lRadius,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 10},
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/k_c_mean' % d,
                  x_field='k', y_field='radius', z_field='n',
                  x_label=r'$k$', y_label=r'$r_\mathrm{cells}$', z_label=r'$n$',
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/k_c_max' % d,
                  x_field='k', y_field='radius', z_field='n',
                  x_label=r'$k$', y_label=r'$r_\mathrm{cells}$', z_label=r'$n$',
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/k_c_box' % d,
                 x_field='k', y_field='radius', z_field='n',
                 x_label=r'$k$', y_label=r'$r_\mathrm{cells}$', z_label=r'$n$',
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/k_tNorm_mean' % d,
                  x_field='k', y_field=None, z_field='n',
                  x_label=r'$k$', y_label=r'$\frac{n_\mathrm{triangulated}}{n_\mathrm{cell}}$', z_label=r'$n$',
                  y_transform=lTNorm,
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/k_tNorm_max' % d,
                  x_field='k', y_field=None, z_field='n',
                  x_label=r'$k$', y_label=r'$\frac{n_\mathrm{triangulated}}{n_\mathrm{cell}}$', z_label=r'$n$',
                  y_transform=lTNorm,
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/k_tNorm_box' % d,
                 x_field='k', y_field=None, z_field='n',
                 x_label=r'$k$', y_label=r'$\frac{n_\mathrm{triangulated}}{\frac{n}{k}}$', z_label=r'$n$',
                 y_transform=lTNorm,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/k_tNorm2_mean' % d,
                  x_field='k', y_field=None, z_field='n',
                  x_label=r'$k$',
                  y_label=r'$\frac{n_\mathrm{triangulated} - n_\mathrm{chunk}}{n_\mathrm{chunk}^\frac{D-1}{D}}$',
                  z_label=r'$n$',
                  y_transform=lTNorm2,
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/k_tNorm2_max' % d,
                  x_field='k', y_field=None, z_field='n',
                  x_label=r'$k$',
                  y_label=r'$\frac{n_\mathrm{triangulated} - n_\mathrm{chunk}}{n_\mathrm{chunk}^\frac{D-1}{D}}$',
                  z_label=r'$n$',
                  y_transform=lTNorm2,
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/k_tNorm2_box' % d,
                 x_field='k', y_field=None, z_field='n',
                 x_label=r'$k$',
                 y_label=r'$\frac{n_\mathrm{triangulated} - n_\mathrm{chunk}}{n_\mathrm{chunk}^\frac{D-1}{D}}$',
                 z_label=r'$n$',
                 y_transform=lTNorm2,
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 descr='%iD' % d
                 )

    #######################################################################################################################

    make_lineplot(data, 'output/%id/k_t_mean' % d,
                  x_field='k', y_field='n_p_tria', z_field='n',
                  x_label=r'$k$', y_label=r'$n_\mathrm{triangulated}$', z_label=r'$n$',
                  y_agg='mean',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 2},
                  descr='%iD' % d
                  )

    make_lineplot(data, 'output/%id/k_t_max' % d,
                  x_field='k', y_field='n_p_tria', z_field='n',
                  x_label=r'$k$', y_label=r'$n_\mathrm{triangulated}$', z_label=r'$n$',
                  y_agg='max',
                  x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                  x_tick_transform=lambda x: math.log2(x),
                  z_tick_transform=lambda x: math.log2(x),
                  x_scale={'value': 'log', 'basex': 2},
                  y_scale={'value': 'log', 'basey': 2},
                  descr='%iD' % d
                  )

    make_boxplot(data, 'output/%id/k_t_box' % d,
                 x_field='k', y_field='n_p_tria', z_field='n',
                 x_label=r'$k$', y_label=r'$n_\mathrm{triangulated}$', z_label=r'$n$',
                 x_tick=r'$2^{%i}$', z_tick=r'$2^{%i}$',
                 x_tick_transform=lambda x: math.log2(x),
                 z_tick_transform=lambda x: math.log2(x),
                 y_scale={'value': 'log', 'basey': 2},
                 descr='%iD' % d
                 )
