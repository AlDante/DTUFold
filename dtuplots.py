import sys
from string import ascii_uppercase, ascii_lowercase


import matplotlib
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np
import py3Dmol
from matplotlib import collections as mcoll
from matplotlib import colors
from matplotlib import cm

import dtucolabfold_utils

pymol_color_list = ["#33ff33", "#00ffff", "#ff33cc", "#ffff00", "#ff9999", "#e5e5e5", "#7f7fff", "#ff7f00",
                    "#7fff7f", "#199999", "#ff007f", "#ffdd5e", "#8c3f99", "#b2b2b2", "#007fff", "#c4b200",
                    "#8cb266", "#00bfbf", "#b27f7f", "#fcd1a5", "#ff7f7f", "#ffbfdd", "#7fffff", "#ffff7f",
                    "#00ff7f", "#337fcc", "#d8337f", "#bfff3f", "#ff7fff", "#d8d8ff", "#3fffbf", "#b78c4c",
                    "#339933", "#66b2b2", "#ba8c84", "#84bf00", "#b24c66", "#7f7f7f", "#3f3fa5", "#a5512b"]

pymol_cmap = matplotlib.colors.ListedColormap(pymol_color_list)
alphabet_list = list(ascii_uppercase + ascii_lowercase)


# @title Procedure make_plots()
def make_plots(msa, jobname, query_sequence, outs, homooligomer, num_models):
    # gather MSA info
    deduped_full_msa = list(dict.fromkeys(msa))
    msa_arr = np.array([list(seq) for seq in deduped_full_msa])
    seqid = (np.array(list(query_sequence)) == msa_arr).mean(-1)
    seqid_sort = seqid.argsort()  # [::-1]
    non_gaps = (msa_arr != "-").astype(float)
    non_gaps[non_gaps == 0] = np.nan
    ##################################################################
    plt.figure(figsize=(14, 4), dpi=100)
    ##################################################################
    plt.subplot(1, 2, 1)
    plt.title("Sequence coverage")
    plt.imshow(non_gaps[seqid_sort] * seqid[seqid_sort, None],
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa_arr != "-").sum(0), color='black')
    plt.xlim(-0.5, msa_arr.shape[1] - 0.5)
    plt.ylim(-0.5, msa_arr.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")

    ##################################################################
    plt.subplot(1, 2, 2)
    plt.title("Predicted lDDT per position")
    for model_name, value in outs.items():
        plt.plot(value["plddt"], label=model_name)
    if homooligomer > 0:
        for n in range(homooligomer + 1):
            x = n * (len(query_sequence) - 1)
            plt.plot([x, x], [0, 100], color="black")
    plt.legend()
    plt.ylim(0, 100)
    plt.ylabel("Predicted lDDT")
    plt.xlabel("Positions")
    plt.savefig(jobname + "_coverage_lDDT.png")
    ##################################################################
    # DEJ: plt.show()

    print("Predicted Alignment Error")
    ##################################################################
    plt.figure(figsize=(3 * num_models, 2), dpi=100)
    for n, (model_name, value) in enumerate(outs.items()):
        plt.subplot(1, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
    plt.savefig(jobname + "_PAE.png")
    # DEJ: plt.show()


# @title Function plot_plddt_legend()

def plot_plddt_legend():
    thresh = ['plDDT:', 'Very low (<50)', 'Low (60)', 'OK (70)', 'Confident (80)', 'Very high (>90)']
    plt.figure(figsize=(1, 0.1), dpi=100)
    ########################################
    for c in ["#FFFFFF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF"]:
        plt.bar(0, 0, color=c)
    plt.legend(thresh, frameon=False,
               loc='center', ncol=6,
               handletextpad=1,
               columnspacing=1,
               markerscale=0.5, )
    plt.axis(False)
    return plt


# @title Function plot_confidence()
def plot_confidence(Ls, plddt, pae, model_num=1, homooligomer=1):
    """

    :param Ls: Length of query sequence
    :param plddt:
    :param pae:
    :param model_num:
    :param homooligomer:
    :return:
    """
    model_name = f"model_{model_num}"
    plt.figure(figsize=(10, 3), dpi=100)
    """Plots the legend for plDDT."""
    #########################################
    plt.subplot(1, 2, 1)
    plt.title('Predicted lDDT')
    # plt.plot(outs[model_name]["plddt"])
    plt.plot(plddt)
    for n in range(homooligomer + 1):
        # x = n * (len(query_sequence))
        x = n * Ls
        plt.plot([x, x], [0, 100], color="black")
    plt.ylabel('plDDT')
    plt.xlabel('position')
    #########################################
    plt.subplot(1, 2, 2)
    plt.title('Predicted Aligned Error')
    # plt.imshow(outs[model_name]["pae"], cmap="bwr", vmin=0, vmax=30)
    plt.imshow(pae, cmap="bwr", vmin=0, vmax=30)
    plt.colorbar()
    plt.xlabel('Scored residue')
    plt.ylabel('Aligned residue')
    #########################################
    return plt


# @title Function show_pdb()
def show_pdb(jobname: str, model_num: int = 1, homooligomer: int = 1, use_amber: bool = True,
             show_mainchains: bool = False,
             show_sidechains: bool = False, color="lDDT"):
    """

    :param jobname:         Filename of output job
    :param model_num:       Model number being processed (1-5)
    :param homooligomer:    Number of homooligomers. 1 for no homooligomers
    :param use_amber:       True if Amber relaxation to be used
    :param show_mainchains: True to display main chains in plot
    :param show_sidechains: True to display sidechains in plot
    :param color:           Colour scheme. LDDT, rainbow or chain.
    :return:
    """
    model_name = f"model_{model_num}"
    if use_amber:
        pdb_filename = f"{jobname}_relaxed_{model_name}.pdb"
    else:
        pdb_filename = f"{jobname}_unrelaxed_{model_name}.pdb"

    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js', )
    view.addModel(open(pdb_filename, 'r').read(), 'pdb')

    if color == "lDDT":
        view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 50, 'max': 90}}})
    elif color == "rainbow":
        view.setStyle({'cartoon': {'color': 'spectrum'}})
    elif color == "chain":
        for n, chain, color in zip(range(homooligomer), list("ABCDEFGH"),
                                   ["lime", "cyan", "magenta", "yellow", "salmon", "white", "blue", "orange"]):
            view.setStyle({'chain': chain}, {'cartoon': {'color': color}})
    if show_sidechains:
        BB = ['C', 'O', 'N']
        view.addStyle({'and': [{'resn': ["GLY", "PRO"], 'invert': True}, {'atom': BB, 'invert': True}]},
                      {'stick': {'colorscheme': f"WhiteCarbon", 'radius': 0.3}})
        view.addStyle({'and': [{'resn': "GLY"}, {'atom': 'CA'}]},
                      {'sphere': {'colorscheme': f"WhiteCarbon", 'radius': 0.3}})
        view.addStyle({'and': [{'resn': "PRO"}, {'atom': ['C', 'O'], 'invert': True}]},
                      {'stick': {'colorscheme': f"WhiteCarbon", 'radius': 0.3}})
    if show_mainchains:
        BB = ['C', 'O', 'N', 'CA']
        view.addStyle({'atom': BB}, {'stick': {'colorscheme': f"WhiteCarbon", 'radius': 0.3}})

    view.zoomTo()
    return view


def plot_ticks(Ls):
    Ln = sum(Ls)
    L_prev = 0
    for L_i in Ls[:-1]:
        L = L_prev + L_i
        L_prev += L_i
        plt.plot([0, Ln], [L, L], color="black")
        plt.plot([L, L], [0, Ln], color="black")
    ticks = np.cumsum([0] + Ls)
    ticks = (ticks[1:] + ticks[:-1]) / 2
    plt.yticks(ticks, alphabet_list[:len(ticks)])


def plot_paes(paes, Ls=None, dpi=100, fig=True):
    """

    :param paes:
    :param Ls:
    :param dpi:
    :param fig:
    :return:
    """
    num_models = len(paes)
    if fig: plt.figure(figsize=(3 * num_models, 2), dpi=dpi)
    for n, pae in enumerate(paes):
        plt.subplot(1, num_models, n + 1)
        plt.title(f"rank_{n + 1}")
        Ln = pae.shape[0]
        plt.imshow(pae, cmap="bwr", vmin=0, vmax=30, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1: plot_ticks(Ls)
        plt.colorbar()
    return plt


def plot_adjs(adjs, Ls=None, dpi=100, fig=True):
    num_models = len(adjs)
    if fig: plt.figure(figsize=(3 * num_models, 2), dpi=dpi)
    for n, adj in enumerate(adjs):
        plt.subplot(1, num_models, n + 1)
        plt.title(f"rank_{n + 1}")
        Ln = adj.shape[0]
        plt.imshow(adj, cmap="binary", vmin=0, vmax=1, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1: plot_ticks(Ls)
        plt.colorbar()
    return plt


def plot_dists(dists, Ls=None, dpi=100, fig=True):
    num_models = len(dists)
    if fig: plt.figure(figsize=(3 * num_models, 2), dpi=dpi)
    for n, dist in enumerate(dists):
        plt.subplot(1, num_models, n + 1)
        plt.title(f"rank_{n + 1}")
        Ln = dist.shape[0]
        plt.imshow(dist, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1: plot_ticks(Ls)
        plt.colorbar()
    return plt


def plot_pseudo_3D(xyz, c=None, ax=None, chainbreak=5,
                   cmapname:str="gist_rainbow", pymol_cmap:matplotlib.colors.ListedColormap=None, line_w=2.0,
                   cmin=None, cmax=None, zmin=None, zmax=None):
    def rescale(a, amin=None, amax=None):
        a = np.copy(a)
        if amin is None: amin = a.min(initial=sys.float_info.max)
        if amax is None: amax = a.max(initial=sys.float_info.min)
        a[a < amin] = amin
        a[a > amax] = amax
        return (a - amin) / (amax - amin)

    # make segments
    xyz = np.asarray(xyz)
    seg = np.concatenate([xyz[:-1, None, :], xyz[1:, None, :]], axis=-2)
    seg_xy = seg[..., :2]
    seg_z = seg[..., 2].mean(-1)
    zord = seg_z.argsort()

    # set colors
    if c is None:
        c = np.arange(len(seg))[::-1]
    else:
        c = (c[1:] + c[:-1]) / 2
    c = rescale(c, cmin, cmax)

    validated_cmap:matplotlib.colors.Colormap

    if pymol_cmap is not None:
        colors = pymol_cmap(c)
    else:
        if cmapname == "gist_rainbow":
            c *= 0.75
        colors = matplotlib.cm.get_cmap(cmapname)

    if chainbreak is not None:
        dist = np.linalg.norm(xyz[:-1] - xyz[1:], axis=-1)
        colors[..., 3] = (dist < chainbreak).astype(np.float)

    # add shade/tint based on z-dimension
    z = rescale(seg_z, zmin, zmax)[:, None]
    tint, shade = z / 3, (z + 2) / 3
    colors[:, :3] = colors[:, :3] + (1 - colors[:, :3]) * tint
    colors[:, :3] = colors[:, :3] * shade

    set_lim = False
    if ax is None:
        fig, ax = plt.subplots()
        fig.set_figwidth(5)
        fig.set_figheight(5)
        set_lim = True
    else:
        fig = ax.get_figure()
        if ax.get_xlim() == (0, 1):
            set_lim = True

    if set_lim:
        xy_min = xyz[:, :2].min() - line_w
        xy_max = xyz[:, :2].max() + line_w
        ax.set_xlim(xy_min, xy_max)
        ax.set_ylim(xy_min, xy_max)

    ax.set_aspect('equal')

    # determine linewidths
    width = fig.bbox_inches.width * ax.get_position().width
    linewidths = line_w * 72 * width / np.diff(ax.get_xlim())

    lines = mcoll.LineCollection(seg_xy[zord], colors=colors[zord], linewidths=linewidths,
                                 path_effects=[matplotlib.patheffects.Stroke(capstyle="round")])

    return ax.add_collection(lines)


def add_text(text, ax):
    return plt.text(0.5, 1.01, text, horizontalalignment='center',
                    verticalalignment='bottom', transform=ax.transAxes)


def plot_protein(protein=None, pos=None, plddt=None, Ls=None, dpi=100, best_view=True, line_w=2.0):
    if protein is not None:
        pos = np.asarray(protein.atom_positions[:, 1, :])
        plddt = np.asarray(protein.b_factors[:, 0])

    # get best view
    if best_view:
        if plddt is not None:
            weights = plddt / 100
            pos = pos - (pos * weights[:, None]).sum(0, keepdims=True) / weights.sum()
            pos = pos @ dtucolabfold_utils.kabsch(pos, pos, weights, return_v=True)
        else:
            pos = pos - pos.mean(0, keepdims=True)
            pos = pos @ dtucolabfold_utils.kabsch(pos, pos, return_v=True)

    if plddt is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_figwidth(6)
        fig.set_figheight(3)
        ax = [ax1, ax2]
    else:
        fig, ax1 = plt.subplots(1, 1)
        fig.set_figwidth(3)
        fig.set_figheight(3)
        ax = [ax1]

    fig.set_dpi(dpi)
    fig.subplots_adjust(top=0.9, bottom=0.1, right=1, left=0, hspace=0, wspace=0)

    xy_min = pos[..., :2].min() - line_w
    xy_max = pos[..., :2].max() + line_w
    for a in ax:
        a.set_xlim(xy_min, xy_max)
        a.set_ylim(xy_min, xy_max)
        a.axis(False)

    if Ls is None or len(Ls) == 1:
        # color N->C
        c = np.arange(len(pos))[::-1]
        plot_pseudo_3D(pos, line_w=line_w, ax=ax1)
        add_text("colored by N→C", ax1)
    else:
        # color by chain
        c = np.concatenate([[n] * L for n, L in enumerate(Ls)])
        if len(Ls) > 40:
            plot_pseudo_3D(pos, c=c, line_w=line_w, ax=ax1)
        else:
            plot_pseudo_3D(pos, c=c, pymol_cmap=pymol_cmap, cmin=0, cmax=39, line_w=line_w, ax=ax1)
        add_text("colored by chain", ax1)

    if plddt is not None:
        # color by pLDDT
        plot_pseudo_3D(pos, c=plddt, cmin=50, cmax=90, line_w=line_w, ax=ax2)
        add_text("colored by pLDDT", ax2)

    return fig
