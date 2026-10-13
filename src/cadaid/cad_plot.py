"""Draw the directed graph with matplotlib."""
import numpy as np
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D


def draw_graph(ax, res):
    ax.clear()
    N = res.N
    ang = np.pi / 2 - 2 * np.pi * np.arange(N) / N
    P = np.c_[np.cos(ang), np.sin(ang)]

    groups = {}
    for k, (s, d) in enumerate(zip(res.src, res.dst)):
        groups.setdefault((min(s, d), max(s, d)), []).append(k)

    tree = set(res.tree)
    for ks in groups.values():
        m = len(ks)
        for i, k in enumerate(ks):
            rad = (i - (m - 1) / 2) * 0.4
            if res.src[k] > res.dst[k]:
                rad = -rad
            p0, p1 = P[res.src[k]], P[res.dst[k]]
            is_tree = k in tree
            ax.add_patch(FancyArrowPatch(
                p0, p1, connectionstyle=f"arc3,rad={rad}",
                arrowstyle="-|>", mutation_scale=14,
                shrinkA=12, shrinkB=12, lw=1.8,
                color="tab:blue" if is_tree else "tab:red",
                linestyle="-" if is_tree else "--"))
            d = p1 - p0
            mid = (p0 + p1) / 2 + 0.5 * rad * np.array([d[1], -d[0]])
            ax.text(mid[0], mid[1], res.names[k], ha="center", va="center",
                    fontsize=11, zorder=4,
                    bbox=dict(fc="white", ec="none", pad=1.0))

    ax.scatter(P[:, 0], P[:, 1], s=500, c="white", edgecolors="black", zorder=3)
    for i, (x, y) in enumerate(P):
        ax.text(x, y, str(i + 1), ha="center", va="center", zorder=5, fontsize=11)

    ax.legend(handles=[
        Line2D([0], [0], color="tab:blue", lw=2, label="tree branch"),
        Line2D([0], [0], color="tab:red", lw=2, ls="--", label="link"),
    ], loc="upper right", fontsize=9)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.axis("off")
