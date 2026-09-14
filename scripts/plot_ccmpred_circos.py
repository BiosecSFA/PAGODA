#!/usr/bin/env python3
"""
ccmpred_circos.py

Load a CCMpred/CCMpredPy Potts model from a msgpack ".raw/.braw" (optionally gzipped),
compute pairwise interaction strength via Frobenius norm, and draw a circos-style chord plot
where stronger interactions get thicker lines.

Typical CCMpredPy output is a MessagePack file like:  model.braw.gz
  - x_pair shape: (L, L, 21, 21)   (often includes gap as the 21st state)
  - x_single shape: (L, 20)

Usage examples:
  python plot_ccmpred_circos.py data/purE_hisI_library/hisI_Sean.raw figures/hisI_circos.png --fasta data/purE_hisI_library/hisI.fasta  --label-every 20 2>&1
  python plot_ccmpred_circos.py data/purE_hisI_library/purE.raw figures/purE_circos.png --fasta data/purE_hisI_library/purE.fasta  --label-every 20 2>&1

"""
import argparse
import gzip
import math
import os
from typing import Tuple, Optional
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Wedge


# ---------------------------
# I/O: read msgpack raw/braw
# ---------------------------
from pycameox.ilp import load_ccmpred


# ---------------------------
# Amino acid color scheme
# ---------------------------
def load_aa_colors(color_scheme_path: str = "data/my_palette_scheme.json"):
    """Load amino acid color scheme from JSON file."""
    try:
        with open(color_scheme_path, 'r') as f:
            scheme = json.load(f)
        return scheme['colors']
    except:
        # Default color scheme if file not found
        return {
            "A": "#8d8488", "C": "#6c0c0a", "D": "#eec12f", "E": "#feb560",
            "F": "#366e48", "G": "#68a6db", "H": "#caebc9", "I": "#626111",
            "K": "#fe8448", "L": "#67582d", "M": "#825900", "N": "#fbe5ce",
            "P": "#ffacc8", "Q": "#e88a1b", "R": "#d44d06", "S": "#b69b8d",
            "T": "#a17677", "V": "#75681f", "W": "#034d78", "Y": "#348d7a"
        }

AA_COLORS = load_aa_colors()


# ---------------------------
# FASTA parsing
# ---------------------------
def read_fasta(fasta_path: str) -> str:
    """Read the first sequence from a FASTA file and return it as a string."""
    sequence = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            sequence.append(line)
    return ''.join(sequence).upper()


# ---------------------------
# Strength computation
# ---------------------------
import torch

def frobenius_strength_matrix(
        potts: "Potts",
        min_sep: int = 0
) -> np.ndarray:
    """
    Compute strength S_ij = ||J_ij||_F for all i,j from a Potts model object.
    The Potts object contains attributes: A (alphabet size, 21), L (length),
    W (coupling matrix), and h (fields).

    Pairs with |i-j| < max(1, min_sep) are left as zero.
    Returns a symmetric (L, L) float64 matrix.
    """
    L = potts.L
    A = potts.A

    # Reshape from (L*A, L*A) to (L, L, A, A)
    x_pair = potts.W.weight.data.reshape(L, A, L, A).transpose(1, 2).detach().cpu().numpy()

    # Prepare output and compute upper triangle, symmetrize
    S = np.zeros((L, L), dtype=np.float64)
    step = max(1, int(min_sep))

    for i in range(L):
        j_start = i + step
        if j_start >= L:
            continue
        for j in range(j_start, L):
            Jij = x_pair[i, j, :, :]
            # Frobenius norm of the (A x A) coupling block
            Sij = float(np.linalg.norm(Jij, ord="fro"))
            S[i, j] = Sij
            S[j, i] = Sij

    return S


def top_pairs_from_strength(
    S: np.ndarray,
    top: int = 200,
    threshold: Optional[float] = None
):
    """
    Return list of (i, j, s) sorted by descending strength.
    If threshold is given, keep pairs with s >= threshold * max_strength.
    Else keep top N.
    """
    L = S.shape[0]
    iu = np.triu_indices(L, k=1)
    vals = S[iu]

    # Filter out zeros
    mask = vals > 0
    i_all = iu[0][mask]
    j_all = iu[1][mask]
    v_all = vals[mask]

    if v_all.size == 0:
        return []

    vmax = float(v_all.max())
    if threshold is not None:
        keep = v_all >= (threshold * vmax)
        i_all, j_all, v_all = i_all[keep], j_all[keep], v_all[keep]
    else:
        # top-N
        order = np.argsort(v_all)[::-1]
        order = order[: min(top, order.size)]
        i_all, j_all, v_all = i_all[order], j_all[order], v_all[order]

    # Sort descending
    order = np.argsort(v_all)[::-1]
    return [(int(i), int(j), float(v)) for i, j, v in zip(i_all[order], j_all[order], v_all[order])]


# ---------------------------
# Circos plotting
# ---------------------------

def _polar_to_cart(theta: float, r: float) -> Tuple[float, float]:
    return r * math.cos(theta), r * math.sin(theta)


def _chord_path(theta1: float, theta2: float, r: float, ctrl_r: float) -> Path:
    """
    Quadratic Bezier-ish chord using two cubic control points pulled inward.
    """
    x1, y1 = _polar_to_cart(theta1, r)
    x2, y2 = _polar_to_cart(theta2, r)

    # Control points: same angles, smaller radius
    cx1, cy1 = _polar_to_cart(theta1, ctrl_r)
    cx2, cy2 = _polar_to_cart(theta2, ctrl_r)

    verts = [
        (x1, y1),    # MOVETO
        (cx1, cy1),  # CURVE4 control 1
        (cx2, cy2),  # CURVE4 control 2
        (x2, y2),    # CURVE4 end
    ]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    return Path(verts, codes)


def plot_circos(
    L: int,
    pairs,
    outpath: str,
    sequence: Optional[str] = None,
    title: str = "",
    label_every: int = 0,
    start_angle_deg: float = 90.0,
    radius: float = 1.0,
    ctrl_radius: float = 0.15,
    min_lw: float = 1.0,
    max_lw: float = 8.0,
    alpha: float = 0.7,
    dpi: int = 300
):
    """
    Draw nodes around a circle and chords between interacting pairs.
    Line width scales with strength.
    Outer circle is colored by amino acid type.
    Chords are colored based on amino acids at both ends.
    """
    if not pairs:
        raise ValueError("No pairs to plot (after filtering). Try lowering --min-sep or --threshold, or increasing --top.")

    strengths = np.array([s for _, _, s in pairs], dtype=np.float64)
    smin, smax = float(strengths.min()), float(strengths.max())
    if smax == smin:
        norm = np.ones_like(strengths)
    else:
        norm = (strengths - smin) / (smax - smin)

    # Angles for each position
    start = math.radians(start_angle_deg)
    thetas = np.linspace(0, 2 * math.pi, L, endpoint=False) + start

    # Angular width for each position
    dtheta = 2 * math.pi / L

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)
    ax.set_aspect("equal")
    ax.axis("off")

    # Draw colored arcs for each amino acid position
    arc_width = 0.05  # Width of the colored arc ring
    default_color = "#888888"  # Gray for unknown amino acids

    if sequence:
        print(f"Drawing {L} positions with sequence colors (seq length: {len(sequence)})")
    else:
        print(f"Drawing {L} positions with default gray (no sequence provided)")

    for i in range(L):
        th_start = math.degrees(thetas[i] - dtheta/2)
        th_end = math.degrees(thetas[i] + dtheta/2)

        # Determine color based on sequence if available
        if sequence and i < len(sequence):
            aa = sequence[i].upper()
            color = AA_COLORS.get(aa, default_color)
        else:
            color = default_color

        # Draw arc as a wedge
        wedge = Wedge(
            center=(0, 0),
            r=radius,
            theta1=th_start,
            theta2=th_end,
            width=arc_width,
            facecolor=color,
            edgecolor='white',
            linewidth=0.5
        )
        ax.add_patch(wedge)

        # Optional labels
        if label_every and (i % label_every == 0):
            th = thetas[i]
            x_lab, y_lab = _polar_to_cart(th, radius * 1.05)
            ax.text(x_lab, y_lab, str(i + 1), ha="center", va="center", fontsize=12)

    # Chords - colored based on amino acids at both ends
    for (k, (i, j, s)) in enumerate(pairs):
        th1, th2 = float(thetas[i]), float(thetas[j])
        lw = min_lw + float(norm[k]) * (max_lw - min_lw)

        # Determine chord color (blend colors from both amino acids)
        if sequence and i < len(sequence) and j < len(sequence):
            aa1 = sequence[i].upper()
            aa2 = sequence[j].upper()
            color1 = AA_COLORS.get(aa1, default_color)
            color2 = AA_COLORS.get(aa2, default_color)
            # Use the first amino acid's color (could also blend)
            chord_color = color1
        else:
            chord_color = "#555555"

        path = _chord_path(th1, th2, radius - arc_width, ctrl_radius)
        patch = PathPatch(path, fill=False, linewidth=lw, alpha=alpha, edgecolor=chord_color)
        ax.add_patch(patch)

    if title:
        ax.set_title(title, fontsize=14)

    # Frame limits
    lim = radius * 1.2
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


# ---------------------------
# CLI
# ---------------------------

def main():
    p = argparse.ArgumentParser(description="Circos plot of Potts couplings from CCMpred raw/braw msgpack file.")
    p.add_argument("raw_file", help="CCMpred/CCMpredPy msgpack raw/braw file (optionally .gz)")
    p.add_argument("out", help="Output image file (png/pdf/svg)")
    p.add_argument("--fasta", type=str, default=None,
                   help="FASTA file with protein sequence to color amino acid positions (optional).")
    p.add_argument("--top", type=int, default=200, help="Number of strongest pairs to plot (ignored if --threshold is set).")
    p.add_argument("--threshold", type=float, default=None,
                   help="Keep pairs with strength >= threshold * max_strength (e.g., 0.85). Overrides --top.")
    p.add_argument("--include-gap", action="store_true",
                   help="Include gap state in Frobenius norm (default: exclude last state if q==21).")
    p.add_argument("--min-sep", type=int, default=0, help="Minimum sequence separation |i-j| to include (e.g., 5).")
    p.add_argument("--label-every", type=int, default=20, help="Label every N-th position around the circle (0 disables).")
    p.add_argument("--title", type=str, default="", help="Plot title.")
    p.add_argument("--start-angle-deg", type=float, default=90.0, help="Angle offset for position 1 (degrees).")
    p.add_argument("--alpha", type=float, default=0.7, help="Chord transparency.")
    p.add_argument("--min-lw", type=float, default=1.0, help="Minimum chord line width.")
    p.add_argument("--max-lw", type=float, default=8.0, help="Maximum chord line width.")
    p.add_argument("--ctrl-radius", type=float, default=0.5,
                   help="Control radius for chords (smaller pulls chords inward more).")
    args = p.parse_args()

    potts = load_ccmpred(args.raw_file)

    # Read sequence from FASTA if provided
    sequence = None
    if args.fasta:
        sequence = read_fasta(args.fasta)
        print(f"Loaded sequence from {args.fasta}: length={len(sequence)}, first 20 AA: {sequence[:20]}")
        if len(sequence) != potts.L:
            print(f"WARNING: Sequence length ({len(sequence)}) != model length ({potts.L})")

    S = frobenius_strength_matrix(potts, min_sep=args.min_sep)
    pairs = top_pairs_from_strength(S, top=args.top, threshold=args.threshold)

    if not args.title:
        base = os.path.basename(args.raw_file)
        title = f"Potts couplings (Frobenius) — {base}  |  plotted={len(pairs)}"
    else:
        title = args.title

    plot_circos(
        L=potts.L,
        pairs=pairs,
        outpath=args.out,
        sequence=sequence,
        title=title,
        label_every=args.label_every,
        start_angle_deg=args.start_angle_deg,
        alpha=args.alpha,
        min_lw=args.min_lw,
        max_lw=args.max_lw,
        ctrl_radius=args.ctrl_radius,
    )

    print(f"Wrote: {args.out} (L={potts.L}, pairs={len(pairs)})")


if __name__ == "__main__":
    main()


