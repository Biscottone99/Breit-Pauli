#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys


def main():
    angle = 75
    files = [
        f"u_9/{angle}.dat",
        f"u_9.5/{angle}.dat",
        f"u_9.75/{angle}.dat",
        f"u_10/{angle}.dat",
        f"u_11/{angle}.dat",
    ]

    name = "Non Coherent Hopping, 1111 system"
    labels = [
        "U = 9 eV",
        "U = 9.5 eV",
        "U = 9.75 eV",
        "U = 10 eV",
        "U = 11 eV",
    ]

    styles = [
        {"color": "paleturquoise", "lw": 2.8, "ls": "-"},
        {"color": "coral", "lw": 2.8, "ls": "-"},
        {"color": "forestgreen", "lw": 2.8, "ls": "-"},
        {"color": "navy", "lw": 2.8, "ls": "-"},
        {"color": "slateblue", "lw": 2.8, "ls": "-"},
    ]

    plt.rcParams.update({
        "figure.figsize": (11, 7),
        "axes.labelsize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 16,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.alpha": 0.25,
    })

    fig, ax = plt.subplots()

    max_x_global = 0
    not_found = []

    for fpath, lab, sty in zip(files, labels, styles):
        f = Path(fpath)
        if not f.exists():
            not_found.append(f.name)
            continue

        try:
            data = np.loadtxt(f, comments="#", usecols=(0, 1))
            data = np.atleast_2d(data)

            x = data[:, 0] * 1e-3  # converti ps → ns
            y = data[:, 1] * 50.0  # scaling polarizzazione

            max_x_global = max(max_x_global, x.max())

            ax.plot(x, y, label=lab, **sty)

        except Exception as e:
            print(f"❌ Errore leggendo '{f}': {e}", file=sys.stderr)

    # Etichette e titolo
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Spin polarization (%)")
    ax.set_title(name, fontsize=24, fontweight="bold", pad=15)

    # Limiti assi
    ax.set_xlim(0, 4000)

    # Legenda centrata sotto il grafico
    ax.legend(
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.18),
        frameon=False
    )
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.subplots_adjust(bottom=0.12)

    nome = f"{name.replace(',', '').replace(' ', '_')}.png"
    fig.savefig(nome, dpi=300, bbox_inches="tight")

    if not_found:
        print("⚠️ File non trovati:", ", ".join(not_found), file=sys.stderr)
    else:
        print(f"✅ Grafico creato e salvato come: {nome}")

    plt.show()


if __name__ == "__main__":
    main()
