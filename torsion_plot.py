import numpy as np
import matplotlib.pyplot as plt
import sys

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

def main():
    name = 'Non Coherent Exponential Hopping, Hubbard Model, 1111 system'

    # Angoli disponibili (in gradi)
    angles = [15, 30, 45, 60, 75, 90]

    # 🔁 Inverti l’ordine (ora da 90° → 15°)
    angles = angles[::-1]

    # Genera automaticamente i file e le label
    files = [f"{a}.dat" for a in angles]
    labels = [f"{a}°" for a in angles]

    # Stili coerenti (abbinati ai 6 file)
    colors = ["paleturquoise", "darkslategray", "lightseagreen", "lightsalmon", "sandybrown", "forestgreen"]

    styles = [
        {"color": col, "lw": 2.8, "ls": "-", "marker": None}
        for col in colors
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

            x = data[:, 0] * 1e-3
            y = data[:, 1] * 50.0

            if x.max() > max_x_global:
                max_x_global = x.max()

            plot_kwargs = {"lw": 2.8}
            plot_kwargs.update({k: v for k, v in sty.items() if v is not None})

            ax.plot(x, y, label=lab, **plot_kwargs)

        except Exception as e:
            print(f"❌ Errore leggendo '{f}': {e}", file=sys.stderr)

    # Etichette e titolo
    ax.set_xlabel("time (fs)")
    ax.set_ylabel("spin polarization %")
    ax.set_title(name, fontsize=24, fontweight="bold")

    # Imposta intervallo x da 0 al massimo globale
    #ax.set_xlim(0, max_x_global)
    ax.set_xlim(0, 10000)
    # Legenda centrata sotto il grafico
    ax.legend(
        ncol=6,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.16),
        frameon=False
    )
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.subplots_adjust(bottom=0.12)
    nome = f"{name}.png"
    #fig.savefig(nome, dpi=300, bbox_inches="tight")

    if not_found:
        print("⚠️ File non trovati:", ", ".join(not_found), file=sys.stderr)
    else:
        print("✅ Grafico creato e salvato.")

    plt.show()


if __name__ == "__main__":
    main()
