import numpy as np
import matplotlib.pyplot as plt

def main():
    # File da leggere
    files = [
        "1111_mono.dat",
        "1111_tot.dat",
        "2110_mono.dat",
        "2110_tot.dat"
    ]

    # Titoli dei singoli riquadri
    titles = [
        "One-electron term",
        "Two-electron term",
        "One-electron term",
        "Two-electron term"
    ]

    # Colori: arancione per 1 e 3, blu per 2 e 4
    colors = ['forestgreen', 'slateblue', 'forestgreen', 'slateblue']

    # Crea figura 2x2 (senza sharey)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    axes = axes.flatten()

    for i, (filename, ax) in enumerate(zip(files, axes)):
        try:
            data = np.loadtxt(filename)
        except Exception as e:
            print(f"Errore nel leggere {filename}: {e}")
            continue

        x = data[:, 0]*1e-3
        y = data[:, 1] * 50  # scala come nel tuo script originale

        ax.plot(x, y, linewidth=3.0, color=colors[i])
        ax.set_title(titles[i], fontsize=12, fontweight='bold')
        ax.grid(True)
        ax.set_xlim(0, 10000)

        # Linea a y=0
        ax.axhline(0, color='black', lw=1)

        # Label solo sui bordi esterni
        if i >= 2:
            ax.set_xlabel("Time (ns)", fontsize=12)
        if i % 2 == 0:
            ax.set_ylabel("Spin polarization (%)", fontsize=12)

    # Titoli centrali per i due sistemi
    fig.text(0.5, 0.94, "System 1", ha='center', va='center', fontsize=14, fontweight='bold')
    fig.text(0.5, 0.48, "System 2", ha='center', va='center', fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()


if __name__ == "__main__":
    main()
