import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def plot_heatmap_real(matrice, name):
    plt.imshow(matrice, cmap='viridis', origin='upper')
    plt.colorbar(label='Valore')
    plt.title(name)
    plt.xlim(0.5, matrice.shape[1] - 0.5)
    plt.ylim(matrice.shape[0] - 0.5, 0.5)  # inverti asse y per mantenere "origin='upper'"
    plt.xlabel('Colonna')
    plt.ylabel('Riga')
    plt.show()

def plot_heatmap_cplx(matrice, nome):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    A = np.real(matrice)
    B = np.imag(matrice)
    im1 = axes[0].imshow(A, cmap='viridis', origin='upper')
    axes[0].set_title("Real")
    plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)

    # Seconda heatmap
    im2 = axes[1].imshow(B, cmap='plasma', origin='upper')
    axes[1].set_title("Imaginary")
    plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
    fig.suptitle(nome)
    # Layout ordinato
    plt.tight_layout()
    plt.show()
def plot_geom(coords, charges):
    """
    Visualizza sfere 3D colorate in base alla carica e collegate in sequenza.

    Parameters
    ----------
    coords : array-like, shape (N, 3)
        Coordinate cartesiane dei punti.
    charges : array-like, shape (N,)
        Valori di carica (usati per il colore).
    """
    coords = np.asarray(coords)
    charges = np.asarray(charges)
    n = len(charges)

    # Colori in base alla carica
    cmap = plt.cm.coolwarm
    norm = plt.Normalize(vmin=np.min(charges), vmax=np.max(charges))
    colors = cmap(norm(charges))

    # Figura 3D
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Sfere (scatter 3D)
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
               s=800, c=colors, edgecolor='k', alpha=0.9)

    # Etichette (numero d’ordine dentro la sfera)
    for i, (x, y, z) in enumerate(coords):
        ax.text(x, y, z, str(i + 1), color='black',
                ha='center', va='center', fontsize=10, weight='bold')

    # Collega solo la sfera i con i+1
    for i in range(n - 1):
        x_line = [coords[i, 0], coords[i + 1, 0]]
        y_line = [coords[i, 1], coords[i + 1, 1]]
        z_line = [coords[i, 2], coords[i + 1, 2]]
        ax.plot(x_line, y_line, z_line, color='gray', linewidth=2)

    # Impostazioni grafiche
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1, 1, 1])  # rapporti uguali
    plt.tight_layout()
    plt.show()

def plot_curve(x, y, name, xlabel, ylabel):
    """
    Plotta una curva 2D con etichette e titolo.

    Parameters
    ----------
    x : array-like
        Dati per l'asse x.
    y : array-like
        Dati per l'asse y.
    name : str
        Titolo del grafico.
    xlabel : str
        Etichetta dell'asse x.
    ylabel : str
        Etichetta dell'asse y.
    """
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, linestyle='-', color='b', linewidth=3)
    plt.title(name)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.tight_layout()
    plt.xlim(np.min(x), np.max(x))
    plt.show()
