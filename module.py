import numpy as np
import numpy as np
import os

import plotting_module as plot

imag = 1j
me = 9.1093837015e-31
gs = 2.00231930436256
e = 1.602176634e-19
e0 = 8.8541878128e-12
cl = 299792458
pf = ((gs * e ** 2) / (8 * np.pi * e0 * me * cl ** 2)) * 10.0e10
pauli_matrix = np.zeros((2, 2, 3), dtype=complex)
pauli_matrix[0, 1, 0] = 1.0
pauli_matrix[1, 0, 0] = 1.0

pauli_matrix[0, 1, 1] = -imag
pauli_matrix[1, 0, 1] = imag

pauli_matrix[0, 0, 2] = 1.0
pauli_matrix[1, 1, 2] = -1.0


def generate_basis(electrons):
    """
    Generate basis states in bit representation, arranged in order of increasing spin.

    Args:
        electrons (int): Number of electrons.

    Returns:
        configs (np.array): configurazioni in forma intera
        spins (np.array): spin associati
        nf (int): numero di configurazioni generate
    """
    # ==== Parametri ====
    nso = electrons * 2  # numero di spinorbitali

    # ==== Limiti del ciclo ====
    n_max = sum(2 ** i for i in range(nso - electrons, nso))
    n_min = sum(2 ** i for i in range(electrons))

    nf = 0
    usefull = []

    # ==== Loop principale ====
    for n in range(n_min, n_max + 1):
        count = 0
        a = 0  # spin-up
        b = 0  # spin-down
        array = ['0'] * nso

        for i in range(nso):
            if (n >> i) & 1:
                array[i] = '1'
                count += 1
                if i % 2 == 0:
                    a += 1
                else:
                    b += 1

        spin = 0.5 * (a - b)

        if count == electrons:
            config = sum(2 ** i for i in range(nso) if array[i] == '1')
            usefull.append((config, spin, ''.join(array)))
            nf += 1

    print(f"Numero totale di configurazioni valide: {nf}")

    # ==== Ordinamento per spin crescente ====
    usefull_sorted = sorted(usefull, key=lambda x: x[1])

    # ==== Creazione cartella se manca ====
    folder = '../hubbard_soc'
    os.makedirs(folder, exist_ok=True)

    # ==== Stampa su file ====
    with open(os.path.join(folder, 'configurations.txt'), 'w') as f3:
        for config, spin, array in usefull_sorted:
            f3.write(f"{array}  spin={spin:+.1f}  config={config}\n")

    # ==== Array NumPy da restituire ====
    configs = np.array([u[0] for u in usefull_sorted])
    spins = np.array([u[1] for u in usefull_sorted])

    return configs, spins, nf


def btest(n, i):
    return ((n >> i) & 1) == 1
def ibset(n, i):
    return n | (1 << i)
def ibclr(n, i):
    return n & ~(1 << i)

def linear_search(arr, val):
    """
    Cerca un valore in una lista non ordinata.

    Args:
        arr (list): Lista in cui cercare.
        val: Valore da cercare.

    Returns:
        int: indice del valore se trovato, -1 altrimenti.
    """
    for i, v in enumerate(arr):
        if v == val:
            return i
    return -1
def tb_to_rs(dim, nso, basis, op_tb):
    op = np.zeros((dim, dim), dtype=complex)
    for j in range(dim):  # indice di colonna
        for iso in range(nso):  # creation operator index
            for jso in range(nso):  # annihilation operator index
                jstate = basis[j]
                if btest(jstate, jso):
                    istate = ibclr(jstate, jso)
                    if not btest(istate, iso):
                        istate = ibset(istate, iso)
                        i = linear_search(basis, istate)
                        if i != -1:
                            # Calcolo fase
                            if jso != iso:
                                step = 1 if iso > jso else -1
                                conta = 0
                                for a in range(jso + step, iso, step):
                                    if btest(istate, a):
                                        conta += 1
                                phase = 1 if conta % 2 == 0 else -1
                            else:
                                phase = 1
                            # Aggiornamento matrice
                            op[i, j] += phase * op_tb[iso, jso]
    return op


def generate_rotated_geometry(nsiti, length, pitch, angle_deg, output_file="geom_rotate.dat"):
    """
    Genera le coordinate di un'elica (a partire da lunghezza, pitch e angolo),
    calcola la matrice di rotazione di Rodrigues per allineare l'asse molecolare
    (primo-ultimo atomo) con l'asse Z del laboratorio, e salva le coordinate ruotate.

    Parametri:
    - nsiti: numero di siti atomici
    - length: distanza (l)
    - pitch: passo lungo z (c)
    - angle_deg: angolo di torsione in gradi (es. 36 per pi/5)

    Ritorna:
    - np.ndarray con le coordinate ruotate (nsiti x 3)
    """
    # 1. Conversione angolo in radianti e calcolo del raggio
    alpha = np.radians(angle_deg)

    if alpha == 0:
        raise ValueError("L'angolo di torsione (alpha) non può essere zero.")

    num = length ** 2 - pitch ** 2
    if num < 0:
        raise ValueError("Errore fisico: la lunghezza (l) deve essere maggiore o uguale al passo (c).")

    r = np.sqrt(num / (2 * (1 - np.cos(alpha))))

    # 2. Generazione punti dell'elica nativa
    points = np.zeros((nsiti, 3))
    for i in range(nsiti):
        # Indice Python parte da 0, quindi equivale a (i-1) in Fortran
        points[i, 0] = r * np.cos(np.pi / 2 - i * alpha)
        points[i, 1] = r * np.cos(i * alpha)
        points[i, 2] = pitch * i

    # 3. Funzioni di supporto annidate per l'algebra vettoriale
    def normalize(v):
        norm = np.linalg.norm(v)
        return v if norm == 0 else v / norm

    def rodrigues_rotation_matrix(v, target_axis):
        v = normalize(v)
        target_axis = normalize(target_axis)
        k = np.cross(v, target_axis)
        sin_theta = np.linalg.norm(k)
        cos_theta = np.dot(v, target_axis)

        # Se i vettori sono già allineati o anti-allineati
        if sin_theta == 0:
            return np.eye(3) if cos_theta > 0 else -np.eye(3)

        k = normalize(k)
        K = np.array([
            [0, -k[2], k[1]],
            [k[2], 0, -k[0]],
            [-k[1], k[0], 0]
        ])

        # Formula di Rodrigues
        R = np.eye(3) + sin_theta * K + (1 - cos_theta) * np.dot(K, K)
        return R

    # 4. Calcolo asse e matrice di rotazione
    # Vettore dal primo all'ultimo sito (in Python l'indice -1 indica l'ultimo elemento)
    vector_ad = points[-1] - points[0]
    z_axis = np.array([0.0, 0.0, 1.0])

    R = rodrigues_rotation_matrix(vector_ad, z_axis)

    # 5. Applica la rotazione a tutti i punti
    rotated_points = np.dot(points, R.T)

    # 6. Scrittura delle coordinate ruotate su file
    with open(output_file, 'w') as f:
        for p in rotated_points:
            f.write(f"{p[0]:15.8f}  {p[1]:15.8f}  {p[2]:15.8f}\n")

    print(f"✅ Geometria elica generata, asse Z inclinato (allineato AD) e salvata in '{output_file}'")

    # Ritorna le coordinate in memoria per l'uso immediato nello script
    return rotated_points

def read_input(filename="input.inp", outputfile="output.txt"):
    """
    Legge i dati da un file di input strutturato e li scrive su un file di output.
    Ritorna un dizionario con tutti i dati letti.
    """
    # === Lettura sicura del file ===
    with open(filename, "r") as f:
        # Prende solo la parte prima del '#' e rimuove gli spazi bianchi
        lines = [line.split('#')[0].strip() for line in f]

    # Rimuove le righe vuote (utile se ci sono righe di soli commenti o spazi)
    lines = [line for line in lines if line]

    # === Lettura variabili scalari in ordine esatto ===
    nsiti = int(lines[0])
    length = float(lines[1])  # Å
    angolo = float(lines[2])
    pitch = float(lines[3])
    t = float(lines[4])  # eV
    PPPflag = int(lines[5])
    mono_flag = int(lines[6])
    bi_flag = int(lines[7])
    multiply_factor = float(lines[8])
    hop_flag = int(lines[9])
    deltat = float(lines[10])  # fs
    points = int(lines[11])
    dynamic_flag = int(lines[12])

    # === Alloca array ===
    u = np.zeros(nsiti)
    esite = np.zeros(nsiti)
    nz = np.zeros(nsiti, dtype=int)
    coord = np.zeros((nsiti, 3))

    # === Lettura per-sito (parte dalla riga 13 in poi) ===
    for i in range(nsiti):
        parts = lines[13 + i].split()
        u[i] = float(parts[0])
        esite[i] = float(parts[1])
        nz[i] = int(parts[2])

    # === Scrittura su file di output ===
    with open(outputfile, "w") as out:
        out.write(f"Numero siti: {nsiti}\n")
        out.write(f"Lunghezza (Å): {length}\n")
        out.write(f"Angolo: {angolo}\n")
        out.write(f"Pitch: {pitch}\n")
        out.write(f"t-hopping (eV): {t}\n")
        out.write(f"Diagonal type (0:Hubbard, 1:PPP): {PPPflag}\n")
        out.write(f"SOC mono_flag: {mono_flag}\n")
        out.write(f"SOC bi_flag: {bi_flag}\n")
        out.write(f"Multiply factor: {multiply_factor}\n")
        out.write(f"Hopping flag (1:NN, 2:exp, 3:incoherent): {hop_flag}\n")
        out.write(f"Δt (fs): {deltat}\n")
        out.write(f"Points: {points}\n")
        out.write(f"Dynamic flag: {dynamic_flag}\n\n")

        out.write("=== Parametri per ogni sito ===\n")
        out.write("  i        u(i)        esite(i)      nz(i)\n")
        for i in range(nsiti):
            out.write(f"{i + 1:3d}  {u[i]:10.5f}  {esite[i]:10.5f}  {nz[i]:3d}\n")

    print(f"\n✅ Input letto da '{filename}' e scritto in '{outputfile}'")

    # === Ritorna i dati come dizionario (utile per altre parti del codice) ===
    return {
        "nsiti": nsiti,
        "length": length,
        "angolo": angolo,
        "pitch": pitch,
        "t": t,
        "PPPflag": PPPflag,
        "mono_flag": mono_flag,
        "bi_flag": bi_flag,
        "multiply_factor": multiply_factor,
        "hop_flag": hop_flag,
        "deltat": deltat,
        "points": points,
        "dyn_flag": dynamic_flag,
        "u": u,
        "esite": esite,
        "nz": nz,
        "coord": coord,
    }
def hubbard_diagonal(nsiti, dimension, basis, esite, u):
    """Generate the diagonal part of the Hubbard Hamiltonian.
               Args:
                nso (int): Number of spin orbitals.
                dimension (int): Dimension of the basis.
                basis (array): Array of basis states in bit representation.
                esite (array): On-site energies.
                u (array): Hubbard U values for each site.

            Returns:
                H (array): Hamiltonian 2D array but with only diagonal elements filled.
            """
    H = np.zeros((dimension, dimension), dtype=complex)
    for i in range(dimension):
        occupazioni = np.zeros(nsiti)
        for j in range(nsiti):
            bool_bit_up = (basis[i] >> (2 * j)) & 1
            occupazioni[j] += bool_bit_up
            bool_bit_down = (basis[i] >> (2 * j + 1)) & 1
            occupazioni[j] += bool_bit_down

        # Inter-site Coulomb interactions
        H[i, i] += occupazioni @ esite
        for j in range(nsiti):
            if occupazioni[j] == 2:
                H[i, i] += u[j]
    return H
def ppp_diagonal(dimension, coord, nsiti, esite, u, nz, basis):
    """Generate the diagonal part of the PPP Hamiltonian.
           Args:
            dimension (int): Dimension of the basis.
            coord (array): Coordinates of the sites.
            nsiti (int): Number of sites.
            esite (array): On-site energies.
            u (array): Hubbard U values for each site.
            nz (array): Number of electrons per site.
            basis (array): Array of basis states in bit representation.

        Returns:
            H (array): Hamiltonian 2D array but with only diagonal elements filled.
        """
    H = np.zeros((dimension, dimension), dtype=complex)
    for i in range(dimension):
        occupazioni = np.zeros(nsiti)
        for j in range(nsiti):
            bool_bit_up = (basis[i] >> (2 * j)) & 1
            occupazioni[j] += bool_bit_up
            bool_bit_down = (basis[i] >> (2 * j + 1)) & 1
            occupazioni[j] += bool_bit_down

        # Inter-site Coulomb interactions
        H[i, i] += occupazioni @ esite
        for j in range(nsiti):
            if occupazioni[j] == 2:
                H[i, i] += u[j]
        for j in range(nsiti - 1):
            for k in range(j + 1, nsiti):
                dist = np.linalg.norm(coord[j] - coord[k])
                vij = 14.397 / dist  # eV·Å
                ppp = (28.794 / (u[k] + u[j]) ** 2) * (nz[j] - occupazioni[j]) * (nz[k] - occupazioni[k])
                H[i, i] += (vij + ppp)
    return H
def hopping_matrix(nsiti, t, length, hop_flag, geom):
    """
    Generate the hopping matrix for the system.

    Args:
        nsiti (int): Number of sites.
        t (float): Hopping parameter.
        length (float): Length of the system.
        hop_flag (int): Flag to determine hopping type.
        geom (ndarray): Array of atomic coordinates, shape (nsiti, 3).

    Returns:
        hopping (ndarray): (2*nsiti, 2*nsiti) hopping matrix.
    """
    hop = np.zeros((nsiti, nsiti))
    distances = np.zeros((nsiti, nsiti))

    # --- Calcolo delle distanze ---
    for i in range(nsiti):
        for j in range(nsiti):
            if i != j:
                dist = np.linalg.norm(geom[i] - geom[j])
                distances[i, j] = dist
    # --- Costruzione matrice hop in base al flag ---
    if hop_flag == 1:
        for i in range(nsiti):
            for j in range(nsiti):
                if distances[i, j] <= 1.001 * length and i != j:
                    hop[i, j] = t

    elif hop_flag == 2:
        for i in range(nsiti):
            for j in range(nsiti):
                if distances[i, j] <= length and i != j:
                    hop[i, j] = t * np.exp(length - distances[i, j])

    elif hop_flag == 3:
        for i in range(nsiti):
            for j in range(nsiti):
                if (i == 2 and j == 3) or (i == 3 and j == 2):
                    hop[i, j] = t
                else:
                    hop[i, j] = 0.1 * t

    elif hop_flag == 4:
        for i in range(nsiti):
            for j in range(nsiti):
                if i == j:
                    continue

                factor = 1.0 if (i, j) in [(1, 2), (2, 1)] else 0.1
                hop[i, j] = factor * t * np.exp(length - distances[i, j])
    # --- Espansione nella matrice completa 2*nsiti x 2*nsiti ---
    hopping = np.zeros((2 * nsiti, 2 * nsiti))
    for i in range(2 * nsiti):
        for j in range(2 * nsiti):
            if i % 2 == j % 2:  # stessi "spin" o sottolivelli
                hopping[i, j] = hop[i // 2, j // 2]
            else:
                hopping[i, j] = 0.0

    return hopping
def is_hermitian(H, tol=1e-12):
    return np.allclose(H, H.conj().T, atol=tol)
def hopping_siti_matrix(coord, nsiti, hop_flag, length, t):
    hop = np.zeros((nsiti, nsiti))
    distances = np.zeros((nsiti, nsiti))

    # --- Calcolo delle distanze ---
    for i in range(nsiti):
        for j in range(nsiti):
            if i != j:
                dist = np.linalg.norm(coord[i] - coord[j])
                distances[i, j] = dist
    # --- Costruzione matrice hop in base al flag ---
    if hop_flag == 1:
        for i in range(nsiti):
            for j in range(nsiti):
                if distances[i, j] <= 1.001 * length and i != j:
                    hop[i, j] = t

    elif hop_flag == 2:
        for i in range(nsiti):
            for j in range(nsiti):
                if distances[i, j] <= length and i != j:
                    hop[i, j] = t * np.exp(length - distances[i, j])

    elif hop_flag == 3:
        for i in range(nsiti):
            for j in range(nsiti):
                if i == 1  or j == 1 or i == nsiti or j == nsiti:
                    hop[i, j] = 0.1 * t
                else:
                    hop[i, j] = t


    elif hop_flag == 4:

        for i in range(nsiti):

            for j in range(nsiti):

                if i == j:
                    continue

                factor = 1.0 if (i, j) in [(1, 2), (2, 1)] else 0.1

                hop[i, j] = factor * t * np.exp(length - distances[i, j])

    return hop
def rotate_matrix(op, eigenvectors, strati, dimension):
    """
    Rotates operator 'op' in the eigenvector basis.

    Args:
        op: (dim, dim) operator matrix
        eigenvectors: (dim, dim) eigenvectors
        strati: number of layers (int)
        dimension: dimension of the system (int)

    Returns:
        op_rotated: (dim, dim, strati) rotated operator
    """
    # Crea array 3D per il risultato
    if strati != 1:
        op_rotated = np.zeros((dimension, dimension, strati), dtype=complex)

        # Rotazione: op_rotated[:,:,s] = V^dagger @ op @ V
        for s in range(strati):
            op_rotated[:, :, s] = np.conj(eigenvectors.T) @ op[:, :, s] @ eigenvectors
    else:
        op_rotated = np.zeros((dimension, dimension), dtype=complex)
        op_rotated = np.conj(eigenvectors.T) @ op @ eigenvectors

    return op_rotated
def compute_soc_mono(basis, coord, nsiti, hop_flag, length, t, nz):
    hop = hopping_siti_matrix(coord, nsiti, hop_flag, length, t)
    mom = np.zeros((nsiti, nsiti, 3), dtype=complex)
    for i in range(nsiti):
        for j in range(nsiti):
            for k in range(3):
                dr = coord[i, k] - coord[j, k]
                mom[i, j, k] = imag * dr * hop[i, j]

    soc_mono = np.zeros((2 * nsiti, 2 * nsiti), dtype=complex)

    for i in range(2 * nsiti):
        sitoi = i // 2
        for j in range(2 * nsiti):
            sitoj = j // 2
            for atomo in range(nsiti):
                if (atomo != sitoj) and (atomo != sitoi):
                    term1 = np.cross(coord[sitoi] - coord[atomo], mom[sitoi, sitoj]) / np.linalg.norm(
                        coord[sitoi] - coord[atomo]) ** 3
                    term2 = np.cross(coord[sitoj] - coord[atomo], mom[sitoi, sitoj]) / np.linalg.norm(
                        coord[sitoj] - coord[atomo]) ** 3
                    if i % 2 == 0:
                        si = 0
                    else:
                        si = 1
                    if j % 2 == 0:
                        sj = 0
                    else:
                        sj = 1
                    soc_mono[i, j] = pf * 0.5 * nz[atomo] * (term1 + term2) @ pauli_matrix[si, sj, :]

    one_electron_soc = tb_to_rs(len(basis), 2 * nsiti, basis, soc_mono)
    return one_electron_soc
def bielectron(basis, nso, op):
    """
    Costruisce la matrice dell'operatore bielettronico su una base determinante.
    basis : lista di determinanti (bitstring interi)
    nso   : numero di spinorbitali
    op    : tensore bielettronico (a,b,c,d) con convenzione ⟨ab|cd⟩
    """
    dim = len(basis)
    op_rs = np.zeros((dim, dim), dtype=complex)

    for n, det_n in enumerate(basis):
        # ciclo su indici dell'operatore a_a† a_b† a_c a_d
        for a in range(nso):
            for b in range(nso):
                for c in range(nso):
                    for d in range(nso):
                        # controlla se c e d sono occupati nel determinante di partenza
                        if btest(det_n, d) and btest(det_n, c):
                            # rimuovi d, c
                            temp = ibclr(ibclr(det_n, d), c)

                            # controlla se a e b sono vuoti
                            if not btest(temp, a) and not btest(temp, b):
                                # aggiungi b, a
                                temp2 = ibset(ibset(temp, b), a)

                                # calcola la fase da permutazione fermionica
                                perm = 0
                                for i in range(min(a, b) + 1, max(a, b)):
                                    if btest(temp, i):
                                        perm += 1
                                for i in range(min(c, d) + 1, max(c, d)):
                                    if btest(det_n, i):
                                        perm += 1
                                phase = 1 if perm % 2 == 0 else -1

                                # trova l'indice m del nuovo determinante
                                m = linear_search(basis, temp2)
                                if m != -1:
                                    op_rs[m, n] += phase * op[a, b, c, d]

    return op_rs
def compute_sso(basis, coord, nsiti, hop_flag, length, t):
    hop = hopping_siti_matrix(coord, nsiti, hop_flag, length, t)
    mom = np.zeros((nsiti, nsiti, 3), dtype=complex)
    for i in range(nsiti):
        for j in range(nsiti):
            for k in range(3):
                dr = coord[i, k] - coord[j, k]
                mom[i, j, k] = imag * dr * hop[i, j]

    soc_bi = np.zeros((2 * nsiti, 2 * nsiti, 2 * nsiti, 2 * nsiti), dtype=complex)

    for a in range(2 * nsiti):
        for b in range(2 * nsiti):
            for c in range(2 * nsiti):
                asito = a // 2
                bsito = b // 2
                csito = c // 2
                if asito != bsito and bsito != csito:
                    term1 = np.cross(coord[asito] - coord[bsito], mom[asito, csito]) / np.linalg.norm(
                        coord[asito] - coord[bsito]) ** 3
                    term2 = np.cross(coord[csito] - coord[bsito], mom[asito, csito]) / np.linalg.norm(
                        coord[csito] - coord[bsito]) ** 3
                    if a % 2 == 0:
                        sa = 0
                    else:
                        sa = 1
                    if c % 2 == 0:
                        sc = 0
                    else:
                        sc = 1
                    soc_bi[a, b, c, b] += pf * 0.5 * (term1 + term2) @ pauli_matrix[sa, sc, :]
    two_electron_soc = bielectron(basis, 2 * nsiti, soc_bi)
    return two_electron_soc
def compute_soo(basis, coord, nsiti, hop_flag, length, t):
    hop = hopping_siti_matrix(coord, nsiti, hop_flag, length, t)
    mom = np.zeros((nsiti, nsiti, 3), dtype=complex)
    for i in range(nsiti):
        for j in range(nsiti):
            for k in range(3):
                dr = coord[i, k] - coord[j, k]
                mom[i, j, k] = imag * dr * hop[i, j]

    soc_bi = np.zeros((2 * nsiti, 2 * nsiti, 2 * nsiti, 2 * nsiti), dtype=complex)

    for a in range(2 * nsiti):
        for b in range(2 * nsiti):
            for c in range(2 * nsiti):
                for d in range(2 * nsiti):
                    asito = a // 2
                    bsito = b // 2
                    csito = c // 2
                    dsito = d // 2
                    if a % 2 == 0:
                        sa = 0
                    else:
                        sa = 1
                    if b % 2 == 0:
                        sb = 0
                    else:
                        sb = 1
                    if c % 2 == 0:
                        sc = 0
                    else:
                        sc = 1
                    if d % 2 == 0:
                        sd = 0
                    else:
                        sd = 1
                    if sa == sc and bsito == dsito:
                        if asito != bsito and csito != bsito:
                            term1 = np.cross(coord[asito] - coord[bsito], mom[asito, csito]) / np.linalg.norm(coord[asito] - coord[bsito]) ** 3
                            term2 = np.cross(coord[csito] - coord[bsito], mom[asito, csito]) / np.linalg.norm(coord[csito] - coord[bsito]) ** 3
                            soc_bi[a, b, c, d] += pf * (term1 + term2) @ pauli_matrix[sb, sd, :]

    two_electron_soc = bielectron(basis, 2 * nsiti, soc_bi)
    return two_electron_soc

def apply_annihilation(state, i):
    """Applica c_i sul determinante. Ritorna (new_state, phase) oppure (None, 0) se annulla."""
    if not btest(state, i):
        return None, 0  # stato vuoto → annulla
    # segno = (-1)^(numero elettroni prima del sito i)
    phase = (-1)**(bin(state & ((1 << i)-1)).count("1"))
    new_state = ibclr(state, i)
    return new_state, phase

def apply_creation(state, i):
    """Applica c_i^† sul determinante. Ritorna (new_state, phase) oppure (None, 0)."""
    if btest(state, i):
        return None, 0  # già occupato → annulla
    phase = (-1)**(bin(state & ((1 << i)-1)).count("1"))
    new_state = ibset(state, i)
    return new_state, phase

def spin_matrices(configs, nso):
    """
    Costruisce S_z, S^+, S^-, S^2 usando operatori scaletta in base di configurazioni.
    """
    # 🔧 fix: configs deve essere una lista per usare .index()
    configs = list(configs)

    dim = len(configs)

    Sz = np.zeros((dim, dim), complex)
    Sp = np.zeros((dim, dim), complex)
    Sm = np.zeros((dim, dim), complex)
    S2 = np.zeros((dim, dim), complex)

    # ----------------------------
    #     S_z (diagonale)
    # ----------------------------
    for a, state in enumerate(configs):
        n_up = 0
        n_down = 0
        for i in range(0, nso, 2):
            if btest(state, i):     n_up += 1
            if btest(state, i+1):   n_down += 1
        Sz[a, a] = 0.5 * (n_up - n_down)

    # ----------------------------
    #     S^+ = sum_i c†_{i↑} c_{i↓}
    # ----------------------------
    for a, state in enumerate(configs):
        for site in range(nso // 2):
            down = 2*site + 1
            up   = 2*site

            state2, ph1 = apply_annihilation(state, down)
            if state2 is None:
                continue

            state3, ph2 = apply_creation(state2, up)
            if state3 is None:
                continue

            try:
                b = configs.index(state3)
                Sp[b, a] += ph1 * ph2
            except ValueError:
                pass

    # ----------------------------
    #     S^- = sum_i c†_{i↓} c_{i↑}
    # ----------------------------
    for a, state in enumerate(configs):
        for site in range(nso // 2):
            down = 2*site + 1
            up   = 2*site

            state2, ph1 = apply_annihilation(state, up)
            if state2 is None:
                continue

            state3, ph2 = apply_creation(state2, down)
            if state3 is None:
                continue

            try:
                b = configs.index(state3)
                Sm[b, a] += ph1 * ph2
            except ValueError:
                pass

    # ----------------------------
    #     S^2 = S^- S^+ + S_z(S_z + 1)
    # ----------------------------
    S2 = Sm @ Sp + Sz @ (Sz + np.identity(dim))

    return Sz, S2