import numpy as np
def generate_basis(electrons):
    """Generate basis states in bit representation.
       arranged in order of increasing spin.

    Args:
        sites (int): Number of lattice sites.
        electrons (int): Number of electrons."""
    # ==== Parametri di input ====
    nso = electrons * 2  # numero di spinorbitali

    # ==== Limiti del ciclo ====
    n_max = sum(2 ** i for i in range(nso - electrons, nso))
    n_min = sum(2 ** i for i in range(electrons))

    nf = 0  # contatore configurazioni
    usefull = []  # lista delle configurazioni utili

    # ==== Loop principale ====
    with open('../hubbard_soc/out1.txt', 'w') as f1, open('../hubbard_soc/out2.txt', 'w') as f2:
        for n in range(n_min, n_max + 1):
            count = 0
            config = 0
            a = 0
            b = 0
            array = ['0'] * nso

            # ciclo sugli orbitali
            for i in range(nso):
                bool_bit = (n >> i) & 1  # equivale a btest(n, i)
                if bool_bit:
                    array[i] = '1'
                    count += 1
                    if i % 2 == 0:
                        a += 1
                    else:
                        b += 1
                else:
                    array[i] = '0'

            spin = (a - b) * 0.5

            if count == electrons:
                config = sum(2 ** i for i in range(nso) if array[i] == '1')


                # Aggiungi alla lista 'usefull' come tuple (config, spin, array)
                usefull.append((config, spin, ''.join(array)))
                nf += 1

    print(f"Numero totale di configurazioni valide: {nf}")

    # ==== Ordinamento per spin crescente ====
    usefull_sorted = sorted(usefull, key=lambda x: x[1])

    # ==== Stampa finale ordinata ====
    with open('../hubbard_soc/configurations.txt', 'w') as f3:
        for config, spin, array in usefull_sorted:
            f3.write(f"{array}  spin={spin:+.1f}  config={config}\n")
    configs = np.array([u[0] for u in usefull])
    spins = np.array([u[1] for u in usefull])
    return configs, spins, nf

def btest(n, i):
    return ((n >> i) & 1) == 1

def ibset(n, i):
    return n | (1 << i)

def ibclr(n, i):
    return n & ~(1 << i)

# Funzione di ricerca binaria (assumendo basis ordinata)
def binary_search(arr, val, left=0, right=None):
    if right is None:
        right = len(arr) - 1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == val:
            return mid
        elif arr[mid] < val:
            left = mid + 1
        else:
            right = mid - 1
    return -1  # non trovato
def tb_to_rs(dim, nso, basis, op_tb):
    op = np.zeros((dim, dim))
    for j in range(dim):  # colonna
        jstate = basis[j]
        for iso in range(nso):  # indice creazione
            for jso in range(nso):  # indice annichilazione
                if btest(jstate, jso):
                    istate = ibclr(jstate, jso)
                    if not btest(istate, iso):
                        istate = ibset(istate, iso)

                        i = binary_search(basis, istate)
                        if i != -1:
                            # calcolo fase
                            if iso == jso:
                                phase = 1.0
                            else:
                                step = -1 if jso > iso else 1
                                conta = 0
                                for a in range(jso + step, iso, step):
                                    if btest(istate, a):
                                        conta += 1
                                phase = 1.0 if conta % 2 == 0 else -1.0

                            op[i, j] += phase * op_tb[iso, jso]
    return op
def read_input(filename="input.inp", outputfile="output.txt"):
    """
    Legge i dati da un file di input (stile Fortran) e li scrive su un file di output.
    Ritorna un dizionario con tutti i dati letti.
    """
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # === Lettura variabili scalari ===
    nsiti       = int(lines[0])
    length      = float(lines[1])  # Å
    t           = float(lines[2])  # eV
    PPPflag     = int(lines[3])
    mono_flag   = int(lines[4])
    bi_flag     = int(lines[5])
    hop_flag     = int(lines[6])
    # === Alloca array ===
    u     = np.zeros(nsiti)
    esite = np.zeros(nsiti)
    nz    = np.zeros(nsiti, dtype=int)
    coord = np.zeros((nsiti, 3))

    # === Lettura per-sito ===
    for i in range(nsiti):
        parts = lines[7 + i].split()
        u[i] = float(parts[0])
        esite[i] = float(parts[1])
        nz[i] = int(parts[2])

    # === Scrittura su file di output ===
    with open(outputfile, "w") as out:
        out.write(f"Numero siti: {nsiti}\n")
        out.write(f"Lunghezza (Å): {length}\n")
        out.write(f"t (eV): {t}\n")
        out.write(f"Diagonal type: {PPPflag}\n")
        out.write(f"mono_flag: {mono_flag}\n")
        out.write(f"bi_flag: {bi_flag}\n")
        out.write(f"hopping flag: {hop_flag}\n\n")
        out.write("=== Parametri per ogni sito ===\n")
        out.write("  i        u(i)        esite(i)      nz(i)\n")
        for i in range(nsiti):
            out.write(f"{i+1:3d}  {u[i]:10.5f}  {esite[i]:10.5f}  {nz[i]:3d}\n")

    print(f"\n✅ Input letto da '{filename}' e scritto in '{outputfile}'")

    # === Ritorna i dati come dizionario (utile per altre parti del codice) ===
    return {
        "nsiti": nsiti,
        "length": length,
        "t": t,
        "PPPflag": PPPflag,
        "mono_flag": mono_flag,
        "bi_flag": bi_flag,
        "u": u,
        "esite": esite,
        "nz": nz,
        "coord": coord,
        "hop_flag": hop_flag
    }
def hubbard_diagonal(nso, dimension, basis, esite, u):
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
    H = np.zeros((dimension,dimension), dtype = complex)
    for i in range(dimension):
        for j in range(nso,2):
            sito = (j+2)//2
            bool_bit = (basis[i] >> j) & 1
            if bool_bit:
                H[i,i] += esite[sito-1]
                # check if the opposite spin orbital is also occupied
            bool_bit_opposite = (basis[i] >> j) & 1
            if bool_bit_opposite:
                H[i,i] += esite[sito-1]
            if bool_bit and bool_bit_opposite:
                H[i,i] += u[sito-1]
    return H

def ppp_diagonal(dimension,coord, nsiti, esite, u, nz, basis):
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
    H = np.zeros((dimension,dimension), dtype = complex)
    for i in range(dimension):
        occupazioni = np.zeros(nsiti)
        for j in range(nsiti):
            bool_bit_up = (basis[i] >> (2*j)) & 1
            occupazioni[j] += bool_bit_up
            bool_bit_down = (basis[i] >> (2*j + 1)) & 1
            occupazioni[j] += bool_bit_down
            if bool_bit_up:
                H[i,i] += esite[j]
            if bool_bit_down:
                H[i,i] += esite[j]
            if bool_bit_up and bool_bit_down:
                H[i,i] += u[j]
        # Inter-site Coulomb interactions

        for j in range(nsiti):
            for k in range(j+1, nsiti):
                dist = np.linalg.norm(coord[j] - coord[k])
                vij = 14.397 / dist  # eV·Å
                ppp = (28.794/(u[i]+u[j]**2)) * (nz[j]-occupazioni[j]) * (nz[k]-occupazioni[k])
                H[i,i] += vij * ppp
    return H
def hopping_matrix(nsiti, t, length, hop_flag, geom):
        """Generate the hopping matrix for the system.
        Args:
        nsiti (int): Number of sites.
        t (float): Hopping parameter.
        length (float): Length of the system.
        hop_flag (int): Flag to determine hopping type.
        """
    hop = np.zeros((nsiti, nsiti))
    distances = np.zeros((nsiti, nsiti))
    for i in range(nsiti):
        for j in range(nsiti):
            if i != j:
                dist = np.linalg.norm(geom[i] - geom[j])
                distances[i,j] = dist
    if hop_flag == 1:
        for i in range(nsiti):
            for j in range(nsiti):
                if distances[i,j] <= length and i != j:
                    hop[i,j] =
    if hop_flag == 2:
        for i in range(nsiti):
            for j in range(nsiti):
                if distances[i,j] <= length and i != j:
                    hop[i,j] = t * np.exp(length - distances[i,j])
    if hop_flag == 3:
        for i in range(nsiti):
            for j in range(nsiti):
                if (i==2 and j==3) or (i==3 and j==2):
                    hop[i,j] = t
                else:
                    hop[i,j] = 0.1 * t
    if hop_flag == 4:
    for i in range(nsiti):
        for j in range(nsiti):
            if (i == 2 and j == 3) or (i == 3 and j == 2):
                hop[i, j] = t * np.exp(length - distances[i, j])
            else:
                hop[i, j] = 0.1 * t*np.exp(length - distances[i, j])
    return hopping
