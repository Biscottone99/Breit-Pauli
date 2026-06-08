import numpy as np
import sys
from scipy.linalg import eigh
import module as mod
import plotting_module as plot
import subprocess
if __name__ == "__main__":
    #================Constants============================================
    imag = 1j
    me = 9.1093837015e-31
    gs = 2.00231930436256
    e = 1.602176634e-19
    e0 = 8.8541878128e-12
    cl = 299792458
    pf = ((gs * e ** 2) / (8 * np.pi * e0 * me * cl ** 2)) * 10.0e10
    check_flag = False
    #================Initialize our system===================================================
    input = mod.read_input("input.inp", "output.txt")
    basis, sz, dimension = mod.generate_basis(input["nsiti"])
    with open("basis.txt", 'w') as f:
        for i in range(dimension):
            f.write(f'State {i}: {basis[i]:0{input["nsiti"]*2}b}\n')
    coords = mod.generate_rotated_geometry(input["nsiti"], input["length"], input["pitch"], input["angolo"], output_file="geom_rotate.dat")
    plot.plot_geom(coords, input["nz"])
    with open('output.txt', 'a') as f:
        f.write("\n")
        f.write('\nAtomic Coordinates:\n')
        for i in range(input["nsiti"]):
            f.write(f'Atom {i + 1}: x={coords[i, 0]:.6f} y={coords[i, 1]:.6f} z={coords[i, 2]:.6f}\n')
        f.write("\n")
    nso = 2 * input["nsiti"]
    hop_nsiti = mod.hopping_siti_matrix(coords, input["nsiti"], input["hop_flag"], input["length"], input["t"])
    with open('output.txt', 'a') as f:
        f.write('Hopping matrix between sites (eV):\n')
        for i in range(input["nsiti"]):
            for j in range(input["nsiti"]):
                f.write(f'{hop_nsiti[i,j]:.6f} ')
            f.write('\n')
        f.write('\n')
#================Generate Hamiltonian===================================================
    hamiltonian = np.zeros((dimension,dimension), dtype = complex)
    if input["PPPflag"] == 0:
        hamiltonian = mod.hubbard_diagonal(input["nsiti"], dimension, basis, input["esite"], input["u"])
    elif input["PPPflag"] == 1:
        hamiltonian = mod.ppp_diagonal(dimension, coords, input["nsiti"], input["esite"], input["u"], input["nz"], basis)
    hopping_so = mod.hopping_matrix(input["nsiti"], input["t"], input["length"], input["hop_flag"],coords)
    hopping = mod.tb_to_rs(dimension, nso, basis, hopping_so)
    if check_flag:
        plot.plot_heatmap_real(hopping_so, 'Hopping spin orbitals')
    hamiltonian -= hopping
    if input["mono_flag"] == 1:
        soc_mono = mod.compute_soc_mono(basis, coords, input["nsiti"], input["hop_flag"], input["length"], input["t"], input["nz"])
        hamiltonian += soc_mono * 10**input["multiply_factor"]
    if input["bi_flag"] == 1:
        soo = mod.compute_soo(basis, coords, input["nsiti"], input["hop_flag"], input["length"], input["t"])
        sso = mod.compute_sso(basis, coords, input["nsiti"], input["hop_flag"], input["length"], input["t"])
        hamiltonian += (soo+sso)*10**input["multiply_factor"]
    #=========================Generate some usefull operators================================
    sz, s2 = mod.spin_matrices(basis, nso)
    number_operator = np.zeros((dimension,dimension,nso), dtype = complex)
    for i in range(dimension):
        for j in range(nso):
            if mod.btest(basis[i],j):
                number_operator[i,i,j] = 1.0
    dipole_moment = np.zeros((dimension, dimension, 3), dtype = complex)
    cariche = np.zeros((dimension, dimension, input["nsiti"]), dtype=complex)

    for i in range(dimension):
        for j in range(0, nso, 2):  # ciclo a coppie (spin up / spin down)
            sito = j // 2

            # carica = nz - (n_up + n_down)
            carica = input["nz"][sito] - (
                    number_operator[i, i, j] +
                    number_operator[i, i, j + 1]
            )

            cariche[i, i, sito] = carica

            # aggiorno dipolo (x,y,z)
            dipole_moment[i, i, 0] += coords[sito, 0] * carica
            dipole_moment[i, i, 1] += coords[sito, 1] * carica
            dipole_moment[i, i, 2] += coords[sito, 2] * carica
    #========================Diagonalize Hamiltonian=========================================
    #Let's check if is hermitian
    if check_flag:
        plot.plot_heatmap_real(np.real(hopping), 'Hopping RS')
        plot.plot_heatmap_cplx(hamiltonian,'Hamiltonian')
    bool = mod.is_hermitian(hamiltonian, 1e-8)
    if(bool):
        print('Hamiltonian is Hermitian')
        plot.plot_heatmap_cplx(hamiltonian,'Hamiltonian')
    else:
        print('Hamiltonian is not Hermitian')
        print('Blocking the program')
        sys.exit()

    eigenvalue, eigenvectors = eigh(hamiltonian)
    eigenvalue -= eigenvalue[0]
    D = np.diag(eigenvalue)
    res = eigenvectors.conj().T @ hamiltonian @ eigenvectors - D
    print("Residuo max:", np.max(np.abs(res)))
    sz_rot = mod.rotate_matrix(sz,eigenvectors,1,dimension)
    s2_rot = mod.rotate_matrix(s2,eigenvectors,1,dimension)
    number_operator_rot = mod.rotate_matrix(number_operator,eigenvectors,nso,dimension)
    rotated_dipole = mod.rotate_matrix(dipole_moment,eigenvectors,3,dimension)
    #========================Print some results=============================================
    with open('output.txt', 'a') as f:
        f.write('Eigenvalues (eV):\n')
        for energy, sz, s2 in zip(eigenvalue, np.real(np.diag(sz_rot)), np.real(np.diag(s2_rot))):
            f.write(f'Energy: {energy:.10f} eV, <Sz>: {sz:.6f}, <S2>: {s2:.6f}\n')
        f.write('\n')

    with open('output.txt', 'a') as f:
        f.write('Number operator:\n')
        for i in range(20):
            row_str = " ".join(f"{x:.6f}" for x in np.real(number_operator_rot[i, i, :]))
            f.write(f"State: {i}, <n>: {row_str}\n")
        f.write('\n')
if input["dyn_flag"] == 0:
    sys.exit()
#========================Dynamic part===========================================================
psi0 = rotated_dipole[0,:,2] #prima riga del vettore dipolo di transizione lungo z
psi0 /= np.linalg.norm(psi0)
with open('output.txt', 'a') as f:
    f.write('Initial State (Dipole along z):\n')
    for i in range(dimension):
        f.write(f'{np.real(psi0[i]):.6f} + {np.imag(psi0[i]):.6f}j\n')
    f.write('\n')

density_matrix = np.outer(np.conj(psi0), psi0)
density_matrix.tofile('rho.bin')
spindensity = ( (number_operator_rot[:,:,6] - number_operator_rot[:,:,7]) - (number_operator_rot[:,:,0] - number_operator_rot[:,:,1]) )
spindensity.tofile('spindensity.bin')
eigenvalue.tofile('eigenvalues.bin')


rho_file = "rho.bin"
eigen_file = "eigenvalues.bin"
props_files = ["spindensity.bin"]  # lista dei file dei tensori di proprietà
n_prop = len(props_files)

subprocess.run(['ifx', 'Unitary.f90', '-o', 'unitary.e', '-qmkl', '-qopenmp'], check=True)

program_input = "\n".join([
    str(dimension),
    rho_file,
    eigen_file,
    str(n_prop),
] + props_files + [
    str(input["deltat"]),
    str(input["points"])
]) + "\n"

subprocess.run(['./unitary.e'], input=program_input, text=True, check=True)
subprocess.run(['mv', 'prop1.dat', 'spinpol_evolution.dat'], check=True)
#========================End of unitary evolution===========================================
# Now we can plot some results from the unitary evolution
data = np.loadtxt('spinpol_evolution.dat', usecols=(0, 1))

time = data[:, 0] * 1e-3     # converto in ns
spinpol = data[:, 1]  * 50  # convert in percentual
plot.plot_curve(time, spinpol, 'Spin Polarization Evolution', 'Time (ns)', 'Spin Polarization (%)')