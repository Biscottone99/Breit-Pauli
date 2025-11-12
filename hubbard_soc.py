import numpy as np
import module as mod


if __name__ == "__main__":
    #================Constants============================================
    imag = 1j
    me = 9.1093837015e-31
    gs = 2.00231930436256
    e = 1.602176634e-19
    e0 = 8.8541878128e-12
    cl = 299792458
    pf = ((gs * e ** 2) / (8 * np.pi * e0 * me * cl ** 2)) * 10.0e10
    #================Initialize our system===================================================
    input = mod.read_input("input.inp", "output.txt")
    basis, sz, dimension = mod.generate_basis(input["nsiti"])
    coords = np.loadtxt('geom.dat')
    nso = 2*input["nsiti"]
    #================Generate Hamiltonian===================================================
    hamiltonian = np.zeros((dimension,dimension), dtype = complex)
    hamiltonian += mod.hubbard_diagonal(nso, dimension, basis, input["esite"], input["u"])

