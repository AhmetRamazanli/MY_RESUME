# Psi4 plugins: DF-MP2, DF-CIS and a CI configuration generator

These are C++ plugins I wrote for the [Psi4](https://psicode.org) quantum chemistry package during my time at the Chemistry Department of Hacettepe University. They are built on Psi4's plugin interface and a tensor library (`tensors.h`) used in our research group. Two-electron integrals are handled with density fitting (DF).

## plugin-mp2.cc

Computes the MP2 (second-order Møller–Plesset perturbation theory) correlation energy on top of a converged SCF calculation. It reads the MO energies and the DF three-index integrals (Q|ia), assembles the (ia|jb) integrals, applies the orbital-energy denominators to get the T2 amplitudes, and prints the SCF, correlation and total MP2 energies. Both RHF and UHF references are handled.

## plugin_cis.cc

Computes CIS (Configuration Interaction Singles) excitation energies. The CIS Hamiltonian is never stored explicitly; sigma vectors (σ = Hc) are formed directly from the DF integrals. For the RHF path the lowest roots are converged with a Davidson solver, and for UHF the spin-resolved (alpha/beta) sigma builds are in place.

A note on attribution: the Davidson diagonalization routine is adapted from existing reference code and is not my own implementation. My work here is the DF integral handling, the diagonal Hamiltonian and initial guess setup, and the sigma-vector construction.

## configmatrix.cc

A small standalone program that enumerates electron configurations for a CI expansion. It computes the number of configurations C(NMO, Nelec), builds a vertex-weight table from binomial coefficients, and decodes each configuration index into its list of occupied orbitals (written to `output.dat`).

## Building

The two plugins compile against Psi4's C++ plugin interface. `configmatrix.cc` is standalone and only needs the tensor library.
