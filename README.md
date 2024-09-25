These codes are written in Fortran 77 with the collaboration of Prof Sonjoy Majumder, Indian Institute of Technology Kharagpur, India.
There are several different codes for different levels of coupled cluster and Hartree Fock calculations
1)gauss_mpi.F: it generates the Hartree Fock orbitals to be used in the coupled cluster calculations as reference states of orbitals.
Prof S. Majumder's group primarily wrote this code. Dr. Sourav Roy had updated it to implement Isotopic shift calculations. Dr Roy
also updates this code for the MPI version. 
2) ccsd_mpi_11sym.F: This code uses orbitals generated in gauss_mpi.F to calculate the correlations imposed by the electrons according to coupled cluster theory 
inside the core.
3) ccea_mpi_11sym.F: This code uses updated orbitals generated with ccsd_mpi_11sym.F code and calculate correlations for valence orbitals.
4) ccsd_mpi_SMS_11sym.F and ccea_mpi_SMS_11sym.F: These codes calculate Isotopic shifts for atomic systems.
5) ccsd_mpi_new.F and ccea_mpi_new.F: This code calculates up to 9 symmetry orbitals.
Dr Sourav Roy updated These codes to Isotopic Shift calculations and MPI versions.
