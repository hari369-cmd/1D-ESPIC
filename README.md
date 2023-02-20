1D-ESPIC developed by Hari

About
    This is an in-house particle-in-cell (PIC) code built in C++ with support for parallelization through     OPENMP and MPI for simulating 1D         electrostatic plasma in a Traveling Magnetic Field Thruster (TMFT) environment. This code supports 2 species: electrons, and ions.

    This code is used in the following papers:
    (1) https://www.jstage.jst.go.jp/article/tastj/18/6/18_18.317/_article/-char/ja/
    (2) https://aip.scitation.org/doi/abs/10.1063/5.0065665


ConstParam_pulse.cpp
    Definition of the crucial simulation constants.

define.hpp
    Global class definitions.

main.cpp
    Contains the class declarations and the main time loop.

SetParticle.cpp
    Loading the particle species into the domain and set initial conditions.

SolveField.cpp
    Initializing the solver grid and the Poisson equation solver.

SolveMotion.cpp
    Time stepping of the particles using the Boris method and the evaluation of particle/field boundary conditions.

Output.cpp
    Outputs the data generated from the simulation.



How to run the code?

1) Make sure MPI and GCC compiler are installed.
2) Run the Makefile in Terminal or CMD.
