
#define _USE_MATH_DEFINES
#include <math.h>
#include "define.hpp"


const double ConstParam::ME = 9.10938356E-31;    // Mass of an electron [kg]
const double ConstParam::MI = 1.672621898E-27;   // Mass of a H+ ion [kg]
const double ConstParam::QE = -1.6021766208E-19; // Charge of an electron [C]
const double ConstParam::QI = +1.6021766208E-19; // Charge of a H+ ion [C]
const double ConstParam::EPSI = 8.85418782E-12;  // Vacuum permittivity
const double ConstParam::K = 1.3806503E-23;      // Boltzmann Constant [J/K]
const int ConstParam::C = 299792458;          // Speed of light [m/s]

const double ConstParam::L = 0.2;                      // Length of the domain (x) [m]
const double ConstParam::DS = M_PI * 0.0125 * 0.0125;  // Fictional area of the domain (y & z) [m^2]
const double ConstParam::T_end = 8e-6;                 // End time for the simulation [s]
const double ConstParam::T_BSS = 4.5e-6;               // Start time of pulse injection [s]

const double ConstParam::BV = 3E+4;       // Velocity of the moving magnetic pulse [m/s]
const double ConstParam::BI = 0.04;       // Magnetic field intensity [T]
const double ConstParam::BW = 0.006;      // Magnetic pulse width [m]
const int ConstParam::BSS = T_BSS / DT;   // Time step of the pulse injection
const double ConstParam::B_MN = 0;        // Background magnetic field [T]

const double ConstParam::CFL = 100;      // Courant–Friedrichs–Lewy condition

const double ConstParam::W_c = BI * QI / ME;   // Cyclotron frequency [1/s]
const double ConstParam::DT = 1.0 * (1 / W_c); // Time discretization [s]


const int ConstParam::NT = T_end / DT;      // Number of calculation steps
const int ConstParam::NG = 1000;            // Number of grid points
const double ConstParam::DX = L / NG;       // Spatial discretization [m]
//const double ConstParam::DT = CFL * DX / C; // Time discretization [s]
const int ConstParam::NAP = 40000;          // Super particles per particle
const int ConstParam::N = 120000;           // Total number of super particles (same number of electrons and ions)
const int ConstParam::NV = 3000;            // Particle initial velocity division number
const int ConstParam::FREQ = 10;            // Output frequency

const double ConstParam::alpha = 0;  // Percentage of electrons at a different temperature
const double ConstParam::VE1 = 0.5;  // Initial electron temperature in y and z velocities (VE1=VE when alpha=0)
const double ConstParam::VE = 0.5;   // Initial electron temperature for x-velocity [eV]
const double ConstParam::VI = 0.1;   // Initial ion temperature [eV]

const int ConstParam::margin = 2;   // Number of cells between the computing edge and the wall that absorbs particles
