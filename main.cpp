#include "define.hpp"
#include <stdio.h>
#include <complex.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <mpi.h>

int main(int argc, char **argv)
{
    int p, id;
    MPI_Init ( NULL, NULL );
    MPI_Comm_size ( MPI_COMM_WORLD, &p );
    MPI_Comm_rank ( MPI_COMM_WORLD, &id );

    // Declaration of the constants
    ConstParam cp;

    // Declaration of output file
    OutPut op;
    //op.Initialize(); // Particle trajectory output initialization

    // Electron generation
    Particle *e;
    e = new Particle[cp.N];

    // Ion generation
    Particle *Ion;
    Ion = new Particle[cp.N];

    // Creattion of the full grid
    GridPoint *gp;
    gp = new GridPoint[cp.NG + 1];

    // Creation of the half grid
    HalfGridPoint *gph;
    gph = new HalfGridPoint[cp.NG];

    // Initialization of electron parameters
    SetParticle spe;
    spe.SetMQ(cp, e, cp.ME, cp.QE);
    spe.SetPosition(cp, e);
    spe.InitVe(cp, e, cp.VE, cp.VE1, cp.ME);
    printf("Electron loading succesful\n");

    // Initialization of ion parameters
    SetParticle spi;
    spi.SetMQ(cp, Ion, cp.MI, cp.QI);
    // Match the initial position of the ion with the electron
    for (int i = 0; i < cp.N; i++)
    {
        Ion[i].x = e[i].x;
    }
    spi.InitVi(cp, Ion, cp.VI, cp.MI);
    printf("Ion loading succesful\n");

    // Custom speed definition (edit only for special cases)
    //for (int i = 0; i < cp.N; i++)
    //{
    //    e[i].vy = 0;
    //    e[i].vz = 0;

    //    Ion[i].vy = 0;
    //    Ion[i].vz = 0;
    //}

    // Initializing the electric field
    SolveField sf;
    sf.InitField(cp, e, gp);
    printf("Field initialization succesful");

    // Coordinates of the half grid
    for (int i = 0; i < cp.NG; i++)
    {
        gph[i].X = gp[i].X + cp.DX * 0.5;
    }

    // Initializing the particle density
    sf.InitRho(cp, gp, gph);
    printf("Particle density initialization succesful\n");

    // Charge density calculation (second order weighting)
    sf.CalcDensity_2(cp, e, gp);
    sf.CalcDensity_2(cp, Ion, gp);

    // Initializing the solver. Matrix definition of Thomas method
    sf.SetSolver(cp);
    printf("Solver initialization succesful\n\n\n");

    // Potential calculation from Poisson equation
    sf.CalcPhi(cp, gp);
    // Calculate electric field from electric potential
    sf.CalcE(cp, gp, gph);

    // Initializing the forces
    SolveMotion sm;

    // Initialization of the electric field that the particle feels
    sm.InitE(cp, e);
    sm.InitE(cp, Ion);

    // Calculate the electric field that the particle feels (fourth order weighting)
    // Can be changes to second order weighting
    sm.weighting_4(cp, gp, e);
    sm.weighting_4(cp, gp, Ion);
    sm.weighting_4vhxe(cp, gp, e);
    sm.weighting_4vhxi(cp, gp, Ion);

    // Calculate the particle velocity in half step before Buneman Boris method
    sm.setV(cp, e);
    sm.setV(cp, Ion);

    // Flag definition
    int nan_flag = 0;

    //MAIN LOOP
    for (int i = 0; i <= cp.NT; i++)
    {
        if (id == 0)
        {
            sm.defineB(cp, e, i);
            sm.defineB(cp, Ion, i);

            sm.accel(cp, e);
            sm.accel(cp, Ion);

            sm.move(cp, e);
            sm.move(cp, Ion);

            sm.boundary_left_adsorp_high(cp, e, Ion);
            sm.boundary_right_adsorp_high(cp, e, Ion);

            sf.InitRho(cp, gp, gph);
            sf.CalcDensity_2(cp, e, gp);
            sf.CalcDensity_2(cp, Ion, gp);
            sf.CalcPhi(cp, gp);
            sf.CalcE(cp, gp, gph);

            sm.InitE(cp, e);
            sm.InitE(cp, Ion);

            sm.weighting_4(cp, gp, e);
            sm.weighting_4(cp, gp, Ion);
            sm.weighting_4vhxe(cp, gp, e);
            sm.weighting_4vhxi(cp, gp, Ion);


            if (i == 0)
            {
                printf("Cell size: %E [m]\n", cp.DX);
                printf("Time step size: %E [s]\n", cp.DT);
                printf("Total number of iterations: %d\n", cp.NT);
                printf("Iteration of pulse injection: %d\n", cp.BSS);
            }

            if (i % cp.FREQ == 0)
            {
                printf("Iteration: %d\n", i);
                op.Particle_Info(cp,e,Ion,i);
			          op.Grid_Info(cp,gp,i);
                op.History_Info(i,cp,e,Ion,gp,gph);

                // Calculation ends when Nan comes out
                for (int j = 0; j < cp.N; j++)
                {
                    if (isnan(e[j].x) || isnan(Ion[j].x))
                        nan_flag = 1;
                }
                if (nan_flag == 1){
                    printf("Particle position out of bounds \n");
                    break;
                }
            }
        }

        /**
        //Particle Trajectory output
        if (i % cp.FREQ == 0)
        {
            op.trajectory(cp ,e, i, id, p);
        }
        **/
    }

    MPI_Finalize();
	return 0;
}
