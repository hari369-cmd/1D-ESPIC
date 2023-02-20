#include "define.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <mpi.h>

void OutPut::Initialize(){

    for (int j = 0; j < 100; j++)
    {
        char filename[32];
        snprintf(filename, sizeof(char) * 32, "trace%i.dat", j);
        trace[j] = fopen(filename, "w");
    }
};

void OutPut::Particle_Info(ConstParam cp,Particle *e,Particle *Ion,int k){

    std::ostringstream a;
    a << k;

    std::string name = "./Electron_Data/Electron_";
    std::string name1 = "./Ion_Data/Ion_";
    std::string filename;
    std::string filename1;

    filename = name + a.str() + ".dat";
    filename1 = name1 + a.str() + ".dat";

    std::ofstream outFile_ved(filename);
    std::ofstream outFile_vid(filename1);

    for (int j = 0; j < cp.N; j++)
	{
        outFile_ved << e[j].x << "\t" << e[j].vhx << "\t" << e[j].vhy << "\t" << e[j].vhz << "\n";
        outFile_vid << Ion[j].x << "\t" << Ion[j].vhx << "\t" << Ion[j].vhy << "\t" << Ion[j].vhz << "\n";
	}

  //fclose(filename);
  //fclose(filename1);
};

void OutPut::Grid_Info(ConstParam cp,GridPoint *gp,int k){

    std::ostringstream a;
    a << k;

    std::string name = "./Grid_Data/Grid_";
    std::string filename;

    filename = name + a.str() + ".dat";

    std::ofstream outFile_grid(filename);

    for (int j = 0; j <= cp.NG; j++)
    {
        outFile_grid << gp[j].X << "\t" << gp[j].Phi << "\t" << gp[j].E << "\t" << gp[j].denseE << "\t" << gp[j].denseI << "\n";
    }

    //fclose(filename);
};

void OutPut::trajectory(ConstParam cp, Particle *e, int k, int rank, int num){
    /**
    int temp = cp.N / part_per_file;

    for(int j = 0; j < temp; j++)
    {
        char filename[32];
        snprintf(filename, sizeof(char) * 32, "trace%i.dat", j);
        f = fopen(filename, "a+");

        for(int m = 0; m < part_per_file; m++)
        {
            int loc = m + j * part_per_file;
            fprintf(f, "%e %e %e %e %e \t", e[loc].x,e[loc].vhx,e[loc].vhy,e[loc].E,e[loc].Bz);
        }
        fprintf(f, "\n");
        fclose(f);
    }
    **/

    int part_per_rank = cp.N / num;
    int file_per_rank = 25;
    int part_per_file = part_per_rank / file_per_rank;
    int start = rank * part_per_rank;
    int end = 0;
    int file_start = rank * file_per_rank;
    int file_end = file_start + file_per_rank;

    for(int j = file_start; j < file_end; j++)
    {
        end = start + part_per_file;
        for(int m = start; m < end; m++)
        {
            fprintf(trace[j],"%e %e %e %e \t", e[m].x,e[m].vhx,e[m].vhy,e[m].Bz);
        }
        fprintf(trace[j], "\n");
        start = end;
    }
};

void OutPut::History_Info(int step,ConstParam cp,Particle *e,Particle *Ion,GridPoint *gp,HalfGridPoint *gph){

	left_frag = 0;
	right_frag = 0;
  energyE= 0;
  energyI=0;

    for (int j = 0; j < cp.N; j++)
	{
		if (e[j].over_frag == -1)
		{
		    left_frag++;
		}
		else if (e[j].over_frag == 1)
		{
			right_frag++;
		}
		energyE += 0.5 * e[j].m * (e[j].vhx * e[j].vhx + e[j].vhy * e[j].vhy + e[j].vhz * e[j].vhz);
		energyI += 0.5 * Ion[j].m * (Ion[j].vhx * Ion[j].vhx + Ion[j].vhy * Ion[j].vhy + Ion[j].vhz * Ion[j].vhz);
	}

    printf("%e %e %e %d %d \n", step * cp.DT, energyE, energyI, left_frag, right_frag);

};
