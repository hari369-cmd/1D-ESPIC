#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <mpi.h>

class ConstParam
{
public:

  static const int C;

  static const double ME, MI, QE, QI, EPSI, K;

  static const double alpha, VE1, VE, VI;

  static const double L, DX, DS, DT, W_c, T_end, T_BSS;

  static const double CFL;

  static const int NT, NG, NAP, N, NV, FREQ;

  static const double BV, BI, BW, B_MN;
  static const int BSS;

  static const int margin;
};

class Particle
{
public:
  double x, vx, vy, vz, vhx, vhy, vhz, E, Bx, By, Bz, m, q, x_ini, kee, kei, kex, kix;
  int over_frag, bound_frag;
  int near;
};

class GridPoint
{
public:
  double X, Rho, Phi, E, denseE, denseI, TKE, elect, ion, TKE1, Vavg1_e, KE, KE1, den_e, den_i, Vavg_e, Vavg_i, Vavg1_i, vtherm;
};

class HalfGridPoint
{
public:
  double X, E, Rho;
};

//set initial particle property
class SetParticle
{
private:
public:
  void SetMQ(ConstParam, Particle *, double, double);
  void SetPosition(ConstParam, Particle *);
  void InitVe(ConstParam, Particle *, double, double, double);
    void InitVi(ConstParam, Particle *, double, double);
};

//calculate charge density with weightning
class SolveField
{
public:
  double dx;
  double A[2048][3];

  void InitField(ConstParam, Particle *, GridPoint *);
  void InitRho(ConstParam, GridPoint *, HalfGridPoint *);

  void SetSolver(ConstParam);

  void CalcPhi(ConstParam, GridPoint *);
  void CalcE(ConstParam, GridPoint *, HalfGridPoint *);
  void SolveMaxwel(ConstParam, GridPoint *, HalfGridPoint *);

  void CalcDensity_2(ConstParam, Particle *, GridPoint *);
  void CalcDensity_4(ConstParam, Particle *, GridPoint *);

private:
  double delta;
};

class SolveMotion
{
public:
  double dx, dt;

  void InitE(ConstParam, Particle *);

  void weighting_2(ConstParam, GridPoint *, Particle *);
  void weighting_4(ConstParam, GridPoint *, Particle *);

  void weighting_2vhxe(ConstParam, GridPoint *, Particle *);
  void weighting_2vhxi(ConstParam, GridPoint *, Particle *);
  void weighting_4vhxe(ConstParam, GridPoint *, Particle *);
  void weighting_4vhxi(ConstParam, GridPoint *, Particle *);

  void average_e(ConstParam, GridPoint *, Particle *);
  void average_i(ConstParam, GridPoint *, Particle *);

  void defineB(ConstParam, Particle *, int);

  void setV(ConstParam, Particle *);
  void accel(ConstParam, Particle *);
  void move(ConstParam, Particle *);

  void boundary_left_mirror(ConstParam, Particle *);
  void boundary_right_mirror(ConstParam, Particle *);

  void boundary_left_adsorp(ConstParam, Particle *, Particle *, Particle *, Particle *);
  void boundary_right_adsorp(ConstParam, Particle *, Particle *, Particle *, Particle *);

  void boundary_left_adsorp_high(ConstParam, Particle *, Particle *);
  void boundary_right_adsorp_high(ConstParam, Particle *, Particle *);

private:
  double v_hm[3], v_m[3], v_p[3], v_o[3], B[3], T[3], S[3], delta;

  int tmp, tmp_rand;
};

class OutPut
{
public:
  FILE *trace[100];
  void Initialize();
  void Grid_Info(ConstParam, GridPoint *, int);
  void Particle_Info(ConstParam,Particle *,Particle *,int);
  void trajectory(ConstParam, Particle *, int, int, int);
  void History_Info(int,ConstParam,Particle *,Particle *,GridPoint *,HalfGridPoint *);
private:
    int part_per_file, file_per_rank, file_start, file_end, part_per_rank, start, end;//, temp;
    int right_frag, left_frag;
    double energyE, energyI;
    //FILE *f;
};
