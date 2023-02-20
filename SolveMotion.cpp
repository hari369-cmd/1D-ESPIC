#include "define.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <random>

void SolveMotion::InitE(ConstParam cp, Particle *p)
{

    for (int i = 0; i < cp.N; i++)
    {
        p[i].E = 0;
    }
}

void SolveMotion::defineB(ConstParam cp, Particle *p, int t)
{
    for (int i = 0; i < cp.N; i++)
    {
        if (cp.BSS <= t)
        {
            p[i].Bx = cp.B_MN; //Initial disturbance, if any, in the x-direction. In this case, zero
            p[i].By = 0;
            p[i].Bz = cp.BI * exp(-pow(-0.1 + p[i].x - (t - cp.BSS) * cp.DT * cp.BV, 2) * (2 / (cp.BW * cp.BW)));
            // A moving Gaussian pulse in Bz
        }
        else
        {
            p[i].Bx = cp.B_MN;
            p[i].By = 0;
            p[i].Bz = 0;
        }
    }
}


void SolveMotion::setV(ConstParam cp, Particle *p)
{
    dt = cp.DT;

    for (int i = 0; i < cp.N; i++)
    {
        B[0] = p[i].Bx;
        B[1] = p[i].By;
        B[2] = p[i].Bz;

        T[0] = p[i].q * -dt / (4 * p[i].m) * B[0];
        T[1] = p[i].q * -dt / (4 * p[i].m) * B[1];
        T[2] = p[i].q * -dt / (4 * p[i].m) * B[2];

        S[0] = 2 * T[0] / (1 + (T[0] * T[0] + T[1] * T[1] + T[2] * T[2]));
        S[1] = 2 * T[1] / (1 + (T[0] * T[0] + T[1] * T[1] + T[2] * T[2]));
        S[2] = 2 * T[2] / (1 + (T[0] * T[0] + T[1] * T[1] + T[2] * T[2]));

        v_hm[0] = p[i].vx;
        v_hm[1] = p[i].vy;
        v_hm[2] = p[i].vz;

        v_m[0] = v_hm[0] + p[i].q / p[i].m * p[i].E * dt * -0.25;
        v_m[1] = v_hm[1];
        v_m[2] = v_hm[2];

        v_o[0] = v_m[0] + (v_m[1] * T[2] - v_m[2] * T[1]);
        v_o[1] = v_m[1] + (v_m[2] * T[0] - v_m[0] * T[2]);
        v_o[2] = v_m[2] + (v_m[0] * T[1] - v_m[1] * T[0]);

        v_p[0] = v_m[0] + (v_o[1] * S[2] - v_o[2] * S[1]);
        v_p[1] = v_m[1] + (v_o[2] * S[0] - v_o[0] * S[2]);
        v_p[2] = v_m[2] + (v_o[0] * S[1] - v_o[1] * S[0]);

        p[i].vhx = v_p[0] + p[i].q / p[i].m * p[i].E * dt * -0.25;
        p[i].vhy = v_p[1];
        p[i].vhz = v_p[2];

    }
}

void SolveMotion::accel(ConstParam cp, Particle *p)
{
    dt = cp.DT;

    for (int i = 0; i < cp.N; i++)
    {
        B[0] = p[i].Bx;
        B[1] = p[i].By;
        B[2] = p[i].Bz;

        T[0] = p[i].q * dt / (2 * p[i].m) * B[0];
        T[1] = p[i].q * dt / (2 * p[i].m) * B[1];
        T[2] = p[i].q * dt / (2 * p[i].m) * B[2];

        S[0] = 2 * T[0] / (1 + (T[0] * T[0] + T[1] * T[1] + T[2] * T[2]));
        S[1] = 2 * T[1] / (1 + (T[0] * T[0] + T[1] * T[1] + T[2] * T[2]));
        S[2] = 2 * T[2] / (1 + (T[0] * T[0] + T[1] * T[1] + T[2] * T[2]));

        v_hm[0] = p[i].vhx;
        v_hm[1] = p[i].vhy;
        v_hm[2] = p[i].vhz;

        v_m[0] = v_hm[0] + p[i].q / p[i].m * p[i].E * dt * 0.5;
        v_m[1] = v_hm[1];
        v_m[2] = v_hm[2];

        v_o[0] = v_m[0] + (v_m[1] * T[2] - v_m[2] * T[1]);
        v_o[1] = v_m[1] + (v_m[2] * T[0] - v_m[0] * T[2]);
        v_o[2] = v_m[2] + (v_m[0] * T[1] - v_m[1] * T[0]);

        v_p[0] = v_m[0] + (v_o[1] * S[2] - v_o[2] * S[1]);
        v_p[1] = v_m[1] + (v_o[2] * S[0] - v_o[0] * S[2]);
        v_p[2] = v_m[2] + (v_o[0] * S[1] - v_o[1] * S[0]);

        p[i].vhx = v_p[0] + p[i].q / p[i].m * p[i].E * dt * 0.5;
        p[i].vhy = v_p[1];
        p[i].vhz = v_p[2];
    }
}

void SolveMotion::move(ConstParam cp, Particle *p)
{
    for (int i = 0; i < cp.N; i++)
    {
        if (p[i].over_frag == 0)
        {
            p[i].x = p[i].x + p[i].vhx * cp.DT;
        }
    }
}

void SolveMotion::boundary_left_mirror(ConstParam cp, Particle *p) //Mainly for the 1st particle which will be at x=-1/2 due to the algorithm initialization
{
    for (int i = 0; i < cp.N; i++)
    {
        if (p[i].x < 0)
        {
            p[i].x = -p[i].x;
            p[i].vhx = -p[i].vhx;
        }
    }
}
void SolveMotion::boundary_right_mirror(ConstParam cp, Particle *p)
{
    for (int i = 0; i < cp.N; i++)
    {
        if (p[i].x > cp.L)
        {
            p[i].x = cp.L - (p[i].x - cp.L);
            p[i].vhx = -p[i].vhx;
        }
    }
}

void SolveMotion::boundary_left_adsorp_high(ConstParam cp, Particle *e, Particle *Ion)
{
    double left_wall = cp.DX * cp.margin;
    //Absorption BC for particles
    for (int i = 0; i < cp.N; i++)
    {
        if (e[i].x <= left_wall)
        {
            e[i].x = left_wall;

            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;

            e[i].over_frag = -1;
        }
        else if (e[i].over_frag == -1)
        {
            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;
        }
    }

    for (int i = 0; i < cp.N; i++) // Finds a corresponding electron which is similarily out of bounds
    {
        if (Ion[i].x <= left_wall)
        {
            for (int j = 0; j < cp.N; j++)
            {
                if (e[j].over_frag == -1)
                {
                    tmp = j;
                    break;
                }
            }

            e[tmp].over_frag = 0;

            e[tmp].x = e[(tmp + i) % cp.N].x_ini;
            e[tmp].vhx = e[tmp].vx;
            e[tmp].vhy = e[tmp].vy;
            e[tmp].vhz = e[tmp].vz;

            Ion[i].x = e[tmp].x;
            Ion[i].vhx = Ion[i].vx;
            Ion[i].vhy = Ion[i].vy;
            Ion[i].vhz = Ion[i].vz;
        }
    }
}

void SolveMotion::boundary_right_adsorp_high(ConstParam cp, Particle *e, Particle *Ion)
{
    double right_wall = cp.L - (cp.DX * cp.margin);

    for (int i = 0; i < cp.N; i++)
    {
        if (e[i].x >= (cp.L - (cp.DX * 200)))
        {
            e[i].bound_frag = 1;
        }
        else
        {
            e[i].bound_frag = 0;
        }

        if (Ion[i].x >= (cp.L - (cp.DX * 200)))
        {
            Ion[i].bound_frag = 1;
        }
        else
        {
            Ion[i].bound_frag = 0;
        }
    }

    for (int i = 0; i < cp.N; i++)
    {
        if (e[i].x > right_wall)
        {
            e[i].x = right_wall;

            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;

            e[i].over_frag = 1;
        }
        else if (e[i].over_frag == 1)
        {
            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;
        }
    }

    for (int i = 0; i <= cp.N; i++)
    {
        if (Ion[i].x > right_wall)
        {
            for (int j = 0; j < cp.N; j++)
            {
                if (e[j].over_frag == 1)
                {
                    tmp = j;
                    break;
                }
            }
            e[tmp].over_frag = 0;

            //Decide injection rate
            e[tmp].x = e[(tmp + i) % cp.N].x_ini;
            e[tmp].vhx = e[tmp].vx;
            e[tmp].vhy = e[tmp].vy;
            e[tmp].vhz = e[tmp].vz;

            Ion[i].x = e[tmp].x;
            Ion[i].vhx = Ion[i].vx;
            Ion[i].vhy = Ion[i].vy;
            Ion[i].vhz = Ion[i].vz;

            //printf("right %e %e %e %e\n", e[tmp].x, e[tmp].vhx, Ion[i].x, Ion[i].vhx);
        }
    }
}

void SolveMotion::boundary_left_adsorp(ConstParam cp, Particle *e, Particle *Ion, Particle *e_dammy, Particle *I_dammy)
{
    double left_wall = cp.DX * cp.margin;
    for (int i = 0; i < cp.N; i++)
    {
        if (e[i].x < left_wall)
        {
            e[i].x = left_wall;

            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;

            e[i].over_frag = -1;
        }
        else if (e[i].over_frag == -1)
        {
            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;
        }
    }

    for (int i = 0; i < cp.N; i++)
    {
        if (Ion[i].x <= left_wall)
        {
            for (int j = 0; j < cp.N; j++)
            {
                if (e[j].over_frag == -1)
                {
                    tmp = j;
                    break;
                }
            }

            e[tmp].over_frag = 0;
            tmp_rand = (tmp + i) % cp.N;

            e[tmp].x = e_dammy[tmp_rand].x;
            e[tmp].vhx = -2E+6;
            e[tmp].vhy = 0;
            e[tmp].vhz = 0;

            Ion[i].x = e[tmp].x;
            Ion[i].vhx = I_dammy[tmp_rand].vx;
            Ion[i].vhy = 0;
            Ion[i].vhz = 0;

            // printf("left %e %e %e %e\n", e[tmp].x, e[tmp].vhx, Ion[i].x, Ion[i].vhx);
        }
    }
}

void SolveMotion::boundary_right_adsorp(ConstParam cp, Particle *e, Particle *Ion, Particle *e_dammy, Particle *I_dammy)
{
    double right_wall = cp.L - cp.DX * cp.margin;
    for (int i = 0; i < cp.N; i++)
    {
        if (e[i].x > right_wall)
        {
            e[i].x = right_wall;

            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;

            e[i].over_frag = 1;
        }
        else if (e[i].over_frag == 1)
        {
            e[i].vhx = 0;
            e[i].vhy = 0;
            e[i].vhz = 0;
        }
    }

    for (int i = 0; i < cp.N; i++)
    {
        if (Ion[i].x >= right_wall)
        {
            for (int j = 0; j < cp.N; j++)
            {
                if (e[j].over_frag == 1)
                {
                    tmp = j;
                    break;
                }
            }
            e[tmp].over_frag = 0;

            tmp_rand = (tmp + i) % cp.N;

            e[tmp].x = e_dammy[tmp_rand].x;

            e[tmp].vhx = 2E+6;
            e[tmp].vhy = 0;
            e[tmp].vhz = 0;

            Ion[i].x = e[tmp].x;

            Ion[i].x = e[tmp].x;
            Ion[i].vhx = I_dammy[tmp_rand].vx;
            Ion[i].vhy = 0;
            Ion[i].vhz = 0;

            //printf("right %e %e %e %e\n", e[tmp].x, e[tmp].vhx, Ion[i].x, Ion[i].vhx);
        }
    }
}

void SolveMotion::weighting_2(ConstParam cp, GridPoint *gp, Particle *p)
{
    dx = cp.DX;
    delta = 0;
    double j = 0.0;

    for (int i = 0; i < cp.N; i++)
    {
        if(p[i].near == 0)
        {
            delta = fabs(p[i].x - gp[0].X) / dx;
            p[i].E += gp[0].E * (0.75 - pow(delta, 2));
            p[i].E += gp[1].E * 0.5 * pow(0.5 + delta, 2);
        }
        else if (p[i].near == cp.NG)
        {
            delta = fabs(gp[cp.NG].X - p[i].x) / dx;
            p[i].E += gp[cp.NG - 1].E * 0.5 * pow(0.5 - delta, 2);
            p[i].E += gp[cp.NG].E * (0.75 - pow(delta, 2));
        }
        else
        {
            delta = fabs(gp[p[i].near].X - p[i].x) / dx;
            p[i].E += gp[p[i].near - 1].E * 0.5 * pow(0.5 - delta, 2);
            p[i].E += gp[p[i].near].E * (0.75 - pow(delta, 2));
            p[i].E += gp[p[i].near + 1].E * 0.5 * pow(0.5 + delta, 2);
        }
    }
}

void SolveMotion::weighting_2vhxe(ConstParam cp, GridPoint *gp, Particle *p)
{
    dx = cp.L / cp.NG;
    delta = 0;
    double j = 0.0;

    // Calculate the weighting electron velocity at grid points. Will be used for finding
    // the average electron velocity at grid points
    for (int i = 0; i < cp.N; i++)
    {
        j = p[i].vhx;
        if(p[i].near != 0 && p[i].near != cp.NG)
        {
            gp[p[i].near].den_e += cp.NAP;
            gp[p[i].near - 1].Vavg1_e += p[i].vhx * 0.5 * pow(0.5 - delta, 2);
            gp[p[i].near].Vavg1_e += p[i].vhx * (0.75 - pow(delta, 2));
            gp[p[i].near + 1].Vavg1_e += p[i].vhx * 0.5 * pow(0.5 + delta, 2);
        }

        else if (p[i].near == 0)
        {
            gp[p[i].near].den_e += cp.NAP;
            gp[p[i].near].Vavg1_e += p[i].vhx * (0.75 - pow(delta, 2));
            gp[p[i].near + 1].Vavg1_e += p[i].vhx * 0.5 * pow(0.5 + delta, 2);
        }

        else
        {
            gp[p[i].near].den_e += cp.NAP;
            gp[p[i].near - 1].Vavg1_e += p[i].vhx * 0.5 * pow(0.5 - delta, 2);
            gp[p[i].near].Vavg1_e += p[i].vhx * (0.75 - pow(delta, 2));
        }

        p[i].vhx = j;
    }
}

void SolveMotion::weighting_2vhxi(ConstParam cp, GridPoint *gp, Particle *p)
{
    dx = cp.L / cp.NG;
    delta = 0;
    double j = 0.0;

    // Calculate the weighting ion velocity at grid points. Will be used for finding
    // the average ion velocity at grid points
    for (int i = 0; i < cp.N; i++)
    {
        j = p[i].vhx;
        if(p[i].near != 0 && p[i].near != cp.NG)
        {
            gp[p[i].near].den_i += 1;
            gp[p[i].near - 1].Vavg1_i += p[i].vhx * 0.5 * pow(0.5 - delta, 2);
            gp[p[i].near].Vavg1_i += p[i].vhx * (0.75 - pow(delta, 2));
            gp[p[i].near + 1].Vavg1_i += p[i].vhx * 0.5 * pow(0.5 + delta, 2);
        }

        else if (p[i].near == 0)
        {
            gp[p[i].near].den_i += 1;
            gp[p[i].near].Vavg1_i += p[i].vhx * (0.75 - pow(delta, 2));
            gp[p[i].near + 1].Vavg1_i += p[i].vhx * 0.5 * pow(0.5 + delta, 2);
        }

        else
        {
            gp[p[i].near].den_i += 1;
            gp[p[i].near - 1].Vavg1_i += p[i].vhx * 0.5 * pow(0.5 - delta, 2);
            gp[p[i].near].Vavg1_i += p[i].vhx * (0.75 - pow(delta, 2));
        }

        p[i].vhx = j;
    }
}

void SolveMotion::average_e(ConstParam cp, GridPoint *gp, Particle *p)
{
    // Find the averaged electron velocity at grid points
    for (int i = 0; i <= cp.NG; i++)
    {
        if (gp[i].den_e == 0)
        {
            gp[i].den_e = 1;
        }
        else
        {
            gp[i].Vavg_e = gp[i].Vavg1_e / gp[i].den_e;
            gp[i].Vavg1_e = 0;
            gp[i].den_e = 0;
        }
    }

    // Electron kinetic energy calculation
    for (int i = 0; i < cp.N; i++)
    {
        p[i].kee = 0.5 * cp.ME * pow((p[i].vhx - gp[p[i].near].Vavg_e),2);
        p[i].kex = 0.5 * cp.ME * p[i].vx * p[i].vx;
    }
}

void SolveMotion::average_i(ConstParam cp, GridPoint *gp, Particle *p)
{
    // Find the averaged ion velocity at grid points
    for (int i = 0; i <= cp.NG; i++)
    {
        if (gp[i].den_i == 0)
        {
            gp[i].den_i = 1;
        }
        else
        {
            gp[i].Vavg_i = gp[i].Vavg1_i / gp[i].den_i;
            gp[i].Vavg1_i = 0;
            gp[i].den_i = 0;
        }
    }

    // Ion kinetic energy calculation
    for (int i = 0; i < cp.N; i++)
    {
        p[i].kei = 0.5 * cp.MI * pow((p[i].vhx - gp[p[i].near].Vavg_i),2);
        p[i].kix = 0.5 * cp.MI * p[i].vx * p[i].vx;
    }
}

void SolveMotion::weighting_4(ConstParam cp, GridPoint *gp, Particle *p)
{
    dx = cp.L / cp.NG;
    delta = 0;
    for (int i = 0; i < cp.N; i++)
    {
        if (p[i].near == 0)
        {
            delta = fabs(p[i].x - gp[0].X) / dx;
            p[i].E += gp[0].E * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            p[i].E += gp[1].E * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            p[i].E += gp[2].E * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        else if (p[i].near == 1)
        {
            delta = fabs(p[i].x - gp[1].X) / dx;
            p[i].E += gp[0].E * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            p[i].E += gp[1].E * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            p[i].E += gp[2].E * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            p[i].E += gp[3].E * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        else if (p[i].near == cp.NG)
        {
            delta = fabs(gp[cp.NG].X - p[i].x) / dx;
            p[i].E += gp[cp.NG - 2].E * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            p[i].E += gp[cp.NG - 1].E * (0.375 + 0.009 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            p[i].E += gp[cp.NG].E * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
        }
        else if (p[i].near == cp.NG - 1)
        {
            delta = fabs(gp[cp.NG - 1].X - p[i].x) / dx;
            p[i].E += gp[cp.NG - 3].E * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            p[i].E += gp[cp.NG - 2].E * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            p[i].E += gp[cp.NG - 1].E * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            p[i].E += gp[cp.NG].E * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
        }
        else
        {
            delta = fabs(gp[p[i].near].X - p[i].x);
            p[i].E += gp[p[i].near - 2].E * (0.375 + 0.009 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            p[i].E += gp[p[i].near - 1].E * (0.375 + 0.009 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            p[i].E += gp[p[i].near].E * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            p[i].E += gp[p[i].near + 1].E * (0.375 + 0.009 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            p[i].E += gp[p[i].near + 2].E * (0.375 + 0.009 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
    }
}

void SolveMotion::weighting_4vhxe(ConstParam cp, GridPoint *gp, Particle *p)
{
    dx = cp.L / cp.NG;
    delta = 0;
    double a;
    for (int i = 0; i < cp.N; i++)
    {
        a = p[i].vhx;

        if (p[i].near == 0)
        {
            delta = fabs(p[i].x - gp[0].X) / dx;
            gp[0].den_e += cp.NAP;
            gp[p[i].near].Vavg_e += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near + 1].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            gp[p[i].near + 2].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        else if (p[i].near == 1)
        {
            delta = fabs(p[i].x - gp[1].X) / dx;
            gp[1].den_e += cp.NAP;
            gp[p[i].near - 1].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near].Vavg_e += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near + 1].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            gp[p[i].near + 2].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        else if (p[i].near == cp.NG)
        {
            delta = fabs(gp[cp.NG].X - p[i].x) / dx;
            gp[cp.NG].den_e += cp.NAP;
            gp[p[i].near - 2].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            gp[p[i].near - 1].Vavg_e += p[i].vhx * (0.375 + 0.009 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near].Vavg_e += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
        }
        else if (p[i].near == cp.NG - 1)
        {
            delta = fabs(gp[cp.NG - 1].X - p[i].x) / dx;
            gp[cp.NG - 1].den_e += cp.NAP;
            gp[p[i].near - 3].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            gp[p[i].near - 2].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near - 1].Vavg_e += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near].Vavg_e += p[i].vhx * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
        }
        else
        {
            delta = fabs(gp[p[i].near].X - p[i].x);
            gp[p[i].near].den_e += cp.NAP;
            gp[p[i].near - 2].Vavg_e += p[i].vhx * (0.375 + 0.009 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            gp[p[i].near - 1].Vavg_e += p[i].vhx * (0.375 + 0.009 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near].Vavg_e += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near + 1].Vavg_e += p[i].vhx * (0.375 + 0.009 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            gp[p[i].near + 2].Vavg_e += p[i].vhx * (0.375 + 0.009 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        p[i].vhx = a;
    }
}

void SolveMotion::weighting_4vhxi(ConstParam cp, GridPoint *gp, Particle *p)
{
    dx = cp.L / cp.NG;
    delta = 0;
    double a;
    for (int i = 0; i < cp.N; i++)
    {
        a = p[i].vhx;

        if (p[i].near == 0)
        {
            delta = fabs(p[i].x - gp[0].X) / dx;
            gp[0].den_i += cp.NAP;
            gp[p[i].near].Vavg_i += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near + 1].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            gp[p[i].near + 2].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        else if (p[i].near == 1)
        {
            delta = fabs(p[i].x - gp[1].X) / dx;
            gp[1].den_i += cp.NAP;
            gp[p[i].near - 1].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near].Vavg_i += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near + 1].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            gp[p[i].near + 2].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        else if (p[i].near == cp.NG)
        {
            delta = fabs(gp[cp.NG].X - p[i].x) / dx;
            gp[cp.NG].den_i += cp.NAP;
            gp[p[i].near - 2].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            gp[p[i].near - 1].Vavg_i += p[i].vhx * (0.375 + 0.009 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near].Vavg_i += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
        }
        else if (p[i].near == cp.NG - 1)
        {
            delta = fabs(gp[cp.NG - 1].X - p[i].x) / dx;
            gp[cp.NG - 1].den_i += cp.NAP;
            gp[p[i].near - 3].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            gp[p[i].near - 2].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near - 1].Vavg_i += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near].Vavg_i += p[i].vhx * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
        }
        else
        {
            delta = fabs(gp[p[i].near].X - p[i].x);
            gp[p[i].near].den_i += cp.NAP;
            gp[p[i].near - 2].Vavg_i += p[i].vhx * (0.375 + 0.009 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2));
            gp[p[i].near - 1].Vavg_i += p[i].vhx * (0.375 + 0.009 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2));
            gp[p[i].near].Vavg_i += p[i].vhx * (0.375 + 0.0096 * delta - 0.12 * delta * delta);
            gp[p[i].near + 1].Vavg_i += p[i].vhx * (0.375 + 0.009 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2));
            gp[p[i].near + 2].Vavg_i += p[i].vhx * (0.375 + 0.009 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2));
        }
        p[i].vhx = a;
    }
}
