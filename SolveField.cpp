#include "define.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>

void SolveField::InitField(ConstParam cp, Particle *p, GridPoint *gp)
{

	gp[0].X = 0;
	for (int i = 1; i <= cp.NG; i++)
	{
		gp[i].X = gp[i - 1].X + cp.DX;
	}
}

void SolveField::InitRho(ConstParam cp, GridPoint *gp, HalfGridPoint *gph)
{

	for (int i = 0; i < cp.NG + 1; i++)
	{
		gp[i].Rho = 0;
		gp[i].Phi = 0;
		gp[i].E = 0;
        gp[i].Vavg_e = 0;
        gp[i].Vavg_i = 0;
		gp[i].denseE = 0;
		gp[i].denseI = 0;
	}
	for (int i = 0; i < cp.NG; i++)
	{
		gph[i].E = 0;
	}
}

void SolveField::SetSolver(ConstParam cp)
{
	A[0][0] = 0;
	A[0][1] = -2;
	A[0][2] = 1;
	A[cp.NG - 2][0] = 1;
	A[cp.NG - 2][1] = -2;
	A[cp.NG - 2][2] = 0;
	for (int i = 1; i < cp.NG - 2; i++)
	{
		A[i][0] = 1;
		A[i][1] = -2;
		A[i][2] = 1;
	}
};
//Potential calculation by Thomas method
void SolveField::CalcPhi(ConstParam cp, GridPoint *gp)
{
	double e[cp.NG - 1], f[cp.NG - 1], a[cp.NG - 1], b[cp.NG - 1], c[cp.NG - 1];
	for (int i = 0; i < cp.NG - 1; i++)
	{
		a[i] = -A[i][0];
		b[i] = A[i][1];
		c[i] = -A[i][2];
	}
	e[0] = 0;
	f[0] = 0;

	for (int i = 1; i < cp.NG; i++)
	{
		e[i] = c[i-1] / (b[i-1] - a[i-1] * e[i - 1]);
		f[i] = (-gp[i].Rho * cp.DX * cp.DX / cp.EPSI + a[i-1] * f[i - 1]) / (b[i-1] - a[i-1] * e[i - 1]);
	}
	gp[0].Phi = 0;
	gp[cp.NG].Phi = 0;

	gp[cp.NG - 1].Phi = f[cp.NG - 1];
	for (int i = cp.NG - 2; i >= 1; i--)
	{
		gp[i].Phi = e[i] * gp[i + 1].Phi + f[i];
	}
}
//Electric field calculation from Electric Potential
void SolveField::CalcE(ConstParam cp, GridPoint *gp, HalfGridPoint *gph)
{

	for (int i = 0; i < cp.NG + 1; i++)
	{
		gph[i].E = -(gp[i + 1].Phi - gp[i].Phi) / cp.DX;
	}

	gp[0].E = gph[0].E;
	for (int i = 1; i < cp.NG; i++)
	{
		gp[i].E = (gph[i - 1].E + gph[i].E) * 0.5;
	}
	gp[cp.NG].E = gph[cp.NG - 1].E;
}

//Solving Poisson equation to find the Electric field
void SolveField::SolveMaxwel(ConstParam cp, GridPoint *gp, HalfGridPoint *gph)
{
	gp[0].E = 0;
	for (int i = 1; i <= cp.NG; i++)
	{
		gp[i].E = gp[i - 1].E + gph[i - 1].Rho * cp.DX / cp.EPSI;
	}
	for (int i = 0; i < cp.NG; i++)
	{
		gph[i].E = (gp[i].E + gp[i + 1].E) * 0.5;
	}
}

void SolveField::CalcDensity_2(ConstParam cp, Particle *p, GridPoint *gp)
{
	double v = cp.DX * cp.DS;
	delta = 0;

	//This for loop provides the values of p[i].near which will lie in-between 0 and NG in a jumbled order.
	//Further reference: Near_code(Number_Density)
    for (int i = 0; i < cp.N; i++)
	{
		if (fmod(p[i].x, cp.DX) <= cp.DX * 0.5)
		{
			p[i].near = p[i].x / cp.DX;
		}
		else
		{
			p[i].near = p[i].x / cp.DX + 1;
		}
	}

	//2nd order
	for (int j = 0; j < cp.N; j++)
	{
		if (p[j].q < 0)
		{
			gp[p[j].near].denseE += cp.NAP / v;
		}
		else
		{
			gp[p[j].near].denseI += cp.NAP / v;
		}

        if (p[j].near == 0)
		{
			delta = fabs(p[j].x - gp[0].X) / cp.DX;
			gp[0].Rho += p[j].q * (0.75 - delta * delta) / v;
			gp[1].Rho += p[j].q * 0.5 * (0.5 + delta) * (0.5 + delta) / v;
            if (p[j].q < 0)
            {
                gp[0].elect = gp[0].Rho;
                gp[1].elect = gp[1].Rho;
            }
            else
            {
                gp[0].ion = gp[0].Rho;
                gp[1].ion = gp[1].Rho;
            }
		}

		else if (p[j].near == cp.NG)
		{
			delta = fabs(gp[cp.NG].X - p[j].x) / cp.DX;
			gp[cp.NG - 1].Rho += p[j].q * 0.5 * (0.5 - delta) * (0.5 - delta) / v;
			gp[cp.NG].Rho += p[j].q * (0.75 - delta * delta) / v;
            if (p[j].q < 0) //Seperating charge densities into (+) and (-) by assigning them to gp.ion and gp.elect correspondingly
            {
                gp[cp.NG].elect = gp[cp.NG].Rho;
                gp[cp.NG - 1].elect = gp[cp.NG - 1].Rho;
            }
            else // Ion
            {
                gp[cp.NG].ion = gp[cp.NG].Rho;
                gp[cp.NG - 1].ion = gp[cp.NG - 1].Rho;
            }
		}

		else
		{
			delta = fabs(p[j].x - gp[p[j].near].X) / cp.DX;
			gp[p[j].near - 1].Rho += p[j].q * 0.5 * (0.5 - delta) * (0.5 - delta) / v;
			gp[p[j].near].Rho += p[j].q * (0.75 - delta * delta) / v;
			gp[p[j].near + 1].Rho += p[j].q * 0.5 * (0.5 + delta) * (0.5 + delta) / v;
            if (p[j].q < 0) // Electron
            {
                gp[p[j].near - 1].elect = gp[p[j].near - 1].Rho;
                gp[p[j].near].elect = gp[p[j].near].Rho;
                gp[p[j].near + 1].elect = gp[p[j].near + 1].Rho;
            }
            else // Ion
            {
                gp[p[j].near - 1].ion = gp[p[j].near - 1].Rho;
                gp[p[j].near].ion = gp[p[j].near].Rho;
                gp[p[j].near + 1].ion = gp[p[j].near + 1].Rho;
            }
		}
	}
}

void SolveField::CalcDensity_4(ConstParam cp, Particle *p, GridPoint *gp)
{
    double v = cp.DX * cp.DS;
    delta = 0;

    for (int i = 0; i < cp.N; i++)
    {
        if (fmod(p[i].x, cp.DX) <= cp.DX * 0.5)
        {
            p[i].near = p[i].x / cp.DX;
        }
        else
        {
            p[i].near = p[i].x / cp.DX + 1;
        }
    }

    //4th order
    for (int j = 0; j < cp.N; j++)
    {
        if (p[j].q < 0)
        {
            gp[p[j].near].denseE += cp.NAP / v;
        }
        else
        {
            gp[p[j].near].denseI += cp.NAP / v;
        }

        if (p[j].near == 0)
        {
            delta = fabs(p[j].x - gp[0].X) / cp.DX;
            gp[0].Rho += p[j].q * (0.375 + 0.0096 * delta - 0.12 * delta * delta) / v;
            gp[1].Rho += p[j].q * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2)) / v;
            gp[2].Rho += p[j].q * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2)) / v;
            if (p[j].q < 0) // Electron
            {
                gp[0].elect = gp[0].Rho;
                gp[1].elect = gp[1].Rho;
                gp[2].elect = gp[2].Rho;
            }
            else //Ion
            {
                gp[0].ion = gp[0].Rho;
                gp[1].ion = gp[1].Rho;
                gp[2].ion = gp[2].Rho;
            }
        }

        else if (p[j].near == 1)
        {
            delta = fabs(p[j].x - gp[1].X) / cp.DX;
            gp[0].Rho += p[j].q * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2)) / v;
            gp[1].Rho += p[j].q * (0.375 + 0.0096 * delta - 0.12 * delta * delta) / v;
            gp[2].Rho += p[j].q * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2)) / v;
            gp[3].Rho += p[j].q * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2)) / v;
            if (p[j].q < 0) // Electron
            {
                gp[0].elect = gp[0].Rho;
                gp[1].elect = gp[1].Rho;
                gp[2].elect = gp[2].Rho;
                gp[3].elect = gp[3].Rho;
            }
            else //Ion
            {
                gp[0].ion = gp[0].Rho;
                gp[1].ion = gp[1].Rho;
                gp[2].ion = gp[2].Rho;
                gp[3].ion = gp[3].Rho;
            }
        }

        else if (p[j].near == cp.NG)
        {
            delta = fabs(gp[cp.NG].X - p[j].x) / cp.DX;
            gp[cp.NG].Rho += p[j].q * (0.375 + 0.0096 * delta - 0.12 * delta * delta) / v;
            gp[cp.NG - 1].Rho += p[j].q * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2)) / v;
            gp[cp.NG - 2].Rho += p[j].q * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2)) / v;
            if (p[j].q < 0) // Electron
            {
                gp[p[j].near].elect = gp[p[j].near].Rho;
                gp[p[j].near - 1].elect = gp[p[j].near - 1].Rho;
                gp[p[j].near - 2].elect = gp[p[j].near - 2].Rho;
            }
            else //Ion
            {
                gp[p[j].near].ion = gp[p[j].near].Rho;
                gp[p[j].near - 1].ion = gp[p[j].near - 1].Rho;
                gp[p[j].near - 2].ion = gp[p[j].near - 2].Rho;
            }
        }

        else if (p[j].near == cp.NG - 1)
        {
            delta = fabs(gp[cp.NG - 1].X - p[j].x) / cp.DX;
            gp[cp.NG].Rho += p[j].q * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2)) / v;
            gp[cp.NG - 1].Rho += p[j].q * (0.375 + 0.0096 * delta - 0.12 * delta * delta) / v;
            gp[cp.NG - 2].Rho += p[j].q * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2)) / v;
            gp[cp.NG - 3].Rho += p[j].q * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2)) / v;
             if (p[j].q < 0) // Electron
            {
                gp[p[j].near].elect = gp[p[j].near].Rho;
                gp[p[j].near - 1].elect = gp[p[j].near - 1].Rho;
                gp[p[j].near - 2].elect = gp[p[j].near - 2].Rho;
                gp[p[j].near - 3].elect = gp[p[j].near - 3].Rho;
            }
            else //Ion
            {
                gp[p[j].near].ion = gp[p[j].near].Rho;
                gp[p[j].near - 1].ion = gp[p[j].near - 1].Rho;
                gp[p[j].near - 2].ion = gp[p[j].near - 2].Rho;
                gp[p[j].near - 3].ion = gp[p[j].near - 3].Rho;
            }
        }

        else
        {

            delta = fabs(p[j].x - gp[p[j].near].X) / cp.DX;
            gp[p[j].near + 2].Rho += p[j].q * (0.375 + 0.0096 * pow(delta - 2, 4) - 0.12 * pow(delta - 2, 2)) / v;
            gp[p[j].near + 1].Rho += p[j].q * (0.375 + 0.0096 * pow(delta - 1, 4) - 0.12 * pow(delta - 1, 2)) / v;
            gp[p[j].near].Rho += p[j].q * (0.375 + 0.0096 * delta - 0.12 * delta * delta) / v;
            gp[p[j].near - 1].Rho += p[j].q * (0.375 + 0.0096 * pow(delta + 1, 4) - 0.12 * pow(delta + 1, 2)) / v;
            gp[p[j].near - 2].Rho += p[j].q * (0.375 + 0.0096 * pow(delta + 2, 4) - 0.12 * pow(delta + 2, 2)) / v;
            if (p[j].q < 0) // Electron
            {
                gp[p[j].near + 2].elect = gp[p[j].near + 2].Rho;
                gp[p[j].near + 1].elect = gp[p[j].near + 1].Rho;
                gp[p[j].near].elect = gp[p[j].near].Rho;
                gp[p[j].near - 1].elect = gp[p[j].near - 1].Rho;
                gp[p[j].near - 2].elect = gp[p[j].near - 2].Rho;
            }
            else //Ion
            {
                gp[p[j].near + 2].ion = gp[p[j].near + 2].Rho;
                gp[p[j].near + 1].ion = gp[p[j].near + 1].Rho;
                gp[p[j].near].ion = gp[p[j].near].Rho;
                gp[p[j].near - 1].ion = gp[p[j].near - 1].Rho;
                gp[p[j].near - 2].ion = gp[p[j].near - 2].Rho;
            }
        }
    }
}
