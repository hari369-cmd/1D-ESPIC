#include "define.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <random>

void SetParticle::SetMQ(ConstParam cp, Particle *p, double MASS, double CHARGE)
{
	for (int i = 0; i < cp.N; i++)
	{
		p[i].q = CHARGE * cp.NAP;
		p[i].m = MASS * cp.NAP;
	}
}

void SetParticle::SetPosition(ConstParam cp, Particle *p)
{
	std::random_device rnd;
	std::mt19937 mt(rnd());
	std::uniform_real_distribution<> rand(cp.margin * cp.DX, cp.L - cp.margin * cp.DX);
	for (int i = 0; i < cp.N; i++)
	{
		p[i].x = rand(mt);
		p[i].x_ini = p[i].x;
		p[i].over_frag = 0;
	}
}

void SetParticle::InitVe(ConstParam cp, Particle *p, double TEMP, double TEMP1, double M)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    double MASS;
    MASS = M;
    double dell = 1 / sqrt(MASS / (TEMP * -cp.QE));
    std::normal_distribution<> dist(0.0, dell);
    double dell1 = 1 / sqrt(MASS / (TEMP1 * -cp.QE));
    std::normal_distribution<> dist1(0.0, dell1);
    
    for (std::size_t n = 0; n < cp.N; n++)
    {
        if(n <= ceil(cp.alpha / 100 * cp.N))
        {
            p[n].vx = dist(engine);
            dist.reset();
        }
        else
        {
            p[n].vx = dist(engine);
            dist.reset();
        }
    }
    
    for (std::size_t n = 0; n < cp.N; n++)
    {
        p[n].vy = dist1(engine);
    }
    dist.reset();
    for (std::size_t n = 0; n < cp.N; n++)
    {
        p[n].vz = dist(engine);
    }
    dist.reset();
}

void SetParticle::InitVi(ConstParam cp, Particle *p, double TEMP, double M)
{
	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());
	double MASS;
	MASS = M;
	double dell = 1 / sqrt(MASS / (TEMP * -cp.QE));
	std::normal_distribution<> dist(0.0, dell);

	for (std::size_t n = 0; n < cp.N; n++)
	{
		p[n].vx = dist(engine);
	}
	dist.reset();
	for (std::size_t n = 0; n < cp.N; n++)
	{
		p[n].vy = dist(engine);
	}
	dist.reset();
	for (std::size_t n = 0; n < cp.N; n++)
	{
		p[n].vz = dist(engine);
	}
    dist.reset();
}
