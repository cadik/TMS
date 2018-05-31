/* --------------------------------------------------------------------------- *
 * Particles.cpp: implementation of the Particle.   *
 * --------------------------------------------------------------------------- */

#include "particles.h"
#include <cmath>
#include <iomanip>
#include "color.h"

//#define LOG_COMPUTE
//#define LOG_FORCE
//#define LOG_SPRING
//#define LOG_INPUT

#define EPSILON 0.000000000001

void ParticlesManager::createSprings()
{
    int countParticles = particlesCount;
	for(int i = 0;i < countParticles; i++)
	{
		for(int j = 0; j < countParticles; j++)
		{
			if(i==j)
			{
				continue;
			}

			Particle* p1 = particles[i];
			Particle* p2 = particles[j];

			double distance = CalculateDistance(p1, p2);
			
			double QRange = maxDistance; // 373
			double GRange = 100.0;
			double normalization = (GRange / QRange);

			Spring spring;
			spring.length = normalization * distance;
            spring.curLength = spring.length;
            spring.outParticle = j;

        #ifdef LOG_SPRING
                std::cerr << color::bgreen << "#### P" << p1->id+1 << "-P" << spring.outParticle+1 <<  
            " distance: " <<distance <<" QRange: " <<maxDistance  << " length: " <<spring.length << color::reset << std::endl;
        #endif

			
			p1->addSpring(spring);
		}
    }
}

double ParticlesManager::computeForce(Particle* particle)
{
    double force = 0.0;
    int springsCount = particle->springs.size();
	for(int j = 0; j < springsCount; j++)
	{
        Spring s = particle->springs.at(j);
        double tempForce = 0.0;

        double Li = particle->gray;
		double Lj = particles[s.outParticle]->gray;

        double actualLength = abs(Li - Lj);

        double first = 1.0 * (1 - (s.length / actualLength));
        double second = (Lj - Li);

        tempForce = first * second;

        force += tempForce;
    }
    return force;
}

void ParticlesManager::compute(double maxTime)
{
    std::cerr << std::setprecision(5) << std::fixed;

    double time = 0.0;
    double dt = 0.001;

    int index = 0;
    while(time < maxTime)
    {
        
        time += dt;
        if(index >= particlesCount)
        {
            index = 0;
        }
            
        Particle* p = particles[index++];

        double force = abs(p->mass) < EPSILON || p->mass == INFINITY ? 0.0 : computeForce(p);
        
        double ratio = force / p->mass;
        ratio *= dt * dt;

        double tempGray = p->gray;
        p->gray = 2.0 * p->gray - p->lastGray + ratio;
        p->lastGray = tempGray;

        if(p->gray >= 100.0)
		{
			p->gray = 100.0;
            p->mass = INFINITY;
		}

        if(p->gray <= 0.0)
		{
			p->gray = 0.0;
            p->mass = INFINITY;
        }
        #ifdef LOG_COMPUTE
		    std::cerr << color::bblue << "P" 
            <<index+1 << " G:" << p->gray << " lastG:" << p->lastGray 
            << " ratio:" << ratio << " f:" << force << " m:" << p->mass<< std::endl;
        #endif
    }
    
    std::cerr << color::bblue << "Time: " << time << color::reset << std::endl;

    /*

    prev_pos = pos
    

    while (pos > 0)
        time += dt
        temp_pos = pos
        pos = pos * 2 - prev_pos + acc * dt * dt
        prev_pos = temp_pos
    end

    println(time)*/
}


/********************************************
 * ******************************************
 * *****************************************/
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void ParticlesManager::computeMaxDistance()
{
	for(int i = 0;i < particlesCount; i++)
	{
		for(int j = 0; j < particlesCount; j++)
		{
			if(i==j)
			{
				continue;
			}

			Particle* p1 = particles[i];
			Particle* p2 = particles[j];
			double dL = pow(p2->L - p1->L,2);
			double da = pow(p2->a - p1->a,2);
			double db = pow(p2->b - p1->b,2); 
 
			
			double distance = sqrt(dL + da + db);
			if(distance > maxDistance)
			{
				maxDistance = distance;
			}
		}
	}
}
double ParticlesManager::CalculateDistance(Particle* p1, Particle* p2)
{
	double dL = pow(p2->L - p1->L,2);
	double da = pow(p2->a - p1->a,2);
	double db = pow(p2->b - p1->b,2); 
	double distance = sqrt(dL + da + db);
	return distance;
}

double ParticlesManager::getSkFactor(Particle* p)
{
	int countParticles = particlesCount;

	double up = 0;
	double bottom = 0;

	for(int i = 0; i < countParticles; i++)
	{
		if(p->id == i)
		{
			continue;
		}

		Particle* pi = particles[i];

		double distance = CalculateDistance(pi,p);
		double w = 1.0 / distance * distance;
		bottom += w;

		up += w*((abs(p->gray - pi->gray))/ distance );
	}
	return up / bottom;
}
/* Initialize mass of particle */
void Particle::computeMass()
{
	mass = 1.0/sqrt(a*a + b*b);
}
/* Method create particles and save in pool of particles */
void ParticlesManager::initialize(cv::Mat labels,cv::Mat centers)
{
	particles = new Particle*[centers.rows];
	particlesCount = 0;

	for(int i=0; i < centers.rows;i++)
	{
		Particle* particle = new Particle;
		particle->id = i;

		particle->L = centers.at<float>(i, 0);
		particle->a = centers.at<float>(i, 1);
		particle->b = centers.at<float>(i, 2);


		particle->gray = particle->L;
		particle->lastGray = particle->gray;

        // Compute mass of each particle
		particle->computeMass();

		particles[i] = particle;
		particlesCount++;
	}
	for( int y = 0; y < rows; y++ )
	{
		for( int x = 0; x < cols; x++ )
		{ 
			int cluster_idx = labels.at<int>(y + x*rows,0);
			pool[y][x] = particles[cluster_idx];
		}
	}

    #ifdef LOG_INPUT

    std::cerr << "############ Input ###########" << std::endl;
	for(int index = 0;index < particlesCount; index++)
	{
			Particle* p = particles[index];
			std::cerr << color::bblue << "P" <<index+1 << " L:" << p->L << " a:" << p->a << " b:" << p->b<< color::reset << std::endl;
	}
    #endif
}
void ParticlesManager::makeDifferentG()
{
	bool working = true;
	int countParticles = particlesCount;
	while(working)
	{
		bool wasChanged = false;
		for(int i = 0;i < countParticles; i++)
		{
			for(int j = 0; j < countParticles; j++)
			{
				if(i==j)
				{
					continue;
				}

				Particle* p1 = particles[i];
				Particle* p2 = particles[j];
				if(abs(p2->gray - p1->gray) < 0.001)
				{
					double deltaGray;
					while(abs(deltaGray = fRand(-2,2)) < 0.1);

					p1->gray += deltaGray;
					p1->lastGray = p1->lastGray;


					wasChanged = true;
				}
			}
		}
		if(!wasChanged)
		{
			working = false;
		}
	}
}

void Particle::addSpring(Spring s)
{
	springs.push_back(s);
}

ParticlesManager::~ParticlesManager()
{
	for(int i = 0; i < rows; ++i) {
		delete [] pool[i];
	}
	delete [] pool;
	delete [] particles;
}

/* Copy data from pool to destination image */
void ParticlesManager::toDestination(double* outData)
{
	for( int y = 0; y < rows; y++ )
	{
		for( int x = 0; x < cols; x++ )
		{ 
			if(isChrominance)
			{
				*outData++ = 0;
			}
			else
			{
				float L =  pool[y][x]->gray;
				*outData++ = L;
			}
			if(isColor || isChrominance)
			{
				float a = pool[y][x]->a;
				*outData++ = a;
				float b = pool[y][x]->b;
				*outData++ = b;
			}
			else
			{
				*outData++ = 0; 
				*outData++ =  0;
			}
			
		}
	}
}

/* Return particle from pool at col and row */
Particle* ParticlesManager::getParticle(int col, int row)
{
    return pool[row][col];
}

/* Alocate required structures */
ParticlesManager::ParticlesManager(int rows, int cols, bool isColor, bool isChrominance)
{
	this->rows = rows;
	this->cols = cols;
	this->isColor = isColor;
	this->isChrominance = isChrominance;

	pool = new Particle**[rows];
	for(int i = 0; i < rows; i++)
	{
		pool[i] = new Particle*[cols];
	}
	maxDistance = 0;

}