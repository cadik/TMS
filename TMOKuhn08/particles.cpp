/* --------------------------------------------------------------------------- *
 * Particles.cpp: implementation of the Particle.   *
 * --------------------------------------------------------------------------- */

#include "particles.h"
#include <cmath>
#include <iomanip>
#include "color.h"

#define LOG_COMPUTE
#define LOG_FORCE
//#define LOG_SPRING
//#define LOG_INPUT

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

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
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

double ParticlesManager::CalculateDistance(Particle* p1, Particle* p2)
{
	double dL = pow(p2->L - p1->L,2);
	double da = pow(p2->a - p1->a,2);
	double db = pow(p2->b - p1->b,2); 
	double distance = sqrt(dL + da + db);
	return distance;
}

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





void ParticlesManager::computeForce()
{
    for(int index = 0;index < particlesCount; index++)
	{
        Particle* particle = particles[index];

        double force = 0;
        int springsCount = particle->springs.size();
        for(int j = 0; j < springsCount; j++)
        {
            Spring s = particle->springs.at(j);
            double tempForce = 0.0;
            double Li = particle->gray;
            double Lj = particles[s.outParticle]->gray;

            double QRange = maxDistance; // 373
                double GRange = 100.0;
                double normalization = (GRange / QRange);

            s.curLength = s.length + abs(Li - Lj);

            double first = 1.0*(1 - (s.length/ s.curLength));
            double second = (Lj - Li);

            tempForce = first * second;
            
            #ifdef LOG_FORCE
                std::cerr << "---- P" << particle->id+1 << "-P" << s.outParticle+1 <<  
                " lenght: " <<s.length << " curLength:" << s.curLength <<  
                " Lj: " << Lj << " Li:" << Li << std::endl;
                
                std::cerr << "#### P" << particle->id+1 << "-P" << s.outParticle+1 <<  
                " first: " <<first <<" second: " <<second << std::endl;

                std::cerr << "#### P" << particle->id+1 << "-P" << s.outParticle+1 <<  
                " tempForce: " <<tempForce << std::endl;
            #endif

            /*if(Li < Lj && (Li + tempForce) > Lj)
            {
                tempForce -= tempForce - Lj; 
            }
            else if(Li > Lj && (Li - tempForce) < Lj)
            {
                tempForce += tempForce + Lj; 
            }*/
            #ifdef LOG_FORCE

                std::cerr << "#### P" << particle->id+1 << "-P" << s.outParticle+1 <<  
                " tempForce: " <<tempForce << std::endl;
            #endif

            force += tempForce;

        }
        particle->force = force;
    }
}


void ParticlesManager::compute(int steps)
{
	std::cerr << std::setprecision(5) << std::fixed;;

	makeDifferentG();

	bool isAllZero = true;
    bool isZero = false;
    bool isHundred = false;
	for(int t = 0; t < steps; t++)
	{
        #ifdef LOG_COMPUTE
		    std::cerr << color::red << "**********************************************************************************************************:step" << t+1 << color::reset << std::endl;
		#endif
        computeForce();
		for(int index = 0;index < particlesCount; index++)
		{
			Particle* p = particles[index];

            double diff = p->force / p->mass;
			double newGray = p->gray * 2.0 - p->lastGray;

            newGray += diff;
			//p->newGray = (force/p->mass) + ((2.0*p->gray) - p->lastGray);
            #ifdef LOG_COMPUTE
			    std::cerr << color::bblue << "P" 
                <<index+1 << " newG:" << newGray << " G:" << p->gray << " lastG:" << p->lastGray 
                << " diff:" << diff << " f:" << p->force << " m:" << p->mass<< std::endl;
            #endif
			/*if(newGray >= 100.0)
			{
                if(isHundred)
                {
                    newGray = p->gray;
                }
                else
                {
				    newGray = 100.0;
                    isHundred = true;
                }
			}

            if(newGray <= 0.0)
			{
                if(isZero)
                {
                    newGray = p->gray;
                }
                else
                {
				    newGray = 0.0;
                    isZero = true;
                }
            }*/
            #ifdef LOG_COMPUTE
			    std::cerr << color::bblue << "P" <<index+1 << " newG:" << newGray << color::reset << std::endl;
            #endif

			p->lastGray = p->gray;
			p->gray = newGray;

		}
		
	}
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