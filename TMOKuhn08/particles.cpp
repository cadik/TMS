/* --------------------------------------------------------------------------- *
 * Particles.cpp: implementation of the Particle.   *
 * --------------------------------------------------------------------------- */

#include "particles.h"
#include <cmath>
#include <iomanip>

void ParticlesPool::computeMaxDistance()
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
			//double distance = abs(p1->L - p2->L);
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


void ParticlesPool::makeDifferentG()
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

double ParticlesPool::getDistance(Particle* p1, Particle* p2)
{
	double dL = pow(p2->L - p1->L,2);
	double da = pow(p2->a - p1->a,2);
	double db = pow(p2->b - p1->b,2); 
	double distance = sqrt(dL + da + db);
	return distance;
}

void ParticlesPool::createSprings()
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

			double distance = getDistance(p1, p2);
			
			double QRange = maxDistance;//373.0; 
			double GRange = 100.0;
			double normalization = (GRange / QRange);

			Spring spring;
			spring.length = normalization * distance;

			spring.outParticle = j;
			p1->addSpring(spring);
		}
	}
}

double ParticlesPool::getSkFactor(Particle* p)
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

		double distance = getDistance(pi,p);
		double w = 1.0 / distance * distance;
		bottom += w;

		up += w*((abs(p->gray - pi->gray))/ distance );
	}
	return up / bottom;
}





double ParticlesPool::computeForce(Particle* particle)
{
	double force = 0;
	int springsCount = particle->springs.size();
	for(int j = 0; j < springsCount; j++)
	{
		Spring s = particle->springs.at(j);
		double tempForce = 0.0;
		double Li = particle->gray;
		double Lj = particles[s.outParticle]->gray;

		double curLength = abs(Li - Lj);

		double first = 1.0*(1 - (s.length/ curLength));
		double second = (Lj - Li);

		tempForce = first * second;
		
		/*std::cerr << "---- P" << particle->id+1 << "-P" << s.outParticle+1 <<  
		" lenght: " <<s.length << " curLength:" << curLength <<  
		" Lj: " << Lj << " Li:" << Li << std::endl;
		
		std::cerr << "#### P" << particle->id+1 << "-P" << s.outParticle+1 <<  
		" first: " <<first <<" second: " <<second << std::endl;*/
        if((Li + tempForce) > Lj)
        {
            tempForce -= tempForce - Lj; 
        }
        else if((Li - tempForce) < Lj)
        {
            tempForce += tempForce + Lj; 
        }

        std::cerr << "#### P" << particle->id+1 << "-P" << s.outParticle+1 <<  
		" tempForce: " <<tempForce << std::endl;
        
		force += tempForce;

	}
    return force;
}


void ParticlesPool::compute(int steps)
{
	std::cerr << std::setprecision(5) << std::fixed;;

	makeDifferentG();

	bool isAllZero = true;
    bool isZero = false;
    bool isHundred = false;
	for(int t = 0; t < steps; t++)
	{
		std::cerr << "**********************************************************************************************************:step" << t+1 << std::endl;
		int countParticles = particlesCount;
		for(int index = 0;index < countParticles; index++)
		{
			std::cerr << "#############################################################################################:index" << index+1 << std::endl;

			Particle* p = particles[index];
			double force = computeForce(p);
            double diff = force / p->mass;
			double newGray = p->gray * 2.0 - p->lastGray;

            newGray += diff;
			//p->newGray = (force/p->mass) + ((2.0*p->gray) - p->lastGray);
			std::cerr << "P" <<index+1 << " newG:" << newGray << " G:" << p->gray << " lastG:" << p->lastGray << " diff:" << diff << " f:" << force << " m:" << p->mass<< std::endl;

			if(newGray >= 100.0)
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
            }
			std::cerr << "P" <<index+1 << " newG:" << newGray << std::endl;


			p->lastGray = p->gray;
			p->gray = newGray;

		}
		
	}
}


void Particle::computeMass()
{
	mass = sqrt(a*a + b*b);
}








/* --------------------------------------------------------------------------- *
 * Base method for creating and deleting object          			           *
 * --------------------------------------------------------------------------- */
void ParticlesPool::initialize(cv::Mat labels,cv::Mat centers)
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

	//std::cerr << "############ Input ###########" << std::endl;
	for(int index = 0;index < particlesCount; index++)
	{
			Particle* p = particles[index];
			//std::cerr << "P" <<index+1 << " L:" << p->L << " a:" << p->a << " b:" << p->b<< std::endl;
	}
}

Spring::~Spring()
{
}
Spring::Spring()
{
	
}

void Particle::addSpring(Spring s)
{
	springs.push_back(s);
}

Particle::Particle()
{
}
Particle::~Particle()
{
}

ParticlesPool::~ParticlesPool()
{
	for(int i = 0; i < rows; ++i) {
		delete [] pool[i];
	}
	delete [] pool;
	delete [] particles;
}

void ParticlesPool::toDestination(double* outData)
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

ParticlesPool::ParticlesPool(int rows, int cols, bool isColor, bool isChrominance)
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