#include "opencv2/opencv.hpp"

class Spring
{
public:
	Spring();
	virtual ~Spring();

	double length;
	int outParticle;
protected:
};

class Particle
{
public:
	Particle();
	virtual ~Particle();

	double mass;
	double L;
	double a;
	double b;

	double gray;
	double lastGray;
	int id;

	void addSpring(Spring s);
	std::vector<Spring> springs;
	void computeMass();

protected:
};


class ParticlesPool
{
public:
	int rows, cols;
	float maxDistance;
	ParticlesPool(int rows, int cols, bool isColor = false, bool isChrominance = false);
	virtual ~ParticlesPool();

	void computeMaxDistance();
	double computeForce(Particle* particle);
	double getDistance(Particle* p1, Particle* p2);

	void initialize(cv::Mat labels,cv::Mat centers);
	void createSprings();
	void compute(int steps);
	void makeDifferentG();

	void toDestination(double* outData);
	double getSkFactor(Particle* p);
	bool isColor;
	bool isChrominance;
	int particlesCount;
	Particle** particles;
	Particle*** pool;
};
