#include "opencv2/opencv.hpp"

class Spring
{
public:
	Spring() {}
	virtual ~Spring() {}
	double length;
	double curLength;
	int outParticle;

protected:
};

class Particle
{
public:
	Particle() {}
	virtual ~Particle() {}

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
};

class ParticlesManager
{
public:
	int rows, cols;
	float maxDistance;
	ParticlesManager(int rows, int cols, bool isColor = false, bool isChrominance = false);
	virtual ~ParticlesManager();

	void computeMaxDistance();
	double computeForce(Particle *particle);
	double CalculateDistance(Particle *p1, Particle *p2);

	void initialize(cv::Mat labels, cv::Mat centers);
	void createSprings();
	void compute(double maxTime);
	void makeDifferentG();

	void toDestination(double *outData);
	double getSkFactor(Particle *p);
	bool isColor;
	bool isChrominance;
	int particlesCount;
	Particle *getParticle(int col, int row);

private:
	Particle **particles;
	Particle ***pool;
};
