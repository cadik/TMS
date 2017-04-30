// original authors: Laszlo Neumann,
//                   Martin Cadik

//(c)Martin Cadik
//03/2007

#define PROFILE

#include <string>
#include "TMO.h"
#include "coloroid.h"
#include "compute/parallel.h"
#include "compute/command_queue.h"
#include "compute/program.h"
#include "vec.h"

class quadtree;

class TMOCadik08 : public TMO {
	public:
	TMOCadik08();
	virtual ~TMOCadik08();
	virtual int Transform();

	protected:
	TMOBool normalize;	//urcuje zda bude vystupni luminance normalizovana do <0, 1>
	TMODouble gamma;	//gamma correction factor
	mutable TMODouble s;		//overprojection step
	TMODouble eps;		//epsilon for gradient correction
	TMOString type;

	private:
	Coloroid model;
	const cl::parallel simd; // OpenCL environment
	const cl::command_queue step;
	mutable cl::program exe;
	const size_t wgs, // OpenCL work group size
	             dim; // sqrt(wgs) -- for 2D kernel invocations
	const unsigned maba;

	double formulaColoroid(const double* const data,
	                       const long y1, const long x1,
	                       const long y2, const long x2,
	                       const long xmax);

	void correctGradCb(std::vector<vec2d>&, const unsigned,
	                   const unsigned, const double) const;
	void correctGradCbLoc(std::vector<vec2d>&, const unsigned,
                              const unsigned, const double) const;

	void correctGradCyc(std::vector<vec2d>&, const unsigned,
	                    const unsigned, const double) const;

	cl::event evalQuadtree(const cl::buffer&, const unsigned,
	                       const unsigned, vec2d* const, cl::event_list = {}) const;

	void correctGradHier(quadtree&, const double) const;

	cl::event reduce(const std::string, const cl::buffer&, const unsigned, double&,
	                 const cl::event_list = {}) const;

	cl::event scan(const std::string type, const cl::buffer& in,
	               const unsigned n, const cl::event_list pending) const;

	void integrate2x(TMOImage&, TMOImage&) const;

	void calibrate(TMOImage& src_image, TMOImage& dst_image);

	// CPU
	void inconsistencyCorrection(std::vector<vec2d>&, const long, const long,
	                             const double);

	void GFintegration(TMOImage& G_image, TMOImage& pDst);

#ifdef PROFILE
	mutable cl_ulong acc{};
	void accumulate(const cl::event&) const;
#endif
}; //TMOCadik08
