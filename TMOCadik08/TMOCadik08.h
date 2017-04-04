//(c)Martin Cadik
//03/2007

#include <string>
#include "TMO.h"
#include "coloroid.h"
#include "compute/parallel.h"
#include "compute/command_queue.h"
#include "compute/program.h"

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

	private:
	Coloroid model;
	// OpenCL objects
	const cl::parallel simd; // OpenCL environment
	const cl::command_queue step;
	mutable cl::program exe;
	const size_t wgs, // OpenCL work group size
	             dim; // sqrt(wgs) -- for 2D kernel invocations

	double formulaColoroid(const double* const data,
	                       const long y1, const long x1,
	                       const long y2, const long x2,
	                       const long xmax);
	void correctGrad(TMOImage&, const double) const;
	cl::event reduce(const std::string, const cl::buffer&, const unsigned, double&,
	                 const cl::event_list = {}) const;
	void integrate2x(TMOImage&, TMOImage&) const;
	void transRange(TMOImage&, const double, const double) const;
	// CPU
	void inconsistencyCorrection(TMOImage& G_image,
	                             const double eps);
	void GFintegration(TMOImage& G_image, TMOImage& pDst);
	void calibrate(TMOImage& src_image, TMOImage& dst_image);
}; //TMOCadik08
