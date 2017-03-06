//(c)Martin Cadik
//03/2007

//Color2Gray in gradient domain - for CAE07 in Banff
//based on CATAM,
//http://www.maths.cam.ac.uk/undergrad/tripos/catam/

//#define ERIKLIB
#include "TMO.h"
#include "compute/platform.h"
#include "compute/device.h"
#include "compute/context.h"
#include "compute/buffer.h"
#include "compute/command_queue.h"
#include "compute/program.h"

class TMOColor2Gray : public TMO
{
	public:
	TMOColor2Gray();
	virtual ~TMOColor2Gray();
	virtual int Transform();
	void GFintegration(TMOImage& G_image, TMOImage& pDst);
	void calibrate(TMOImage& src_image, TMOImage& dst_image);

	protected:
	TMOBool normalize;			//urcuje zda bude vystupni luminance normalizovana do <0, 1>
	TMODouble gamma;			//gamma correction factor
	TMODouble s;				//overprojection step
	TMODouble eps;				//epsilon for gradient correction

	private:
	int xmax, ymax;

	const cl::platform host;
	const cl::device gpu;
	const cl::context env;
	const cl::command_queue step;
	cl::program exe;
	const size_t lms, // local mem. size
	             wgs; // work group size

	//double formula(double* data, long y1, long x1, long y2, long x2, double* thresholddata);	
	double formulaColoroid(const double* const data,
	                       const long y1, const long x1,
	                       const long y2, const long x2,
	                       double* const thresholdData);	
	void inconsistencyCorrection(TMOImage& G_image,
	                             TMOImage& divG_image,
	                             const double eps);
	void correct_grad(TMOImage&, const double);
	void integrate(TMOImage&, TMOImage&);
	cl::event reduce_max(const cl::buffer&, unsigned, double&,
	                     const std::vector<cl::event> = {});
}; //TMOColor2Gray
