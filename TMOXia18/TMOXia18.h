#include "TMO.h"

#include <tensorflow/core/public/session.h>
#include <tensorflow/core/protobuf/meta_graph.pb.h>


using namespace std;
using namespace tensorflow;

class TMOXia18 : public TMO
{
public:
	TMOXia18();
	virtual ~TMOXia18();
	virtual int Transform();

protected:
    /*
     * Run a color to grayscale conversion or grayscale to color prediction.
     */
    Status run(Tensor &image,
               vector<Tensor> &result,
               int direction,
               Session* session);

    /*
     * The function fills the tensor with data from the source (pSrc).
     */
    void fillTensor(Tensor *t,
                    int batch,
                    int depth,
                    int windowPositionY,
                    int windowPositionX);

    /*
     * The function fills the destination (pDst) with data from the tensor.
     * The tensor is obtained by prediction.
     */
    void fillOutput(Tensor *t,
                    int batchY,
                    int batchX,
                    int batchCountY,
                    int batchCountX);

    /*
     * Function handles the work with the tensor and prediction.
     */
    Status predict(Session* session);

    /*
     * The function loads the model for the prediction.
     */
    Session* loadModel();

    /*
     * Color to grayscale (value: "gray") or grayscale to color direction (value: "color").
     */
    TMOString direction;
};
