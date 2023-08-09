/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2023                                              *
*                                                                              *
*                       Implementation of the TMOXia18 class                   *
*                                                                              *
*******************************************************************************/

#include "TMOXia18.h"

#include <filesystem>

#define TMS_SOURCE_DEPTH 3
#define WINDOW_SIZE 256
#define WINDOW_CROP_ONE_SIZE 224
#define WINDOW_CROP_SIZE 192

typedef std::vector<std::pair<std::string, tensorflow::Tensor>> tensorDict;

using std::filesystem::current_path;
using std::filesystem::path;

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOXia18::TMOXia18() {
    SetName(L"Xia18");
    SetDescription(L"Invertible Grayscale");

    direction.SetName(L"c");
    direction.SetDescription(
            L"Convert to grayscale or to color. value: color,gray (default gray)");
    direction.SetDefault("gray");
    direction = "gray";
    this->Register(direction);
}

TMOXia18::~TMOXia18() {
}

/*
 * Run a color to grayscale conversion or grayscale to color prediction.
 */
Status TMOXia18::run(Tensor &image,
                     vector<Tensor> &result,
                     int direction,
                     Session* session) {
    if (direction == 0) {
        tensorDict feedDict2 = {
                {"input/batch", image}
        };

        return session->Run(feedDict2, {"encode/latent/Tanh"}, {}, &result);
    } else {
        tensorDict feedDict2 = {
                {"encode/latent/Tanh", image}
        };

        return session->Run(feedDict2, {"decode/Tanh"}, {}, &result);
    }
}

/*
 * The method fills the tensor with data from the source (pSrc).
 */
void TMOXia18::fillTensor(Tensor *t,
                          int batch,
                          int depth,
                          int windowPositionY,
                          int windowPositionX) {
    auto inputTensorMapped = t->tensor<float, 4>();

    int yEnd = std::min(this->iHeight, WINDOW_SIZE);
    int xEnd = std::min(this->iWidth, WINDOW_SIZE);

    for (int y = 0; y < yEnd; y++) {
        const double* sourceRow = (double*)this->pSrc->GetData() +
                (windowPositionY + y) * this->iWidth * TMS_SOURCE_DEPTH;

        for (int x = 0; x < xEnd; x++) {
            const double* sourcePixel = sourceRow + (windowPositionX + x) * TMS_SOURCE_DEPTH;

            for (int c = 0; c < depth; c++) {
                const double* sourceValue = sourcePixel + c;
                inputTensorMapped(batch, y, x, c) = (float) *sourceValue;
            }
        }
    }

    if (this->iHeight < WINDOW_SIZE) {
        for (int y = this->iHeight; y < WINDOW_SIZE; y++) {
            for (int x = 0; x < xEnd; x++) {
                for (int c = 0; c < depth; c++) {
                    inputTensorMapped(batch, y, x, c) = 0.0f;
                }
            }
        }
    }

    if (this->iWidth < WINDOW_SIZE) {
        for (int y = 0; y < yEnd; y++) {
            for (int x = this->iWidth; x < WINDOW_SIZE; x++) {
                for (int c = 0; c < depth; c++) {
                    inputTensorMapped(batch, y, x, c) = 0.0f;
                }
            }
        }
    }
}

/*
 * The method fills the destination (pDst) with data from the tensor.
 * The tensor is obtained by prediction.
 */
void TMOXia18::fillOutput(Tensor *t,
                          int batchY,
                          int batchX,
                          int batchCountY,
                          int batchCountX) {
    int yEnd = (batchY + 1 < batchCountY) ? batchY * WINDOW_CROP_SIZE + WINDOW_CROP_ONE_SIZE : this->iHeight;
    int xEnd = (batchX + 1 < batchCountX) ? batchX * WINDOW_CROP_SIZE + WINDOW_CROP_ONE_SIZE : this->iWidth;
    int yStart = (batchY > 0) ? (batchY - 1) * WINDOW_CROP_SIZE + WINDOW_CROP_ONE_SIZE : 0;
    int xStart = (batchX > 0) ? (batchX - 1) * WINDOW_CROP_SIZE + WINDOW_CROP_ONE_SIZE : 0;

    int windowTransferRangeY = (batchY + 1 < batchCountY) ? batchY * WINDOW_CROP_SIZE : yEnd - WINDOW_SIZE;
    int windowTransferRangeX = (batchX + 1 < batchCountX) ? batchX * WINDOW_CROP_SIZE : xEnd - WINDOW_SIZE;

    // Change to 0 if the value is negative.
    windowTransferRangeY = std::max(windowTransferRangeY, 0);
    windowTransferRangeX = std::max(windowTransferRangeX, 0);

    auto outputTensorMapped = t->tensor<float, 4>();

    for (int y = yStart; y < yEnd; y++) {
        for (int x = xStart; x < xEnd; x++) {
            for (int c = 0; c < TMS_SOURCE_DEPTH; c++) {
                this->pDst->GetData()[this->iWidth * y * TMS_SOURCE_DEPTH + x * TMS_SOURCE_DEPTH + c] =
                    (double) outputTensorMapped(batchY * batchCountX + batchX,
                                                  y - windowTransferRangeY,
                                                  x - windowTransferRangeX, c);
            }
        }
    }
}

/*
 * The method handles the work with the tensor and prediction.
 */
Status TMOXia18::predict(Session* session) {
    int depth = (direction.GetString().compare("color") == 0) ? 1 : 3;

    int batchCountX = static_cast<int>(ceil((this->iWidth - WINDOW_SIZE) / ((double)WINDOW_CROP_SIZE))) + 1;
    int batchCountY = static_cast<int>(ceil((this->iHeight - WINDOW_SIZE) / ((double)WINDOW_CROP_SIZE))) + 1;

    // Change to 0 if the value is negative.
    batchCountX = std::max(batchCountX, 1);
    batchCountY = std::max(batchCountY, 1);

    Tensor *t = new Tensor(DT_FLOAT, {batchCountX * batchCountY, WINDOW_SIZE, WINDOW_SIZE, depth});

    int windowPositionY = 0;
    int windowPositionX;

    for (int y = 0; y < batchCountY; y++) {
        windowPositionX = 0;

        for (int x = 0; x < batchCountX; x++) {
            fillTensor(t, y * batchCountX + x, depth, windowPositionY, windowPositionX);

            if (x + 1 < batchCountX - 1) {
                windowPositionX += WINDOW_CROP_SIZE;
            } else {
                windowPositionX = this->iWidth - WINDOW_SIZE;
            }
        }

        if (y + 1 < batchCountY - 1) {
            windowPositionY += WINDOW_CROP_SIZE;
        } else {
            windowPositionY = this->iHeight - WINDOW_SIZE;
        }
    }

    int directionInt = (direction.GetString().compare("color") == 0) ? 1 : 0;

    std::vector<Tensor> outputs;
    auto status = run(*t, outputs, directionInt, session);

    if (status.ok()) {
        for (int y = 0; y < batchCountY; y++) {
            for (int x = 0; x < batchCountX; x++) {
                fillOutput(&(outputs[0]), y, x, batchCountY, batchCountX);
            }
        }
    }

    return status;
}

/*
 * The method loads the model for the prediction.
 *
 * The method is inspired from the following source:
 * Source: https://github.com/tensorflow/tensorflow/issues/32624#issuecomment-537285073
 * Author: amitpandey2194 (https://github.com/amitpandey2194)
 */
Session* TMOXia18::loadModel() {
    // Set up input paths.
    path tmsPath = current_path().parent_path() / "TMOXia18" / "model";

    const string pathToGraph =  (tmsPath / "model.cpkt-202080.meta").string();
    const string checkpointPath =  (tmsPath / "model.cpkt-202080").string();

    // Prepare session.
    auto session = NewSession(SessionOptions());

    if (session == nullptr) {
        throw runtime_error("Could not create Tensorflow session.");
    }

    Status status;

    // Read in the protobuf graph were exported.
    MetaGraphDef graphDef;
    status = ReadBinaryProto(Env::Default(), pathToGraph, &graphDef);
    if (!status.ok()) {
        throw runtime_error("Error reading graph definition from " + pathToGraph + ": " + status.ToString());
    }

    // Add the graph to the session.
    status = session->Create(graphDef.graph_def());
    if (!status.ok()) {
        throw runtime_error("Error creating graph: " + status.ToString());
    }

    // Read weights from the saved checkpoint.
    const std::string restoreOpName = graphDef.saver_def().restore_op_name();
    const std::string filenameTensorName =
            graphDef.saver_def().filename_tensor_name();

    Tensor filenameTensor(DT_STRING,
                          TensorShape());
    filenameTensor.scalar<tstring>()() = checkpointPath;

    tensorDict feedDict = {{filenameTensorName, filenameTensor}};
    status = session->Run(feedDict, {}, {restoreOpName}, nullptr);

    if (!status.ok()) {
        throw runtime_error("Error loading checkpoint from " + checkpointPath + ": " + status.ToString());
    }

    return session;
}

/*
 * Main method for TMOXia18 operator.
 */
int TMOXia18::Transform() {
    if (direction.GetString().compare("color") != 0 &&
        direction.GetString().compare("gray") != 0) {

        // Prediction failed.
        throw runtime_error("Error: Invalid option value (color or gray): " + direction.GetString());
    }

    Session* session = loadModel();
    auto status = predict(session);

    if (!status.ok()) {
        // Prediction failed.
        throw runtime_error("Error: Prediction failed: " + status.ToString());
    }

    return 0;
}
