#include "TMOExports.h"
#include <assert.h>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"
#include "opencv2/videoio/videoio.hpp"

#undef EPS

class TMOLIB_API TMOVideo
{
protected:
      int frameWidth;
      int frameHeight;
      cv::VideoCapture captureObject;
      cv::VideoWriter writerObject;
      char *vName;
      char *vNameOut;
      int fps;
      int fourcc; ///video codek
      int totalNumberOfFrames;

public:
      virtual int OpenVideo(const char *filename);
      virtual cv::VideoCapture getVideoCaptureObject() { return captureObject; }
      virtual cv::VideoWriter getVideoWriterObject() { return writerObject; }
      virtual int GetMatVideoFrame(cv::VideoCapture c, int frameNumber, cv::Mat &frame);
      virtual int createOutputVideo(const TMOVideo &ref);
      virtual int setMatFrame(cv::VideoWriter out, cv::Mat frame);
      virtual int setNameOut(char *suffix);
      virtual int cvMatToTMOImage(TMOImage &img, cv::VideoCapture c, cv::Mat &frame);
      virtual int getTMOImageVideoFrame(cv::VideoCapture c, int frameNumber, TMOImage &img);
      virtual int TMOImageToCvMat(TMOImage &img, cv::Mat &frame);
      virtual int setTMOImageFrame(cv::VideoWriter out, TMOImage &img);
      virtual int GetHeight() { return frameHeight; }
      virtual int GetWidth() { return frameWidth; }
      virtual int GetTotalNumberOfFrames() { return totalNumberOfFrames; }
      virtual int createOutputVideoByName(const char *filename, int width, int height);
      virtual int setVName(char *name);
      virtual char *getVNameOut() { return vNameOut; }

      TMOVideo(const char *filename);
      TMOVideo();
      ~TMOVideo();
};