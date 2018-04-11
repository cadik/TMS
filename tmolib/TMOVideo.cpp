#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TMOImage.h"
#include "TMOVideo.h"
#include "TMORadiance.h"
//-OpenEXR-
#ifndef LINUX
 #define OPENEXR_DLL
#else
 #define _stricmp strcasecmp
#endif

///constrictor
TMOVideo::TMOVideo(const char *filename)
{
	
	OpenVideo(filename);
}
TMOVideo::TMOVideo()
{

}

TMOVideo::~TMOVideo()
{
  captureObject.release();
  writerObject.release();
}

/**
 *Sets the name of the output file, borrowed from TMOImage 
 **/
int TMOVideo::setNameOut(char* suffix)
{
  int length = std::strlen(vName);
  int suffixlength = std::strlen(suffix);
  int point;
  int i;
  char *filename = new char[length + suffixlength + 1 + 128];


  for (point = length - 1; point > 0; point--)
	  if (vName[point] == '.') break;
  if (!point) point = length; 
  for (i = 0; i < point; i++) filename[i] = vName[i];
  for (i = 0; i < suffixlength; i++) filename[point+i] = suffix[i];
  filename[point+i] = 0;
  
  strcat(filename, ".avi");
  vNameOut = filename;
  return 0;
}

/**
 * Creates the video writer object so that frames can be wrriten into a video
 * @ref TMOVideo object that contains the data from the source video
 * */
int TMOVideo::createOutputVideo(const TMOVideo &ref)
{
  

  cv::VideoWriter out = cv::VideoWriter(ref.vNameOut,CV_FOURCC('M','J','P','G'), ref.fps, cv::Size(ref.frameWidth,ref.frameHeight));
  if(!out.isOpened()) throw -1;
  writerObject = out;
  return 0;
}
/**
 * Creates video writer object from name
 * @filename name of the video to be created
 * @width width of the video
 * @height height of the video
 * */
int TMOVideo::createOutputVideoByName(const char *filename, int width, int height)
{
  cv::VideoWriter out = cv::VideoWriter(filename,CV_FOURCC('M','J','P','G'), 30, cv::Size(width,height));
  if(!out.isOpened()) throw -1;
  writerObject = out;
}
/**
 * Sets the frame to an output video, the frame is of type cv::Mat
 * @out video writer object of the destination video
 * @frame cv::Mat frame which is to be written
 * */
int TMOVideo::setMatFrame(cv::VideoWriter out, cv::Mat frame)
{
  cv::normalize(frame,frame,0,255,cv::NORM_MINMAX,CV_8UC3);
  out.write(frame);
  return 0;
}
/**
 * Sets the frame to an output video, the frame is of type TMOImage
 * @out video writer object of the destination video
 * @frame TMOImage frame which is to be written
 * */
int TMOVideo::setTMOImageFrame(cv::VideoWriter out, TMOImage &img)
{
   cv::Mat frame;
   frame.create(img.GetWidth(), img.GetHeight(), CV_32FC3);
  TMOImageToCvMat(img, frame);
  cv::normalize(frame,frame,0,255,cv::NORM_MINMAX,CV_8UC3);
  out.write(frame);

}

/**
 * Opens the desired video and sets neccesary video parameters
 * @filename name of the file to be isOpened
 **/
int TMOVideo::OpenVideo(const char *filename)
{
  int length = strlen(filename), i;
  int counter=0;
  cv::Mat tmp;

 // if (vName) delete[] vName;
  vName = new char[length + 6];
  strcpy(vName, filename);

  cv::VideoCapture cap;
  
  
  if(!cap.open(vName)) throw -1;
  
  for(int i=0;i<cap.get(CV_CAP_PROP_FRAME_COUNT );i++)
  {
    cap.read(tmp);
    if(tmp.rows == 0 || tmp.cols==0) break;   ///chcecking for number of frames
    counter++;
  }
  cap.set(CV_CAP_PROP_POS_FRAMES, 0); ///set ro beginning

  captureObject=cap;
  totalNumberOfFrames = counter;//cap.get(CV_CAP_PROP_FRAME_COUNT ); ///this get turned to be unreliable in wmv files
  frameHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT );
  frameWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH ); 
  fps = cap.get(CV_CAP_PROP_FPS);
  return 0;
}
/**
 * Conversion from TMOImage to cv::math
 * @img TMOImage to be converted
 * @frame cv::Mat in which the conversion will be stored
 * */
int TMOVideo::TMOImageToCvMat(TMOImage &img, cv::Mat &frame)
{
  int height =img.GetHeight();
  int width = img.GetWidth();
  cv::Mat mergedMat=cv::Mat::zeros( height,width, CV_64FC3);
 
  cv::Mat blueMat=cv::Mat::zeros( height,width, CV_64FC1);
  cv::Mat greenMat=cv::Mat::zeros( height,width, CV_64FC1);
  cv::Mat redMat=cv::Mat::zeros( height,width, CV_64FC1);
   
  std::vector<cv::Mat> channels;
  double* data = img.GetData();
  double p;
   for(int j=0; j<height;j++)
  {
    for(int i=0;i<width;i++)
    {
      redMat.at<double>(j,i)= *data++;
      greenMat.at<double>(j,i)= *data++;
      blueMat.at<double>(j,i)= *data++;
    }
  }

 channels.push_back(blueMat);
 channels.push_back(greenMat);
  channels.push_back(redMat); 
  cv::merge(channels,mergedMat);

 frame = mergedMat;

  blueMat.release();
  greenMat.release();
  redMat.release();
 // delete [] data;
  return 0;
}
/**
 * Conversion from cv::Mat to TMOImage
 * @img TMOImage where the conversion will be stored
 * @c video capture object of source video
 * @frame cv::Mat containing the data from source vieo
 * */
int TMOVideo::cvMatToTMOImage( TMOImage &img, cv::VideoCapture c, cv::Mat &frame)
{
 
  int width = c.get(CV_CAP_PROP_FRAME_WIDTH ); 
  int height = c.get(CV_CAP_PROP_FRAME_HEIGHT );
  double* frameData = new double[width*height*3];
  cv::Vec3d colors;
  

  int b=0;
  for(int j=0; j<height;j++)
  {
    for(int i=0;i<width;i++)
    {
      colors = frame.at<cv::Vec3f>(j,i);
     frameData[b++] =colors[2]; //red
     frameData[b++] = colors[1]; // green
     frameData[b++] = colors[0]; // blue
     
    }
  }
 
  img.New(width,height,TMO_RGB);
  
   
 img.SetData(frameData);

  
}
/**
 * Gets the desired video frame from the video
 * @c video capture object of the source video
 * @frameNumber number of the frame to be fetched
 * @img TMOImage which will be returned containg all the frame data
 * */
int TMOVideo::getTMOImageVideoFrame(cv::VideoCapture c, int frameNumber, TMOImage &img)
{
  
  
  cv::Mat frame;
   
  if (frameNumber < totalNumberOfFrames) 
  {
    c.set(CV_CAP_PROP_POS_FRAMES,frameNumber);
    c.read(frame);
    
     cv::normalize(frame,frame,0.0,1.0,cv::NORM_MINMAX,CV_32FC3);
    int dd=frame.type();
    TMOVideo::cvMatToTMOImage(img,c,frame);
  }
  else{
    std::cerr << "Error, frame number is bigger than the total number of frames." << std::endl;
	return -1;
  }
  return 0;
  
}
/**
 * Gets desired video frame from the video
 * @c video capture object of the source video
 * @frameNumber number of the frame to be fetched
 * @frame cv::mat which will be returned containg all the frame data
 * */
int TMOVideo::GetMatVideoFrame(cv::VideoCapture c, int frameNumber, cv::Mat &frame)
{
  
  if (frameNumber < totalNumberOfFrames) 
  {
    c.set(CV_CAP_PROP_POS_FRAMES,frameNumber);
    c.read(frame);
     cv::normalize(frame,frame,0.0,1.0,cv::NORM_MINMAX,CV_32FC3);
  
  }
  else{
    std::cerr << "Error, frame number is bigger than the total number of frames." << std::endl;
	return -1;
  }
  return 0;
    
}
int TMOVideo::setVName(char* name)
{
  vName = name;
  return 0;
}

