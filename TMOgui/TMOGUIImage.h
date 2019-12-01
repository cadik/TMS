// TMOGUIImage.h: interface for the TMOGUIImage class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TMOGUIIMAGE_H
#define TMOGUIIMAGE_H
#include <q3vbox.h>
//Added by qt3to4:
#include <QResizeEvent>
#include <QCustomEvent>
#include <QLabel>
#include <string>
#include "TMOGUIBitmap.h"
#include "TMOGUIAdjustValues.h"
#include "TMOGUIResource.h"

//using namespace std;

class QPushButton;
class Q3ScrollView;
class QStatusBar;
class QLabel;
class TMOGUIProgressBar;
class TMOGUIOutput;
class TMOGUIAdjust;
class TMOGUITransformation;

class TMOGUIImage : public Q3VBox  
{
	Q_OBJECT
public:
	static const wchar_t* GetString(const QChar* pChars);
	virtual bool IsModified() {return pDst;};
	virtual int SetImage(TMOImage *pDest);
	virtual int Swap();
	virtual TMOGUITransformation* Transform();
	virtual TMOImage* GetImage();
	virtual int Open(const char* filename);
	virtual int New(TMOGUIImage *image);
	virtual int Extract(TMOGUIImage *pSrcImage, int iComponent);
	virtual int MergeComponents(TMOGUIImage* pRed, TMOGUIImage* pGreen, TMOGUIImage* pBlue);
	virtual int ImageOperation(TMOGUIImage* pRed, TMOGUIImage* pGreen, int iOperation);
	virtual int New(int iWidth, int iHeight, double *pColors, int iPlacement);	
	virtual int Terminate();
	virtual void fitToScreen(QSize size);
	virtual void fitToWidth(QSize size);
	virtual void fitToHeight(QSize size);
	virtual void zoomIn();
	virtual void zoomOut();
	virtual bool CanUndo(void);
	virtual void SetImageZoomLabel();
	virtual void showtools();
	static QString GetName(QString filename);
	TMOGUIImage(TMOGUIProgressBar *pInitBar, QWidget * parent=0, const char * name=0);
	virtual ~TMOGUIImage();
	TMOGUIAdjust* pAdjust;
	TMOGUIBitmap* pImage;
	TMOGUIOutput* pOutput;
	QPushButton* pToolsButton;

protected:
	virtual void customEvent( QCustomEvent * e );
	virtual void resizeEvent ( QResizeEvent * );
	QSize size;
	Q3ScrollView* pScrollView;
	QLabel* pZoom;
	TMOGUIProgressBar* pProgress;
	TMOGUIProgressBar *pInitProgress;
	TMOGUIAdjustValues values;
	TMOGUIAdjustValues filters;
	TMOGUITransformation *pTransform;
	TMOImage* pSrc;
	TMOImage* pDst;
	Q3HBox* pHBox; 
	QStatusBar *pStatus;	
	QLabel* pTransformLabel;
	QWidget* pParent;
	bool bMaximized;
	bool bTools;
	bool bTransforming;
	int iCounter;

public slots:
	int SetImageSize(int iWidth, int iHeight);
	bool close();
	void setsize();	
	void canceltransform();
	
signals:
	void closeFile();
	void finishTransform();
};

#endif //TMOGUIIMAGE_H

