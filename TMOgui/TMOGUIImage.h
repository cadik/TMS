// TMOGUIImage.h: interface for the TMOGUIImage class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TMOGUIIMAGE_H
#define TMOGUIIMAGE_H
#include <QWidget>
//Added by qt3to4:
#include <QResizeEvent>
#include <QEvent>
#include <QLabel>
#include <string>
#include "TMOGUIBitmap.h"
#include "TMOGUIAdjustValues.h"
#include "TMOGUIResource.h"

//using namespace std;

class QPushButton;
class QScrollArea;
class QStatusBar;
class QLabel;
class TMOGUIProgressBar;
class TMOGUIOutput;
class TMOGUIAdjust;
class TMOGUITransformation;

class TMOGUIImage : public QWidget
{
	Q_OBJECT
public:
	static const wchar_t* GetString(const QChar* pChars);
	virtual bool IsModified() {return pDst;};
	virtual int SetImage(TMOImage *pDest);
    virtual int Swap(bool display);
	virtual TMOGUITransformation* Transform();
	virtual TMOImage* GetImage();
	virtual int Open(const char* filename);
	virtual int New(TMOGUIImage *image);
	virtual int Extract(TMOGUIImage *pSrcImage, int iComponent);
	virtual int MergeComponents(TMOGUIImage* pRed, TMOGUIImage* pGreen, TMOGUIImage* pBlue);
	virtual int ImageOperation(TMOGUIImage* pRed, TMOGUIImage* pGreen, int iOperation);
	virtual int New(int iWidth, int iHeight, double *pColors, int iPlacement);	
    virtual int NewSmall(TMOGUIImage* pSrcImage);
    virtual int Terminate();
    virtual void fitHisto();
	virtual void fitToScreen(QSize size);
	virtual void fitToWidth(QSize size);
	virtual void fitToHeight(QSize size);
	virtual void zoomIn();
	virtual void zoomOut();
	virtual bool CanUndo(void);
	virtual void SetImageZoomLabel();
	virtual void showtools();
    virtual void deleteDest();
    virtual void hideAll(bool);
	static QString GetName(QString filename);
    TMOGUIImage(TMOGUIProgressBar *pInitBar, QWidget * parent=0, const char * name=0, bool isPreview=0);
	virtual ~TMOGUIImage();
	TMOGUIAdjust* pAdjust;
	TMOGUIBitmap* pImage;
	TMOGUIOutput* pOutput;
	QPushButton* pToolsButton;
    QString* imageName;
    bool bPreview;
protected:
    virtual void customEvent( QEvent * e ); // QEvent
	virtual void resizeEvent ( QResizeEvent * );
	QSize size;
    QScrollArea* pScrollView;
	QLabel* pZoom;
	TMOGUIProgressBar* pProgress;
	TMOGUIProgressBar *pInitProgress;
    TMOGUIAdjustValues *values;
    TMOGUIAdjustValues *filters;
	TMOGUITransformation *pTransform;
	TMOImage* pSrc;
	TMOImage* pDst;
    //QHBox* pHBox;
    QWidget* pHBox;
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

