// TMOGUIImage.cpp: implementation of the TMOGUIImage class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIImage.h"
#include "TMOGUIProgressBar.h"
#include "TMOGUIOutput.h"
#include "TMOGUIToneSlider.h"
#include "TMOGUIAdjust.h"
#include "TMOGUIInfo.h"
#include "TMOGUIInfoTool.h"
#include "TMOGUITransformation.h"
#include "TMOGUIHisto.h"
#include "../tmolib/TMO.h"
#include <qstatusbar.h>
#include <qapplication.h>
#include <qpushbutton.h>
#include <qscrollview.h>
#include <qlabel.h>
#include <qtooltip.h>
#include <qfont.h>
#include <qcursor.h>
#include <qmessagebox.h>
#include <math.h>
#include "TMOGUIImage.h"

//#include <iostream>

#define ZOOMSTEP 10
#define MIN_WIDTH 124
#define MIN_HEIGHT 136 // 220 X 101
#define ADJUST_MIN_EXT 36
#define MIN_ADJUSTWIDTH 340
#define MIN_ADJUSTHEIGHT 112
#define MIN_PROGRESSWIDTH 224
wchar_t buffer[256];

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIImage::TMOGUIImage(TMOGUIProgressBar *pInitBar, QWidget* parent, const char * name):
	QVBox(parent, name)
{
	pInitProgress = pInitBar;
	pImage = 0;
	bTools = false;
	pSrc = new TMOImage;
	pDst = 0;
	pParent = parent;
	bMaximized = false;
	bTransforming = false;
	iCounter = 0;
	
	setIcon(*TMOResource::pResource->IconMain->pixmap());

	// Scrollview
	pScrollView = new QScrollView(this, "Scrollview",  WRepaintNoErase|WResizeNoErase|WNorthWestGravity);

	// Image
	pImage = new TMOGUIBitmap(pScrollView->viewport(), "Bitmap");
	pScrollView->addChild(pImage);

	// Adjustment switch
	pToolsButton = new QPushButton("//...\\\\", this, "ToolsButton");
	QFont newFont = pToolsButton->font();
	newFont.setBold(true);
	pToolsButton->setFont(newFont);
	pToolsButton->setFixedHeight(9);
	pToolsButton->setToggleButton(true);
	pToolsButton->setFlat(false);
	pToolsButton->setCursor(QCursor::PointingHandCursor);
	QToolTip::add(pToolsButton, "Open Histogram");
	pToolsButton->hide();

	// Adjustments
	pAdjust = new TMOGUIAdjust(this, "Adjustments");
	pAdjust->hide();

	// Statusbar
	pStatus = new QStatusBar(this, "Statusbar");
	pHBox = new QHBox(pStatus, "Statushbox");
	pTransformLabel = new QLabel("Mapping...", pHBox, "TMOflag");
	pTransformLabel->hide();
	pZoom = new QLabel("100%", pHBox, "Zoomlabel");
	pZoom->hide();
	
	pStatus->addWidget(pHBox, 0, true);

	// Transform
	pTransform = new TMOGUITransformation(this);

	// Progressbar
	pProgress = new TMOGUIProgressBar(pStatus, "Progress");
	pProgress->SetLabel("");
	pTransform->Assign(pSrc);
	pSrc->SetProgress(pTransform->ProgressBar);
	connect(pProgress, SIGNAL(cancelsignal()), this, SLOT(canceltransform()));
	
	// Output
	pOutput = new TMOGUIOutput(TMOGUIOutput::pInfo, "Output");
	pOutput->Assign(pSrc);
	pSrc->SetWriteLine(pTransform->WriteLine);
	
	setCaption(GetName(name));
	setMinimumSize(MIN_PROGRESSWIDTH, MIN_HEIGHT);
	size = parent->size();
	//show();
	
	connect (&values, SIGNAL(render()), pImage, SLOT(valueschanged()));
	
	pImage->update();
	qApp->processEvents();
}

TMOGUIImage::~TMOGUIImage()
{
	if (pTransform) delete pTransform;
	if (pSrc) delete pSrc;
	if (pImage) delete pImage;
	if (pOutput) 
	{
		delete pOutput;
		pOutput = 0;
	}
}


QString TMOGUIImage::GetName(QString filename)
{
	QString s = filename;
	int pos = s.length()-1;

	while (pos)
	{
		if (s[pos--]=='/') break;
	}
	if (pos) s = s.mid(pos+2);

	return s;
}

bool TMOGUIImage::close()
{
	setFocus();
	emit closeFile();
	return true;
}

int TMOGUIImage::Open(const char *filename)
{
	QString s;
	pSrc->ProgressBar(0, 100);
	s = QString("Loading ") + GetName(QString(filename));
	pInitProgress->SetLabel(s);
	//adjustSize();
	try
	{
		pSrc->Open(filename);
	}
	catch ( int e )
	{
		pSrc->ProgressBar(100, 100);
		if (e == -19) 
		{
			return 1;
		}
		QString emsg = "";
		if( e == -2)emsg = "wrong file format\n";

		QMessageBox::critical( 0, "Tone Mapping Studio",
			QString("Failed to load file : \n\n") + name() + "\n" + emsg);// Dialog appearing
		return 2;
	}
	s = QString("File loaded from : ") + filename;
	try
	{
		pSrc->WriteLine(GetString(s.unicode()));
	}
	catch (...)
	{
		// Cannot write to info
	}
	pInitProgress->SetLabel("Converting");
	try
	{
		pSrc->Convert(TMO_RGB);
	}
	catch (int e)
	{
		if (e == -19) 
		{
			return 1;
		}
		QMessageBox::critical( 0, "Tone Mapping Studio",
			QString("Cannot convert file : \n\n") + name() + "\n");// Dialog appearing
		return 2;
	}
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Computing");
	pAdjust->Create(pSrc, &values);
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Displaying");
	pImage->Create(pSrc,&values,&filters,size,pProgress);
	pSrc->ProgressBar(100, 100);
	pAdjust->pToneSlider->setToWidth();
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	setsize();
	pToolsButton->show();
	pZoom->show();
	pInitProgress->SetProgress(100,100);
	pInitProgress = 0;
	show();
	return 0;	
}

int TMOGUIImage::New(TMOGUIImage* pSrcImage)
{
	QString s;
	s = QString("Creating ") + GetName(name());
	pInitProgress->SetLabel(s);
	pSrc->ProgressBar(0, 100);
	adjustSize();
	try
	{
		pSrc->New(*pSrcImage->GetImage(), TMO_RGB, true);
	}
	catch ( int e )
	{
		if (e == -19) 
		{
			pSrc->ProgressBar(0, 0);
			return 1;
		}
		QMessageBox::critical( 0, "Tone Mapping Studio",
			QString("Failed to create file file : \n\n") + name() + "\n");// Dialog appearing
		return 2;
	}
	s = QString("File duplicated from : ") + pSrcImage->name();
	pSrc->WriteLine(GetString(s.unicode()));

	pInitProgress->SetLabel("Computing");
	pAdjust->Create(pSrc, &values);
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Displaying");
	pImage->Create(pSrc,&values,&filters,size,pProgress);
	pSrc->ProgressBar(100, 100);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	setsize();
	pToolsButton->show();
	pZoom->show();
	
	pSrc->SetProgress(pTransform->ProgressBar);
	pInitProgress = 0;
	show();
	return 0;	
}

int TMOGUIImage::New(int iWidth, int iHeight, double *pColors, int iPlacement)
{
	QString s;
	s = QString("Creating ") + GetName(name());
	pInitProgress->SetLabel(s);
	pSrc->ProgressBar(0, 100);
	adjustSize();
	try
	{
		pSrc->New(iWidth, iHeight);
	}
	catch ( int e )
	{
		if (e == -19) 
		{
			// Dialog appearing
			pSrc->ProgressBar(0, 0);
			return 1;
		}
		QMessageBox::critical( 0, "Tone Mapping Studio",
			QString("Failed to load file : \n\n") + name() + "\n");// Dialog appearing
		return 2;
	}

	int iCenterX, iCenterY, iLength;

	switch(iPlacement)
	{
	case 0:
		iLength = sqrt((double)iWidth * iWidth + iHeight * iHeight);
		iCenterX = 0;
		iCenterY = 0;
		break;
	case 1:
		iLength = sqrt((double)iWidth * iWidth * 0.25 + iHeight * iHeight);
		iCenterX = iWidth >> 1;
		iCenterY = 0;
		break;
	case 2:
		iLength = sqrt((double)iWidth * iWidth + iHeight * iHeight);
		iCenterX = iWidth - 1;
		iCenterY = 0;
		break;
	case 3:
		iLength = sqrt((double)iWidth * iWidth + iHeight * iHeight * 0.25);
		iCenterX = 0;
		iCenterY = iHeight >> 1;
		break;
	case 4:
		iLength = sqrt((double)iWidth * iWidth * 0.25 + iHeight * iHeight * 0.25);
		iCenterX = iWidth >> 1;
		iCenterY = iHeight >> 1;
		break;
	case 5:
		iLength = sqrt((double)iWidth * iWidth + iHeight * iHeight * 0.25);
		iCenterX = iWidth - 1;
		iCenterY = iHeight >> 1;
		break;
	case 6:
		iLength = sqrt((double)iWidth * iWidth + iHeight * iHeight);
		iCenterX = 0;
		iCenterY = iHeight - 1;
		break;
	case 7:
		iLength = sqrt((double)iWidth * iWidth * 0.25 + iHeight * iHeight);
		iCenterX = iWidth >> 1;
		iCenterY = iHeight - 1;
		break;
	case 8:
		iLength = sqrt((double)iWidth * iWidth + iHeight * iHeight);
		iCenterX = iWidth - 1;
		iCenterY = iHeight - 1;
		break;
	}

	double dAspect, *pPixel;
	for (int i = 0; i < pSrc->GetHeight(); i++)
		for (int j = 0; j < pSrc->GetWidth(); j++)
		{
			pPixel = pSrc->GetPixel(j, i);
			if (iPlacement > 8)
			{
				pPixel[0] = pColors[0];
				pPixel[1] = pColors[1];
				pPixel[2] = pColors[2];
			}
			else
			{
				dAspect = sqrt((double)(iCenterX - j)*(iCenterX - j) + (iCenterY - i)*(iCenterY - i)) / iLength;
				pPixel[0] = dAspect * pColors[0] + (1 - dAspect) * pColors[3];
				pPixel[1] = dAspect * pColors[1] + (1 - dAspect) * pColors[4];
				pPixel[2] = dAspect * pColors[2] + (1 - dAspect) * pColors[5];
			}
		}

	s = QString("File created.");
	pSrc->WriteLine(GetString(s.unicode()));
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Computing");
	pAdjust->Create(pSrc, &values);
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Displaying");
	pImage->Create(pSrc,&values,&filters,size,pProgress);
	pSrc->ProgressBar(100, 100);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	setsize();
	pToolsButton->show();
	pZoom->show();
	
	pInitProgress = 0;
	show();
	return 0;	
}

void TMOGUIImage::SetImageZoomLabel()
{
	QString s;
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
}

void TMOGUIImage::fitToScreen(QSize size)
{
	QString s;
	double ratio = 1.0;
	double w = pImage->pSrcPixmap->width();
	double h = pImage->pSrcPixmap->height();
	double wWorkspace = size.width() - 10;
	double hWorkspace = size.height() - 63;
	double difW = (w - wWorkspace) / w;
	double difH = (h - hWorkspace) / h;
	if (difW >= difH)
	{
		if ((w - wWorkspace) >= 0)
			ratio = 1.0 - difW;
	}
	else
	{
		if ((h - hWorkspace) >= 0)
			ratio = 1.0 - difH;
	}
	pImage->Zoom(ratio);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	setsize();
}

void TMOGUIImage::fitToWidth(QSize size)
{
	QString s;
	double ratio = 1.0;
	double w = pImage->pSrcPixmap->width();
	double wWorkspace = size.width() - 10;
	ratio = 1.0 - ((w - wWorkspace) / w);
	pImage->Zoom(ratio);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	pImage->SetSize();
	size.setWidth(size.width() - 5);
	size.setHeight(size.height() - 23);
	resize(size);
}

void TMOGUIImage::fitToHeight(QSize size)
{
	QString s;
	double ratio = 1.0;
	double h = pImage->pSrcPixmap->height();
	double hWorkspace = size.height() - 55;
	ratio = 1.0 - ((h - hWorkspace) / h);
	pImage->Zoom(ratio);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	pImage->SetSize();
	size.setWidth(size.width() - 5);
	size.setHeight(size.height() - 23);
	resize(size);
}

void TMOGUIImage::zoomIn()
{
	QString s;
	int currentZoom = pImage->GetRatio();
	currentZoom += ZOOMSTEP;
	if(currentZoom > 1000)
		currentZoom = 1000;
	pImage->Zoom((double)currentZoom / 100);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	pImage->SetSize();
}

void TMOGUIImage::zoomOut()
{
	QString s;
	int currentZoom = pImage->GetRatio();
	currentZoom -= ZOOMSTEP;
	if(currentZoom < 10)
		currentZoom = 10;
	pImage->Zoom((double)currentZoom / 100);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	pImage->SetSize();
}

void TMOGUIImage::setsize()
{
	QSize s = QVBox::size();

	if (bMaximized) 
	{
		pImage->SetSize();
		//PositionImage();
		return;
	}
	pImage->SetSize();
	s.setWidth(pImage->width()+4);
	if (!bTools) 
		s.setHeight(pImage->height()+ADJUST_MIN_EXT);
	else 
	{
		if (s.width() - 4 < MIN_ADJUSTWIDTH) s.setWidth(MIN_ADJUSTWIDTH + 4);
		s.setHeight(pImage->height()+ADJUST_MIN_EXT+pAdjust->height());
	}

    if (s.width() < MIN_WIDTH) s.setWidth(MIN_WIDTH);
	if (s.height() < MIN_HEIGHT) s.setHeight(MIN_HEIGHT);

	resize(s);
}

void TMOGUIImage::resizeEvent ( QResizeEvent * re )
{
	if (re->size() == pParent->size()) bMaximized = true;
	else bMaximized = false;
	QVBox::resizeEvent(re);
}

void TMOGUIImage::showtools()
{
	if (bTools)
	{
		bTools = false;
		pToolsButton->setText("//...\\\\");
		QToolTip::add(pToolsButton, "Open Histogram");
		pAdjust->hide();
	}
	else
	{
		bTools = true;
		pToolsButton->setText("\\\\...//");
		QToolTip::add(pToolsButton, "Close Histogram");
		pAdjust->show();
		QSize sizeAdjust = pAdjust->sizeHint();
		if (sizeAdjust.width() > width()) resize(sizeAdjust.width() + 4, height());	
	}
}

int TMOGUIImage::Swap()
{
	TMOImage *temp = pDst;
	pDst = pSrc;
	pSrc = temp;
	pTransform->Assign(pSrc);
	pOutput->Assign(pSrc);
	pAdjust->SetImage(pSrc);
	pImage->SetImage(pSrc);
	pSrc->WriteLine(L"");
	pSrc->WriteLine(L"Undo last command.");
	return 0;
}

TMOImage* TMOGUIImage::GetImage()
{
	return pSrc;
}

int TMOGUIImage::SetImage(TMOImage *pDest)
{
	if (pDst) delete pDst;
	pDst = pSrc;
	pSrc = pDest;
	pSrc->SetProgress(pTransform->ProgressBar);
	pSrc->SetWriteLine(pTransform->WriteLine);
	pInitProgress = 0;
	pTransform->Assign(pSrc);
	pOutput->Assign(pSrc);
	pAdjust->SetImage(pSrc);
	pImage->SetImage(pSrc);
	pImage->SetSize();
	return 0;
}

const wchar_t* TMOGUIImage::GetString(const QChar* pChars)
{
	for (int i= 0; i < 256; i++)
	{
		buffer[i] = pChars[i];
	}
	return buffer;
}

int TMOGUIImage::SetImageSize(int iWidth, int iHeight)
{
	QString sName;
	TMOImage *pNewImage;
	
	pNewImage = new TMOImage;
	pNewImage->New(iWidth, iHeight);

	double dXRat = (double)pSrc->GetWidth() / iWidth;
	double dYRat = (double)pSrc->GetHeight() / iHeight;
	double *pSrcPixel, *pDstPixel;
	for (int i = 0; i < iHeight; i++)
		for (int j = 0; j < iWidth; j++)
		{
			pSrcPixel = pSrc->GetPixel(dXRat * j, dYRat * i);
			pDstPixel = pNewImage->GetPixel(j, i);
			pDstPixel[0] = pSrcPixel[0];
			pDstPixel[1] = pSrcPixel[1];
			pDstPixel[2] = pSrcPixel[2];		

		}
		
	SetImage(pNewImage);
	return 0;
}

int TMOGUIImage::Extract(TMOGUIImage *pSrcImage, int iComponent)
{
	QString s;
	s = GetName(pSrcImage->name()) + "\n";
	pImage->AddString(s);
	pImage->AddString("Loading ...");
	pSrc->ProgressBar(0, 100);
	adjustSize();
	TMOImage *pLum = 0; 
	if (iComponent == 3) pLum = new TMOImage;
	else pLum = pSrcImage->GetImage();
	try
	{
		if (iComponent == 3) pLum->New(*pSrcImage->GetImage(), TMO_Yxy, true);
		pSrc->New(pLum->GetWidth(), pLum->GetHeight());
	}
	catch ( int e )
	{
		if (e == -19) 
		{
			// Dialog appearing
			pSrc->ProgressBar(0, 0);
			return 1;
		}
	}

	double *pPixel;
	if (pSrc && pLum)
	{
		for (int i = 0; i < pSrc->GetHeight(); i++)
			for (int j = 0; j < pSrc->GetWidth(); j++)
			{
				pPixel = pSrc->GetPixel(j, i);
				if ((iComponent == 0) || (iComponent == 3)) pPixel[0] = pLum->GetPixel(j, i)[iComponent];
				else pPixel[0] = 0.0;
				if ((iComponent == 1) || (iComponent == 3)) pPixel[1] = pLum->GetPixel(j, i)[iComponent];
				else pPixel[1] = 0.0;
				if ((iComponent == 2) || (iComponent == 3)) pPixel[2] = pLum->GetPixel(j, i)[iComponent];
				else pPixel[2] = 0.0;
			}
	}

	s = QString("File created from : ") + pSrcImage->name();
	pSrc->WriteLine(GetString(s.unicode()));
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Computing");
	pAdjust->Create(pSrc, &values);
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Displaying");
	pImage->Create(pSrc,&values,&filters,size,pProgress);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	setsize();
	pToolsButton->show();
	pZoom->show();
	
	pInitProgress = 0;
	show();
	return 0;
}

int TMOGUIImage::MergeComponents(TMOGUIImage* pRed, TMOGUIImage* pGreen, TMOGUIImage* pBlue)
{
	QString s;
	s = name();
	pImage->AddString(s);
	pImage->AddString("Creating ...");
	pSrc->ProgressBar(0, 100);
	adjustSize();
	int iWidth = 0, iHeight = 0;

	iWidth = pRed->GetImage()->GetWidth();
	iHeight = pRed->GetImage()->GetHeight();
	if (pGreen->GetImage()->GetWidth() > iWidth) iWidth = pGreen->GetImage()->GetWidth();
	if (pGreen->GetImage()->GetHeight() > iHeight) iHeight = pGreen->GetImage()->GetHeight();
	if (pBlue->GetImage()->GetWidth() > iWidth) iWidth = pBlue->GetImage()->GetWidth();
	if (pBlue->GetImage()->GetHeight() > iHeight) iHeight = pBlue->GetImage()->GetHeight();

	try
	{
		pSrc->New(iWidth, iHeight);
	}
	catch ( int e )
	{
		if (e == -19) 
		{
			// Dialog appearing
			pSrc->ProgressBar(0, 0);
			return 1;
		}
	}

	double *pPixel;
	TMOImage *pR = pRed->GetImage(), *pG = pGreen->GetImage(), *pB = pBlue->GetImage();
	if (pR && pG && pB)
	{
		for (int i = 0; i < iHeight; i++)
			for (int j = 0; j < iWidth; j++)
			{
				pPixel = pSrc->GetPixel(j, i);
				if ((pR->GetWidth() > j) && (pR->GetHeight() > i)) pPixel[0] = pR->GetPixel(j, i)[0];
				else pPixel[0] = 0;
				if ((pG->GetWidth() > j) && (pG->GetHeight() > i)) pPixel[1] = pG->GetPixel(j, i)[1];
				else pPixel[1] = 0;
				if ((pB->GetWidth() > j) && (pB->GetHeight() > i)) pPixel[2] = pB->GetPixel(j, i)[2];
				else pPixel[2] = 0;
			} 
	}

	s = QString("File created from : \n") + pRed->name() + "\n"  + pGreen->name() + "\n"  + pBlue->name() + "\n";
	pSrc->WriteLine(GetString(s.unicode()));
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Computing");
	pAdjust->Create(pSrc, &values);
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Displaying");
	pImage->Create(pSrc,&values,&filters,size,pProgress);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	setsize();
	pToolsButton->show();
	pZoom->show();
	
	pInitProgress = 0;
	show();
	return 0;
}

int TMOGUIImage::ImageOperation(TMOGUIImage* pRed, TMOGUIImage* pGreen, int iOperation)
{
	QString s;
	s = name();
	pImage->AddString(s);
	pImage->AddString("Creating ...");
	pSrc->ProgressBar(0, 100);
	adjustSize();
	int iWidth = 0, iHeight = 0;

	iWidth = pRed->GetImage()->GetWidth();
	iHeight = pRed->GetImage()->GetHeight();
	if (pGreen->GetImage()->GetWidth() > iWidth) iWidth = pGreen->GetImage()->GetWidth();
	if (pGreen->GetImage()->GetHeight() > iHeight) iHeight = pGreen->GetImage()->GetHeight();

	try
	{
		pSrc->New(iWidth, iHeight);
	}
	catch ( int e )
	{
		if (e == -19) 
		{
			// Dialog appearing
			pSrc->ProgressBar(0, 0);
			return 1;
		}
	}

	double *pPixel;
	double *pRPixel, *pGPixel;
	double pBlack[3] = {0.0, 0.0, 0.0};
	double pWhite[3] = {1.0, 1.0, 1.0};
	TMOImage *pR = pRed->GetImage(), *pG = pGreen->GetImage();
	if (pR && pG)
	{
		for (int i = 0; i < iHeight; i++)
			for (int j = 0; j < iWidth; j++)
			{
				pPixel = pSrc->GetPixel(j, i);
				if ((pR->GetWidth() > j) && (pR->GetHeight() > i)) pRPixel = pR->GetPixel(j, i);
				else pRPixel = 0;
				if ((pG->GetWidth() > j) && (pG->GetHeight() > i)) pGPixel = pG->GetPixel(j, i);
				else pGPixel = 0;

				switch(iOperation)
				{
				case 0:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pGPixel = pBlack;
					pPixel[0] = pRPixel[0] + pGPixel[0];
					pPixel[1] = pRPixel[1] + pGPixel[1];
					pPixel[2] = pRPixel[2] + pGPixel[2];
					break;
				case 1:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pGPixel = pBlack;
					pPixel[0] = pRPixel[0] - pGPixel[0];
					pPixel[1] = pRPixel[1] - pGPixel[1];
					pPixel[2] = pRPixel[2] - pGPixel[2];
					break;
				case 2:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pGPixel = pBlack;
					pPixel[0] = pRPixel[0] * pGPixel[0];
					pPixel[1] = pRPixel[1] * pGPixel[1];
					pPixel[2] = pRPixel[2] * pGPixel[2];
					break;
				case 3:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pGPixel = pWhite;
					pPixel[0] = pRPixel[0] / pGPixel[0];
					pPixel[1] = pRPixel[1] / pGPixel[1];
					pPixel[2] = pRPixel[2] / pGPixel[2];
					break;
				case 4:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pRPixel = pWhite;
					pPixel[0] = pow(pRPixel[0], pGPixel[0]);
					pPixel[1] = pow(pRPixel[1], pGPixel[1]);
					pPixel[2] = pow(pRPixel[2], pGPixel[2]);
					break;
				case 5:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pGPixel = pBlack;
					pPixel[0] = pRPixel[0] > pGPixel[0] ? pRPixel[0] : pGPixel[0];
					pPixel[1] = pRPixel[1] > pGPixel[1] ? pRPixel[1] : pGPixel[1];
					pPixel[2] = pRPixel[2] > pGPixel[2] ? pRPixel[2] : pGPixel[2];
					break;
				case 6:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pGPixel = pBlack;
					pPixel[0] = pRPixel[0] < pGPixel[0] ? pRPixel[0] : pGPixel[0];
					pPixel[1] = pRPixel[1] < pGPixel[1] ? pRPixel[1] : pGPixel[1];
					pPixel[2] = pRPixel[2] < pGPixel[2] ? pRPixel[2] : pGPixel[2];
					break;
				case 7:
					if (!pRPixel) pRPixel = pBlack;
					if (!pGPixel) pGPixel = pBlack;
					pPixel[0] = (pRPixel[0] + pGPixel[0]) / 2;
					pPixel[1] = (pRPixel[1] + pGPixel[1]) / 2;
					pPixel[2] = (pRPixel[2] + pGPixel[2]) / 2;
					break;

				}
			} 
	}

	s = QString("File created from : \n") + pRed->name() + "\n"  + pGreen->name() + "\n";
	pSrc->WriteLine(GetString(s.unicode()));

	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Computing");
	pAdjust->Create(pSrc, &values);
	pSrc->ProgressBar(100, 100);
	pInitProgress->SetLabel("Displaying");
	pImage->Create(pSrc,&values,&filters,size,pProgress);
	pZoom->setText(s.setNum(pImage->GetRatio())+'%');
	setsize();
	pToolsButton->show();
	pZoom->show();
	
	pInitProgress = 0;
	show();
	return 0;
}

TMOGUITransformation* TMOGUIImage::Transform()
{
	if (bTransforming) return 0;
	pProgress->SetProgress(0,100);
	pTransformLabel->show();
	bTransforming = true;
	return pTransform;
}

void TMOGUIImage::customEvent( QCustomEvent * e )
{
	if ( e->type() == QEvent::User )
	{	// finish transform
		pTransformLabel->hide();
		pProgress->SetProgress(0,0);
		bTransforming = false;
		pAdjust->setlinear();
		pAdjust->pHisto->setlinear();
		pAdjust->pToneSlider->setlinear();
		SetImage(pTransform->GetDest());
		pAdjust->setblack("0.0");
		pAdjust->setwhite("1.0");
		pAdjust->pToneSlider->resetsliders();
		pAdjust->pToneSlider->update();
		emit finishTransform();
	}
	else if ( e->type() == QEvent::User + 1 )
	{	// progress
		// int iValue = (int)e->data();
		int iValue = 0; 
		if(e->data() != NULL){
		  iValue = (*((int*)(e->data())));
		}
		if (pInitProgress)
		{
			if (pInitProgress->SetProgress(iValue, 100)) pTransform->Cancel();
		}
		else
		{
			if (pProgress->SetProgress(iValue, 100)) Terminate();
		}
		if (!qApp->hasPendingEvents()) pTransform->refresh->wakeOne();
		if (!iValue) update();
	}
	else if ( e->type() == QEvent::User + 2 )
	{	// writeline
		pOutput->WriteLine(pSrc, (wchar_t*)e->data());
	}
}

int TMOGUIImage::Terminate()
{
	if (iCounter) return 0;

	if (bTransforming)
	{
		QWaitCondition wait;
		pTransform->Cancel();
		pTransformLabel->setText("Closing...");
		pTransformLabel->show();
		while (iCounter < 10) 
		{
			qApp->processEvents();
			if (bTransforming) 
			{
				wait.wait(100);
			}
			else 
			{
				iCounter = 0;
				pTransformLabel->setText("Mapping...");
				return 0;
			}
			iCounter++;
		}
		iCounter = 0;
		pTransform->terminate();
		pTransform->wait();
		pTransformLabel->hide();
		pTransformLabel->setText("Mapping...");
		bTransforming = false;
		pProgress->SetProgress(100,100);
		delete pTransform;
		pTransform = new TMOGUITransformation(this);
	}
	return 0;
}

void TMOGUIImage::canceltransform()
{
	QApplication::postEvent( this, new QCustomEvent((QEvent::Type)(QEvent::User + 1), (void*)0) );
}
bool TMOGUIImage::CanUndo(void)
{
	return pDst;
}
