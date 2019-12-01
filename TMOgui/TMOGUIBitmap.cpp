// TMOGUIBitmap.cpp: implementation of the TMOGUIBitmap class.
//
//////////////////////////////////////////////////////////////////////
#include "TMOGUIBitmap.h"
#include "TMOGUIInfoTool.h"
#include "TMOGUIProgressBar.h"
#include "TMOGUIAdjustValues.h"
#include "TMOGUIResource.h"
#include "../tmolib/TMO.h"
#include <qpainter.h>
#include <qpixmap.h>
#include <qimage.h>
#include <qlabel.h>
#include <qcursor.h>
#include <qapplication.h>
//Added by qt3to4:
#include <QContextMenuEvent>
#include <QPaintEvent>
#include <QMouseEvent>
#include <QEvent>
#include <math.h>
#include <float.h>


#ifdef LINUX
#define _isnan isnan
#endif


#define MAXLENGTH 4000.0
#define INITIAL_BITMAP_WIDTH 220
#define INITIAL_BITMAP_HEIGHT 109

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIBitmap::TMOGUIBitmap(QWidget * parent, const char * name):
	QWidget(parent, name, Qt::WNoAutoErase | Qt::WResizeNoErase)
{
	setFixedSize(INITIAL_BITMAP_WIDTH, INITIAL_BITMAP_HEIGHT);
	dRatio = 1;
	iInitial = 0;
	pPixmap = 0;
	pSrcPixmap = 0;
	pSrc = 0;
	bRendering = false;	
	iTool = 0;
	this->setMouseTracking(true);
}

TMOGUIBitmap::~TMOGUIBitmap()
{
	if (pSrcPixmap) delete pSrcPixmap;
	if (pPixmap) delete pPixmap;
}

void TMOGUIBitmap::paintEvent ( QPaintEvent * pe)
{
	QPainter paint(this);
	paint.setClipRect(pe->rect());
	if (pPixmap) paint.drawPixmap(0,0,*pPixmap);
	else 
	{
		paint.drawPixmap(0, 0, *TMOResource::pResource->PixmapLoading->pixmap());
		paint.drawPixmap(19, 72, *TMOResource::pResource->PixmapProgress->pixmap(), 0, 0, 186 * iInitial / 100, 18);
//		paint.fillRect(pe->rect(), QBrush(QColor(255,255,255)));
//		paint.drawText(0,0,INITIAL_BITMAP_WIDTH,INITIAL_BITMAP_HEIGHT,AlignLeft, sInitial);
	}
	if(iTool)
		iTool->DrawTool(paint);
	paint.end();
}

int TMOGUIBitmap::DrawIt(QPainter* pPainter, double *pMargins, bool bZoom)
{
	double dAspect;
	QRect viewport = pPainter->viewport();
	viewport.setLeft(viewport.left() * pMargins[1]  / 100.0);
	viewport.setRight(viewport.right() * (100.0 - 2.0 * pMargins[2])  / 100.0);
	viewport.setBottom(viewport.bottom() * (100.0 - 2.0 *pMargins[3])  / 100.0);
	viewport.setTop(viewport.top() * pMargins[0]  / 100.0);

	dAspect = (double)pPixmap->width() / viewport.right();
	if ((double)viewport.bottom() * dAspect < (double)pPixmap->height()) 
		dAspect = (double)pPixmap->height() / viewport.bottom();

	if (pPixmap) 
	{
		pPainter->setViewport(viewport);
		if (bZoom) pPainter->setWindow(0, 0, dAspect * viewport.right(), dAspect * viewport.bottom());
		pPainter->drawPixmap(0,0,*pPixmap);
	}
	pPainter->end();
	return 0;
}

int TMOGUIBitmap::Render(bool bRepaint)
{
	int i,j, iWidth = pSrc->GetWidth(), iHeight = pSrc->GetHeight();
	double *offset;
	double p[3];
	QRgb* sl;
	QImage pImage(pSrc->GetWidth(), pSrc->GetHeight(),32);

	if (bRendering) 
	{
		bRendering = false;
		return 0;
	}
	else bRendering = true;

	if (bRepaint)
	{
		iAbove = iBelow = 0;
		iTotal = iHeight * iWidth;
		offset = pSrc->GetOffset(0);
		pValues->dRAverage = pValues->dGAverage = pValues->dBAverage = 0;
		for (i = 0; i < iHeight; i++)
		{
			if (i%10 == 0) pSrc->ProgressBar(i, iHeight);
			if (!bRendering) 
			{
				i = 0;
				offset = pSrc->GetOffset(0);
				bRendering = true;
			}
			sl = (QRgb*)pImage.scanLine(i);
			for(j = 0; j < iWidth; j++)
			{
				p[0] = *offset++;                     // Added
				p[1] = *offset++;                     // Added
				p[2] = *offset++;                     // Added

				// Filters
				bool bFiltering = false;
				if (pFilters->dAverage > 0.0)		// Filters enabled
				{
					if (pFilters->dBAverage > 0.0)	// Component filtering
					{
						if (pFilters->dRAverage < 0.0)		// Inside filters
						{
							if ((p[0] >= pFilters->dRMinimum)&&
								(p[1] >= pFilters->dGMinimum)&&
								(p[2] >= pFilters->dBMinimum)&&
								(p[0] < pFilters->dRMaximum)&&
								(p[1] < pFilters->dGMaximum)&&
								(p[2] < pFilters->dBMaximum))
							{
								bFiltering = true;
							}
						}
						else
						{
							if (((p[0] < pFilters->dRMinimum)&&
								 (p[1] < pFilters->dGMinimum)&&
								 (p[2] < pFilters->dBMinimum))||
								((p[0] >= pFilters->dRMaximum)&&
								 (p[1] >= pFilters->dGMaximum)&&
								 (p[2] >= pFilters->dBMaximum)))
							{
								bFiltering = true;
							}
						}
					}
					else
					{
						double dIntensity = 0.299 * p[0] + 0.587 * p[1] + 0.114 * p[2];
						if (pFilters->dRAverage < 0.0)		// Inside filters
						{
							if ((dIntensity >= pFilters->dMinimum)&&
								(dIntensity < pFilters->dMaximum))
								bFiltering = true;
						}
						else
						{
							if ((dIntensity < pFilters->dMinimum)||
								(dIntensity >= pFilters->dMaximum))
								bFiltering = true;
						}
					}

				}
				
				p[0] = (p[0] - pValues->dRMinimum) / (pValues->dRMaximum - pValues->dRMinimum);
				p[1] = (p[1] - pValues->dGMinimum) / (pValues->dGMaximum - pValues->dGMinimum);
				p[2] = (p[2] - pValues->dBMinimum) / (pValues->dBMaximum - pValues->dBMinimum);
			        	
/*				
				p[0] = p[0] < 0.0 ? 0 : pow(p[0], pValues->dRGamma);
				p[1] = p[1] < 0.0 ? 0 : pow(p[1], pValues->dGGamma);
				p[2] = p[2] < 0.0 ? 0 : pow(p[2], pValues->dBGamma);
*/				

				p[0] = pow(p[0], pValues->dRGamma);
				p[1] = pow(p[1], pValues->dGGamma);
				p[2] = pow(p[2], pValues->dBGamma);
			        
				if(_isnan(p[0]))p[0]=0;
				if(_isnan(p[1]))p[1]=0;
				if(_isnan(p[2]))p[2]=0;

				if (p[0] + p[1] + p[2] <= 0.)
				{
					iBelow++;
				}
			else	
				
				if (p[0] + p[1] + p[2] > 3.)
				{
					iAbove++;
				}


				p[0] = p[0] > 0.0 ? p[0] : 0.0;
				p[0] = p[0] < 1.0 ? p[0] : 1.0;

				
				p[1] = p[1] > 0.0 ? p[1] : 0.0;
				p[1] = p[1] < 1.0 ? p[1] : 1.0;
				
				
				p[2] = p[2] > 0.0 ? p[2] : 0.0;
				p[2] = p[2] < 1.0 ? p[2] : 1.0;

				pValues->dRAverage += p[0];
				pValues->dGAverage += p[1];
				pValues->dBAverage += p[2];

				if (bFiltering)
				{
					p[0] = pFilters->dExtreme * pFilters->dRGamma + 
						(1.0 - pFilters->dExtreme) * p[0];
					p[1] = pFilters->dExtreme * pFilters->dGGamma + 
						(1.0 - pFilters->dExtreme) * p[1];
					p[2] = pFilters->dExtreme * pFilters->dBGamma + 
						(1.0 - pFilters->dExtreme) * p[2];
				}

				p[0] *= 255.99999;
				p[1] *= 255.99999;
				p[2] *= 255.99999;

				sl[j] = qRgb((int)p[0], (int)p[1], (int)p[2]);
			}
		}
	//	std::cout << "iBelow " << iBelow << std::endl;
		pValues->dRAverage /= iTotal;
		pValues->dGAverage /= iTotal;
		pValues->dBAverage /= iTotal;
		pSrc->ProgressBar(0, 0);
		pSrcPixmap->convertFromImage(pImage, Qt::AutoColor|Qt::OrderedDither|Qt::OrderedAlphaDither);
	}
	if ((pSrc->GetWidth()!=s.width())||(pSrc->GetHeight()!=s.height())) 
	{
		QPixmap tempPixmap(*pSrcPixmap);
        QMatrix m;
        m.scale((double)s.width() / pSrcPixmap->width(),
                (double)s.height() / pSrcPixmap->height());
        *pPixmap = tempPixmap.xForm(m);
    }
	else *pPixmap = *pSrcPixmap;
	bRendering = false;
	emit rendered();
	return 0;
}


int TMOGUIBitmap::Create(TMOImage *src, TMOGUIAdjustValues* pVals, TMOGUIAdjustValues* pFilt, QSize window, TMOGUIProgressBar* pProgressBar)
{
	int aspect;	
	pSrc = src;
	pValues = pVals;
	pFilters = pFilt;
	pProgress = pProgressBar;
	if ((window.width()*.9 < src->GetWidth())||(window.height()*.9 < src->GetHeight()))
	{
		if ((window.width() / src->GetWidth())<(window.height() / src->GetHeight())) 
			aspect = 10 * window.width()*.9 / src->GetWidth();
		else
			aspect = 10 * window.height()*.9 / src->GetHeight();
		dRatio = (double)aspect / 10;
	}
	s.setWidth(pSrc->GetWidth()*dRatio);
	s.setHeight(pSrc->GetHeight()*dRatio);
	pPixmap = new QPixmap(s);
	pSrcPixmap = new QPixmap(s);
	Render();
	update();
	setFixedSize(s);
	return 0;
}

int TMOGUIBitmap::Zoom(double ratio)
{
	int length;
	int retval = 0;

	if (s.width() > s.height()) length = pSrc->GetWidth();
	else length = pSrc->GetHeight();
	if ((ratio > .1) && (length * ratio < MAXLENGTH)) dRatio = ratio;
	else 
	{
		if (ratio <= .1) dRatio = .1;
		else dRatio = MAXLENGTH / length;
		retval = 1;
	}
	s.setWidth(pSrc->GetWidth() * dRatio);
	s.setHeight(pSrc->GetHeight() * dRatio);
	Render(false);
	return retval;
}

int TMOGUIBitmap::GetRatio()
{
	return 100 * dRatio + .5;
}

int TMOGUIBitmap::AddString(const QString &s)
{
	sInitial.append(s+"\n");
	update();
	qApp->processEvents();
	return 0;
}

void TMOGUIBitmap::valueschanged()
{
	Render(true);
	update();
}

int TMOGUIBitmap::SetSize(const QSize* size)
{
	if (size) s = *size;
	setFixedSize(s);
	return 0;
}

QSize* TMOGUIBitmap::GetSize()
{
	return &s;
}

int TMOGUIBitmap::SetImage(TMOImage* pImage)
{
	pSrc = pImage;
	s.setWidth(pSrc->GetWidth() * dRatio);
	s.setHeight(pSrc->GetHeight() * dRatio);
	Render(true);
	update();
	return 0;
}

TMOImage* TMOGUIBitmap::GetImage(TMOImage* pDst)
{
	int i,j, iWidth = pSrc->GetWidth(), iHeight = pSrc->GetHeight();
	double *offset;
	double p[3], *dest;
	QRgb* sl;
	QImage pImage(pSrc->GetWidth(), pSrc->GetHeight(),32);

	offset = pSrc->GetOffset(0);
	dest = pDst->GetOffset(0);
	for (i = 0; i < iHeight; i++)
	{
		if (i%10 == 0) pSrc->ProgressBar(i, iHeight);
		sl = (QRgb*)pImage.scanLine(i);
		for(j = 0; j < iWidth; j++)
		{
			p[0] = *offset++;                     // Added
			p[1] = *offset++;                     // Added
			p[2] = *offset++;                     // Added

			p[0] = pow((p[0] - pValues->dRMinimum) / (pValues->dRMaximum - pValues->dRMinimum), pValues->dRGamma);
			p[1] = pow((p[1] - pValues->dGMinimum) / (pValues->dGMaximum - pValues->dGMinimum), pValues->dGGamma);
			p[2] = pow((p[2] - pValues->dBMinimum) / (pValues->dBMaximum - pValues->dBMinimum), pValues->dBGamma);
				
			*dest++ = p[0];
			*dest++ = p[1];
			*dest++ = p[2];
		}
	}
	pSrc->ProgressBar(0,0);
	return pDst;
}

int TMOGUIBitmap::SetInitial(int i)
{
	int oldval = i;
	iInitial = i;
	update();
	qApp->processEvents();
	return oldval;
}

void TMOGUIBitmap::mouseMoveEvent( QMouseEvent * e )
{
	if(iTool)
		iTool->MouseAction(e);
}

void TMOGUIBitmap::contextMenuEvent( QContextMenuEvent * )
{   
	if(iTool)
		iTool->CreateContextMenu();
}

void TMOGUIBitmap::mousePressEvent ( QMouseEvent * e )
{	
	if(iTool)
		iTool->MouseAction(e);
}

void TMOGUIBitmap::leaveEvent ( QEvent * )
{
	if(iTool)
		iTool->LeaveImage();
}

void TMOGUIBitmap::ActivateTool ( TMOGUIInfoTool * ptr )
{
	iTool = ptr;
	this->setCursor(QCursor(Qt::crossCursor));
	iTool->setToolPtr(this);	
}

void TMOGUIBitmap::DeactivateTool ()
{
	iTool = 0;
	this->setCursor(QCursor(Qt::arrowCursor));
}
