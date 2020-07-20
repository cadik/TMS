// TMOGUIHisto.cpp: implementation of the TMOGUIHisto class.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"
#include "TMOGUIHisto.h"
#include "TMOGUIProgressBar.h"
#include "TMOGUIAdjustValues.h"
#include <qpainter.h>
#include <qpixmap.h>
//Added by qt3to4:
#include <QResizeEvent>
#include <QPaintEvent>
#include <QGraphicsView>
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIHisto::TMOGUIHisto(QWidget* parent, const char * name):
    QWidget(parent)
{
    setAttribute(Qt::WA_NoBackground);

    //QGridLayout* layout = new QGridLayout;
    //this->setLayout(layout);
	bLog = true;
	iMode = 0;
	iMaxCount = 0;
	setFixedHeight(64);
    pBackBuffer = nullptr;
	pSrc = 0;
    dScale = 1;
}

TMOGUIHisto::~TMOGUIHisto()
{
	if (pBackBuffer) delete pBackBuffer;
}
void TMOGUIHisto::resizeEvent ( QResizeEvent * re)
{
    if (pBackBuffer != nullptr) *pBackBuffer = pBackBuffer->scaled(re->size()); // TODO check resize
	QWidget::resizeEvent(re);
}

void TMOGUIHisto::paintEvent ( QPaintEvent * pe)
{
    QRect   rcBounds=pe->rect();
    QPainter p;//(pBackBuffer);
    //QPainter out(this);


    QSize s = rcBounds.size();

    QPen lightGray(QColor(128, 128, 128));
    QPen penGray(QColor(220, 220, 220));//, 30, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin); //TODO check
    QPen penWhite(QColor(255, 255, 255));//, 30, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin); //TODO check
    QPen penRed(QColor(255, 0, 0));
    QPen penGreen(QColor(0, 255, 0));
    QPen penBlue(QColor(0, 0, 255));

    QPainterPath pathRed;
    QPainterPath pathGreen;
    QPainterPath pathBlue;

	int i, l, iBlack, iWhite;


    p.begin(pBackBuffer);
    p.setClipRect(rcBounds);
    p.setRenderHint(QPainter::Antialiasing);
    p.setCompositionMode(QPainter::CompositionMode::CompositionMode_SourceOver);

	switch (iMode)
	{
    case 0: //Lum
		iBlack = mapfrom(0.299 * pValues->dRMinimum + 0.587 * pValues->dGMinimum + 0.114 * pValues->dBMinimum) * s.width();
		iWhite = mapfrom(0.299 * pValues->dRMaximum + 0.587 * pValues->dGMaximum + 0.114 * pValues->dBMaximum) * s.width();
        for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
                p.setPen(penGray);
			else 
                p.setPen(penWhite);
            l = 63.0 * pow((double)pLuminance[i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);

            p.drawLine(i,0,i,63-l);
            p.setPen(lightGray);
            p.drawLine(i,63-l,i,63);

		}
        p.setPen(penRed);
        pathRed.moveTo(0,63.0 * (1 - pow((double)pComponents[0][0] / iMaxCount, dScale)));
		for (i = 1; i < HISTOGRAM_WIDTH; i++) 
            pathRed.lineTo(i * s.width() / HISTOGRAM_WIDTH, 63.0 * (1 - pow((double)pComponents[0][i] / iMaxCount, dScale)));
        p.drawPath(pathRed);

        p.setPen(penGreen);
        pathGreen.moveTo(0,63.0 * (1 - pow((double)pComponents[1][0] / iMaxCount, dScale)));
		for (i = 1; i < HISTOGRAM_WIDTH; i++) 
            pathGreen.lineTo(i * s.width() / HISTOGRAM_WIDTH, 63.0 * (1 - pow((double)pComponents[1][i] / iMaxCount, dScale)));
        p.drawPath(pathGreen);

        p.setPen(penBlue);
        pathBlue.moveTo(0,63.0 * (1 - pow((double)pComponents[2][0] / iMaxCount, dScale)));
		for (i = 1; i < HISTOGRAM_WIDTH; i++) 
            pathBlue.lineTo(i * s.width() / HISTOGRAM_WIDTH, 63.0 * (1 - pow((double)pComponents[2][i] / iMaxCount, dScale)));
        p.drawPath(pathBlue);
        break;
    case 1: //Red
		iBlack = mapfrom(pValues->dRMinimum) * s.width();
		iWhite = mapfrom(pValues->dRMaximum) * s.width();
		for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
                p.setPen(penGray);
			else 
                p.setPen(penWhite);
			l = 63.0 * pow((double)pComponents[0][i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);
            /*p.clipPath().moveTo(i,0);
            p.clipPath().lineTo(i,63-l);*/
            p.drawLine(i, 0, i, 63-l);
            p.setPen(penRed);
            //p.clipPath().lineTo(i,63);
            p.drawLine(i, 63-l, i, 63);
		}
		break;
    case 2://Green
		iBlack = mapfrom(pValues->dGMinimum) * s.width();
		iWhite = mapfrom(pValues->dGMaximum) * s.width();
		for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
                p.setPen(penGray);
			else 
                p.setPen(penWhite);
			l = 63.0 * pow((double)pComponents[1][i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);
            /*p.clipPath().moveTo(i,0);
            p.clipPath().lineTo(i,63-l);*/
            p.drawLine(i, 0, i, 63-l);
            p.setPen(penGreen);
            //p.clipPath().lineTo(i,63);
            p.drawLine(i, 63-l, i, 63);
		}
		break;
    case 3://Blue
		iBlack = mapfrom(pValues->dBMinimum) * s.width();
		iWhite = mapfrom(pValues->dBMaximum) * s.width();
		for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
                p.setPen(penGray);
			else 
                p.setPen(penWhite);
			l = 63.0 * pow((double)pComponents[2][i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);
            /*p.clipPath().moveTo(i,0);
            p.clipPath().lineTo(i,63-l);*/
            p.drawLine(i, 0, i, 63-l);
            p.setPen(penBlue);
            //p.clipPath().lineTo(i,63);
            p.drawLine(i, 63-l, i, 63);
		}
		break;
	}
    p.end();

    //bitBlt(this, 0, 0, pBackBuffer, 0, 0, s.width(), s. height(), 1 );
    p.begin(this);
    p.drawPixmap(0, 0, *pBackBuffer, 0, 0, s.width(), s.height());
    p.end();

    //pBackBuffer->grabWidget(this, 0, 0, s.width(), s.height());


}

void TMOGUIHisto::compute()
{
	int i, j, length = pSrc->GetWidth() * pSrc->GetHeight(), index;
	double* offset = pSrc->GetData(), l, p[3];

	pValues->dMinimum = 0;
	if (length) pValues->dMaximum = 0.299 * offset[0]; 
	
	for (i = 0; i < length * 3; i++, offset++) 
	{
		//pValues->dMinimum = pValues->dMinimum > *offset ? *offset : pValues->dMinimum;
		pValues->dMaximum = pValues->dMaximum < *offset ? *offset : pValues->dMaximum;
	}
    if(pValues->dMaximum <= 1 && !bLog) pValues->dMaximum = 1;
	if (bLog) pValues->dMaximum = log(pValues->dMaximum);
	//if (bLog) pValues->dMinimum = log(pValues->dMinimum);
    double abs_max=pValues->dMaximum<0 ? pValues->dMinimum=pValues->dMaximum,pValues->dMaximum: pValues->dMaximum;


	do
	{
		iMaxCount = 0;
		for (i = 0; i < HISTOGRAM_WIDTH; i++) 
		{
			pLuminance[i] = 0;
			for (j = 0; j < COMPONENTS_COUNT; j++) pComponents[j][i] = 0;
		}
		if (bLog) 
		{
		 //pValues->dMinimum -= pValues->dMaximum;
		 pValues->dMinimum -= abs_max;
		}
		offset = pSrc->GetData();
		for (i = 0; i < length; i++)
		{
			if (i%1000 == 0) if (pSrc->ProgressBar(i, length)) 
			{
				pSrc->ProgressBar(0, 0);
				return;
			}
			p[0] = *offset++;
			p[1] = *offset++;
			p[2] = *offset++;
			l = 0.299 * p[0];
			l += 0.587 * p[1];
			l += 0.114 * p[2];
			if (bLog) 
				if (l > 0) l = log(l);
				else 
				{
					l = 0;
				}
			index = HISTOGRAM_WIDTH * (l - pValues->dMinimum) / (pValues->dMaximum - pValues->dMinimum);
			index = index > 0 ? index : 0;
			index = index < HISTOGRAM_WIDTH ? index : HISTOGRAM_WIDTH - 1;
			if ((int)++(pLuminance[index]) > iMaxCount) 
			{
				if (bLog) pValues->dExtreme = exp(l - pValues->dMinimum);
				else pValues->dExtreme = l - pValues->dMinimum;
				iMaxCount = pLuminance[index];
			}
			for (j = 0; j < COMPONENTS_COUNT; j++)
			{
				if (bLog) p[j] = p[j] > 0 ? log(p[j]) : 0;
				index = HISTOGRAM_WIDTH * (p[j] - pValues->dMinimum) / (pValues->dMaximum - pValues->dMinimum);
				index = index > 0 ? index : 0;
				index = index < HISTOGRAM_WIDTH ? index : HISTOGRAM_WIDTH - 1;
				pComponents[j][index]++;
			}
		}

	} while (bLog&&(pLuminance[0] == iMaxCount)&&(pLuminance[0] < .9999 * length));
	pSrc->ProgressBar(0, 0);
}

int TMOGUIHisto::Create(TMOImage *pImage, TMOGUIAdjustValues* pVals)
{
	pValues = pVals;
    pBackBuffer = new QPixmap(size());
    pBackBuffer->fill();
    pSrc = pImage;
    double min, max, avg;
    pSrc->GetMinMaxAvg(&min, &max, &avg);
    if(max <= 1){
        bLog = false;
    }
	compute();
    return bLog;
}

int TMOGUIHisto::SetWindow(TMOImage* pImage)
{
	pSrc = pImage;
	compute();
	update();
	return 0;
}

void TMOGUIHisto::setscale(int iScale)
{
	dScale = exp(iScale * .05);
	update();
}

void TMOGUIHisto::setr()
{
	iMode = 1;
	update();
}

void TMOGUIHisto::setg()
{
	iMode = 2;
	update();
}

void TMOGUIHisto::setb()
{
	iMode = 3;
	update();
}

void TMOGUIHisto::setl()
{
	iMode = 0;
	update();
}

void TMOGUIHisto::setlinear()
{
	bLog = false;
	compute();
	update();
}

void TMOGUIHisto::setlog()
{
	bLog = true;
	compute();
	update();
}

void TMOGUIHisto::reset(int Mode)
{
	switch(Mode)
	{
	case 1:
		pValues->dRMinimum = 0;
		pValues->dRMaximum = 1;
		pValues->dRGamma =1;
		break;
	case 2:
		pValues->dGMinimum = 0;
		pValues->dGMaximum = 1;
		pValues->dGGamma =1;
		break;
	case 3:
		pValues->dBMinimum = 0;
		pValues->dBMaximum = 1;
		pValues->dBGamma =1;
		break;
	default:
		pValues->dRMinimum = 0;
		pValues->dRMaximum = 1;
		pValues->dRGamma =1;
		pValues->dGMinimum = 0;
		pValues->dGMaximum = 1;
		pValues->dGGamma =1;
		pValues->dBMinimum = 0;
		pValues->dBMaximum = 1;
		pValues->dBGamma =1;
		break;
	}
	update();
}

double TMOGUIHisto::mapfrom(double v)
{
	if (bLog) return (log(v) - pValues->dMinimum) / (pValues->dMaximum - pValues->dMinimum);
	else return (v - pValues->dMinimum) / (pValues->dMaximum - pValues->dMinimum);
}

double TMOGUIHisto::mapto(double v)
{
	if (bLog) return exp(pValues->dMinimum + v * (pValues->dMaximum - pValues->dMinimum));
	else return pValues->dMinimum + v * (pValues->dMaximum - pValues->dMinimum);
}
