// TMOGUIHisto.cpp: implementation of the TMOGUIHisto class.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"
#include "TMOGUIHisto.h"
#include "TMOGUIProgressBar.h"
#include "TMOGUIAdjustValues.h"
#include <qpainter.h>
#include <qpixmap.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIHisto::TMOGUIHisto(QWidget* parent, const char * name):
	QWidget(parent, name, WRepaintNoErase | WResizeNoErase)
{
	bLog = true;
	iMode = 0;
	iMaxCount = 0;
	setFixedHeight(64);
	pBackBuffer = 0;
	pSrc = 0;
	dScale = 1;
}

TMOGUIHisto::~TMOGUIHisto()
{
	if (pBackBuffer) delete pBackBuffer;
}
void TMOGUIHisto::resizeEvent ( QResizeEvent * re)
{
	if (pBackBuffer) pBackBuffer->resize(re->size());
	QWidget::resizeEvent(re);
}

void TMOGUIHisto::paintEvent ( QPaintEvent * pe)
{
	QPainter p(pBackBuffer);
	QSize s = size();
	int i, l, iBlack, iWhite;

	p.setClipRect(pe->rect());

	switch (iMode)
	{
	case 0:
		iBlack = mapfrom(0.299 * pValues->dRMinimum + 0.587 * pValues->dGMinimum + 0.114 * pValues->dBMinimum) * s.width();
		iWhite = mapfrom(0.299 * pValues->dRMaximum + 0.587 * pValues->dGMaximum + 0.114 * pValues->dBMaximum) * s.width();
		for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
				p.setPen(QColor(220, 220, 220));
			else 
				p.setPen(QColor(255, 255, 255));
			l = 63.0 * pow((double)pLuminance[i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);
			p.moveTo(i,0);
			p.lineTo(i,63-l);
			p.setPen(QColor(128, 128, 128));
			p.lineTo(i,63);
		}
		p.setPen(QColor(255, 0, 0));
		p.moveTo(0,63.0 * (1 - pow((double)pComponents[0][0] / iMaxCount, dScale)));
		for (i = 1; i < HISTOGRAM_WIDTH; i++) 
			p.lineTo(i * s.width() / HISTOGRAM_WIDTH, 63.0 * (1 - pow((double)pComponents[0][i] / iMaxCount, dScale)));
		p.setPen(QColor(0, 255, 0));
		p.moveTo(0,63.0 * (1 - pow((double)pComponents[1][0] / iMaxCount, dScale)));
		for (i = 1; i < HISTOGRAM_WIDTH; i++) 
			p.lineTo(i * s.width() / HISTOGRAM_WIDTH, 63.0 * (1 - pow((double)pComponents[1][i] / iMaxCount, dScale)));
		p.setPen(QColor(0, 0, 255));
		p.moveTo(0,63.0 * (1 - pow((double)pComponents[2][0] / iMaxCount, dScale)));
		for (i = 1; i < HISTOGRAM_WIDTH; i++) 
			p.lineTo(i * s.width() / HISTOGRAM_WIDTH, 63.0 * (1 - pow((double)pComponents[2][i] / iMaxCount, dScale)));
		break;
	case 1:
		iBlack = mapfrom(pValues->dRMinimum) * s.width();
		iWhite = mapfrom(pValues->dRMaximum) * s.width();
		for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
				p.setPen(QColor(220, 220, 220));
			else 
				p.setPen(QColor(255, 255, 255));
			l = 63.0 * pow((double)pComponents[0][i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);
			p.moveTo(i,0);
			p.lineTo(i,63-l);
			p.setPen(QColor(255, 0, 0));
			p.lineTo(i,63);
		}
		break;
	case 2:
		iBlack = mapfrom(pValues->dGMinimum) * s.width();
		iWhite = mapfrom(pValues->dGMaximum) * s.width();
		for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
				p.setPen(QColor(220, 220, 220));
			else 
				p.setPen(QColor(255, 255, 255));
			l = 63.0 * pow((double)pComponents[1][i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);
			p.moveTo(i,0);
			p.lineTo(i,63-l);
			p.setPen(QColor(0, 255, 0));
			p.lineTo(i,63);
		}
		break;
	case 3:
		iBlack = mapfrom(pValues->dBMinimum) * s.width();
		iWhite = mapfrom(pValues->dBMaximum) * s.width();
		for (i = 0; i < s.width(); i++)
		{
			if (i < iBlack || i > iWhite) 
				p.setPen(QColor(220, 220, 220));
			else 
				p.setPen(QColor(255, 255, 255));
			l = 63.0 * pow((double)pComponents[2][i * HISTOGRAM_WIDTH / s.width()] / iMaxCount, dScale);
			p.moveTo(i,0);
			p.lineTo(i,63-l);
			p.setPen(QColor(0, 0, 255));
			p.lineTo(i,63);
		}
		break;
	}
	bitBlt(this, 0, 0, pBackBuffer, 0, 0, s.width(), s. height(), CopyROP);	
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
	if (bLog) pValues->dMaximum = log(pValues->dMaximum);
	//if (bLog) pValues->dMinimum = log(pValues->dMinimum);
	double abs_max=pValues->dMaximum<0 ? pValues->dMinimum=pValues->dMaximum,-pValues->dMaximum: pValues->dMaximum;


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
	pSrc = pImage;
	compute();
	return 0;
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
