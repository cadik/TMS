// TMOGUIToneSlider.cpp: implementation of the TMOGUIToneSlider class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIToneSlider.h"
#include "TMOGUIAdjustValues.h"
#include <qpainter.h>
#include <qlayout.h>
#include <qpixmap.h>
#include <math.h>
#include "gamma.h"

//#include <iostream>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIToneSlider::TMOGUIToneSlider(QWidget* parent, const char * name):
	QWidget(parent, name, WRepaintNoErase | WResizeNoErase)
{
	setFixedHeight(20);
	iBlack = iWhite = iGamma = 0;
	iMode = 0;
	pValues = 0;
	iCaptured = 0;
	pBackBuffer = 0;
	bLog = true;
}

TMOGUIToneSlider::~TMOGUIToneSlider()
{
	if (pBackBuffer) delete pBackBuffer;
}

int TMOGUIToneSlider::Create(TMOGUIAdjustValues* pVals)
{
	pValues = pVals;
	pBackBuffer = new QPixmap(1,20);
	setl();
//	std::cout << "TMOGUIToneSlider::Create " << iWhite << std::endl;
	return 0;
}

void TMOGUIToneSlider::paintEvent ( QPaintEvent * pe)
{
	QPainter p(pBackBuffer);
	int col;
//	double val;

	p.setClipRect(pe->rect());

	for (int i = 0; i < s.width(); i++)
	{
		/*val = 256.0 * exp(5 * ((double)i / s.width() - 1));
		if (bLog) col = (int)val;
		else col = i * 256 / s.width();
		*/
		col = 256 * mapto((double)i / (s.width() - 1));
		if (col > 255) col = 255;
		switch(iMode)
		{
		case 1:
			p.setPen(QColor(col, 0, 0));
			break;
		case 2:
			p.setPen(QColor(0, col, 0));
			break;
		case 3:
			p.setPen(QColor(0, 0, col));
			break;
		default:
			p.setPen(QColor(col, col, col));
			break;
		}
		p.drawLine(i,0,i,9);
	}
	p.fillRect(0,10,s.width(),10,QBrush(QColor(255,255,255)));
	DrawMarker(&p, iBlack, QColor(0,0,0));
	DrawMarker(&p, iGamma, QColor(128,128,128));
	DrawMarker(&p, iWhite, QColor(255,255,255));
	p.drawLine(iBlack, 9, iBlack, 7);
	p.drawLine(iWhite, 9, iWhite, 7);
	p.fillRect(iBlack + 1, 3, iWhite - iBlack, 4, QBrush(QColor(255,255,255)));
	p.drawRect(iBlack, 2, iWhite - iBlack + 1, 5);
	bitBlt(this, 0, 0, pBackBuffer, 0, 0, s.width(), s. height(), CopyROP);	
}

void TMOGUIToneSlider::resizeEvent ( QResizeEvent * re )
{
	s = re->size();
	if (pBackBuffer) pBackBuffer->resize(s);

	resetsliders();
	QWidget::resizeEvent(re);
}

void TMOGUIToneSlider::resetsliders()
{
	if (pValues)
	{
		switch (iMode)
		{
		case 0:
			iBlack = mapfrom(0.299 * pValues->dRMinimum + 0.587 * pValues->dGMinimum + 0.114 * pValues->dBMinimum) * s.width();
			iWhite = mapfrom(0.299 * pValues->dRMaximum + 0.587 * pValues->dGMaximum + 0.114 * pValues->dBMaximum) * s.width();
			if(iWhite<0)iWhite=0;
			if(iBlack<0)iBlack=0; //pAn - aby cernej neutek doprava
			iGamma = iBlack + (iWhite - iBlack)*(pValues->dRGamma + pValues->dGGamma + pValues->dBGamma)/6;
			break;
		case 1:
			iBlack = mapfrom(pValues->dRMinimum) * s.width();
			iWhite = mapfrom(pValues->dRMaximum) * s.width();
			if(iWhite<0)iWhite=0;
			if(iBlack<0)iBlack=0;
			iGamma = iBlack + (iWhite - iBlack)*pValues->dRGamma/2;
			break;
		case 2:
			iBlack = mapfrom(pValues->dGMinimum) * s.width();
			iWhite = mapfrom(pValues->dGMaximum) * s.width();
			if(iWhite<0)iWhite=0;
			if(iBlack<0)iBlack=0;
			iGamma = iBlack + (iWhite - iBlack)*pValues->dGGamma/2;
			break;
		case 3:
			iBlack = mapfrom(pValues->dBMinimum) * s.width();
			iWhite = mapfrom(pValues->dBMaximum) * s.width();
			if(iWhite<0)iWhite=0;
			if(iBlack<0)iBlack=0;
			iGamma = iBlack + (iWhite - iBlack)*pValues->dBGamma/2;
			break;
		}
	}
}
void TMOGUIToneSlider::setr()
{
	iBlack = mapfrom(pValues->dRMinimum) * s.width();
	iWhite = mapfrom(pValues->dRMaximum) * s.width();
	if(iWhite<0)iWhite=0;
	if(iBlack<0)iBlack=0;
	iMode = 1;
	correct_iGamma();
	update();
}

void TMOGUIToneSlider::setg()
{
	iBlack = mapfrom(pValues->dGMinimum) * s.width();
	iWhite = mapfrom(pValues->dGMaximum) * s.width();
	if(iWhite<0)iWhite=0;
	if(iBlack<0)iBlack=0;
	iMode = 2;
	correct_iGamma();
	update();
}

void TMOGUIToneSlider::setb()
{
	iBlack = mapfrom(pValues->dBMinimum) * s.width();
	iWhite = mapfrom(pValues->dBMaximum) * s.width();
	if(iWhite<0)iWhite=0;
	if(iBlack<0)iBlack=0;
	iMode = 3;
	correct_iGamma();
	update();
}

void TMOGUIToneSlider::setl()
{
	iBlack = mapfrom(0.299 * pValues->dRMinimum + 0.587 * pValues->dGMinimum + 0.114 * pValues->dBMinimum) * s.width();
	iWhite = mapfrom(0.299 * pValues->dRMaximum + 0.587 * pValues->dGMaximum + 0.114 * pValues->dBMaximum) * s.width();
	if(iWhite<0)iWhite=0;
	if(iBlack<0)iBlack=0;
	iMode = 0;
	correct_iGamma();
	update();
//	std::cout << "TMOGUIToneSlider::setl " << iWhite << " " << iBlack << std::endl;
}

void TMOGUIToneSlider::setlinear()
{
	bLog = false;
	switch (iMode)
	{
	case 0:
		setl();
		break;
	case 1:
		setr();
		break;
	case 2:
		setg();
		break;
	case 3:
		setb();
		break;
	}
	update();
}

void TMOGUIToneSlider::setlog()
{
	bLog = true;
	switch (iMode)
	{
	case 0:
		setl();
		break;
	case 1:
		setr();
		break;
	case 2:
		setg();
		break;
	case 3:
		setb();
		break;
	}
	update();
}

int TMOGUIToneSlider::DrawMarker(QPainter* p, int x, QColor col)
{
	QPointArray pa(3);

	pa.setPoint(0, x-4, 19);
	pa.setPoint(1, x, 10);
	pa.setPoint(2, x+4, 19);
	p->setBrush(QBrush(col));
	p->setPen(QColor(128,128,128));
	p->drawPolygon(pa);
	return 0;
}

void TMOGUIToneSlider::mousePressEvent ( QMouseEvent * e ) 
{
	if (e->y() > 9)
	{
		if ((e->x() > iWhite - 5)&&(e->x() < iWhite + 5)) 
		{
			setMouseTracking(true);
			iCaptured = 2;
			return;
		}
		if ((e->x() > iBlack - 5)&&(e->x() < iBlack + 5)) 
		{
			setMouseTracking(true);
			iCaptured = 1;
			return;
		}
		if ((e->x() > iGamma - 5)&&(e->x() < iGamma + 5)) 
		{
			setMouseTracking(true);
			iCaptured = 3;
			return;
		}
	}
	else
	{
		if ((e->x() > iBlack - 1) && (e->x() < iWhite + 1)) 
		{
			iGrab = e->x();	
			iCaptured = 4;
			setMouseTracking(true);
		}
	}
}

void TMOGUIToneSlider::mouseReleaseEvent ( QMouseEvent * e ) 
{
	if (iCaptured)
	{
		iCaptured = 0;
		setMouseTracking(false);
		pValues->RenderValues();
	}
	return;
}
void TMOGUIToneSlider::mouseDoubleClickEvent ( QMouseEvent * e ) 
{
	if (e->y() > 9)
	{
		if ((e->x() > iWhite - 5)&&(e->x() < iWhite + 5)) 
		{
			iWhite = s.width() - 1;
		}
		else if ((e->x() > iBlack - 5)&&(e->x() < iBlack + 5)) 
		{
			iBlack = 0;
		}
		else if ((e->x() > iGamma - 5)&&(e->x() < iGamma + 5)) 
		{
			iGamma = (iWhite + iBlack) / 2;
			UpdateValues();
			pValues->RenderValues();
			return;
		}
	}
	else if ((e->x() > iBlack - 1) && (e->x() < iWhite + 1)) 
	{
		iWhite = s.width() - 1;
		iBlack = 0;
		iGamma = (iWhite + iBlack) / 2;
		UpdateValues();
		pValues->RenderValues();
		return;
	}
	UpdateValues(false);
	correct_iGamma();
	pValues->RenderValues();
}

void TMOGUIToneSlider::setToWidth()
{
		iWhite = s.width() - 1;
		iBlack = 0;
		iGamma = (iWhite + iBlack) / 2;
		UpdateValues();
		pValues->RenderValues();

}


void TMOGUIToneSlider::mouseMoveEvent ( QMouseEvent * e ) 
{
	switch (iCaptured)
	{
	case 1:
		iBlack = e->x();
		iBlack = iBlack < 0 ? 0 : iBlack;
		iBlack = iBlack <= iWhite ? iBlack : iWhite - 1;
		break;
	case 2:
		iWhite = e->x();
		iWhite = iWhite <= iBlack ? iBlack + 1 : iWhite;
		iWhite = iWhite < s.width() ? iWhite : s.width() - 1;
		break;
	case 3:
		iGamma = e->x();
		iGamma = iGamma <= iBlack ? iBlack + 1 : iGamma;
		iGamma = iGamma < iWhite ? iGamma : iWhite - 1;
		switch(iMode)
		{
		case 0:
		        
			pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
			if ((pValues->dRGamma > .49) && (pValues->dRGamma < .51)) pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = .5;
			pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = log(pValues->dRGamma) / log(0.5);

			if (pValues->dRGamma > INV_GAMMA_MIN) 
			{
				pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = INV_GAMMA_MIN;
				correct_iGamma();
			}
			else if (pValues->dRGamma < INV_GAMMA_MAX) 
			{
				pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = INV_GAMMA_MAX;
				correct_iGamma();
			}
			update();
			pValues->UpdateValues();
			break;
		case 1:
			pValues->dRGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
			if ((pValues->dRGamma > .49) && (pValues->dRGamma < .51)) pValues->dRGamma = .5;
			pValues->dRGamma = log(pValues->dRGamma) / log(0.5);
			if (pValues->dRGamma > INV_GAMMA_MIN) 
			{
				pValues->dRGamma = INV_GAMMA_MIN;
				correct_iGamma();
			}
			else if (pValues->dRGamma < INV_GAMMA_MAX) 
			{
				pValues->dRGamma = INV_GAMMA_MAX;
				correct_iGamma();
			}
			update();
			pValues->UpdateValues();
			break;
		case 2:
			pValues->dGGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
			if ((pValues->dGGamma > .49) && (pValues->dGGamma < .51)) pValues->dGGamma = .5;
			pValues->dGGamma = log(pValues->dGGamma) / log(0.5);
			if (pValues->dGGamma > INV_GAMMA_MIN) 
			{
				pValues->dGGamma = INV_GAMMA_MIN;
				correct_iGamma();
			}
			else if (pValues->dGGamma < INV_GAMMA_MAX) 
			{
				pValues->dGGamma = INV_GAMMA_MAX;
				correct_iGamma();
			}
			update();
			pValues->UpdateValues();
			break;
		case 3:
			pValues->dBGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
			if ((pValues->dBGamma > .49) && (pValues->dBGamma < .51)) pValues->dBGamma = .5;
			pValues->dBGamma = log(pValues->dBGamma) / log(0.5);
			if (pValues->dBGamma > INV_GAMMA_MIN) 
			{
				pValues->dBGamma = INV_GAMMA_MIN;
				correct_iGamma();
			}
			else if (pValues->dBGamma < INV_GAMMA_MAX) 
			{
				pValues->dBGamma = INV_GAMMA_MAX;
				correct_iGamma();
			}
			update();
			pValues->UpdateValues();
			break;
		}
		return;
	case 4:
		iWhite = iWhite - iGrab + e->x();
		iBlack = iBlack - iGrab + e->x();
		iGrab = e->x();
		if (iBlack < 0)
		{
			iWhite = iWhite - iBlack;
			iBlack = 0;
		}
		if (iWhite >= s.width())
		{
			iBlack += s.width() - iWhite - 1;
			iWhite = s.width() - 1;
		}
		break;
	}
	UpdateValues(false);
        correct_iGamma();	
}

void TMOGUIToneSlider::correct_iGamma()//posouva ukazatel gammy po zmene B/W
{
	double temp;
	switch(iMode)
	{
	 	case 1:
		 temp = exp(pValues->dRGamma * log(0.5));
		break;
	 	
		case 2:
		 temp = exp(pValues->dGGamma * log(0.5));
		break;
	 	
		case 3:
		 temp = exp(pValues->dBGamma * log(0.5));
		break;
		
		case 0:
		default:
		 temp =  (exp(log(.5) / (1/pValues->dRGamma + 1/pValues->dGGamma + 1/pValues->dBGamma) * 3));
	}
	iGamma = iBlack + (iWhite - iBlack) * temp; 
}


bool TMOGUIToneSlider::InRect(const QPoint& p, int x1, int x2, int y1, int y2)
{
	return (p.x() >= x1)&&(p.x() <= x2)&&(p.y() >= y1)&&(p.y() <= y2);
}

int TMOGUIToneSlider::UpdateValues(bool bGamma)
{
	if (pValues)
	{
		switch (iMode)
		{
		case 0:
			pValues->dRMinimum = pValues->dGMinimum = pValues->dBMinimum = mapto((double)iBlack / (s.width() - 1));
			pValues->dRMaximum = pValues->dGMaximum = pValues->dBMaximum = mapto((double)iWhite / (s.width() - 1));
			if (bGamma)
			{
				pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
				if ((pValues->dRGamma > .49) && (pValues->dRGamma < .51)) pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = .5;
				pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = log(0.5) / log(pValues->dRGamma);
			}
			break;
		case 1:
			pValues->dRMinimum = mapto((double)iBlack / (s.width() - 1));
			pValues->dRMaximum = mapto((double)iWhite / (s.width() - 1));
			if (bGamma)
			{
				pValues->dRGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
				if ((pValues->dRGamma > .49) && (pValues->dRGamma < .51)) pValues->dRGamma = .5;
				pValues->dRGamma = log(0.5) / log(pValues->dRGamma);
			}
			break;
		case 2:
			pValues->dGMinimum = mapto((double)iBlack / (s.width() - 1));
			pValues->dGMaximum = mapto((double)iWhite / (s.width() - 1));
			if (bGamma)
			{
				pValues->dGGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
				if ((pValues->dGGamma > .49) && (pValues->dGGamma < .51)) pValues->dGGamma = .5;
				pValues->dGGamma = log(0.5) / log(pValues->dGGamma);
			}
			break;
		case 3:
			pValues->dBMinimum = mapto((double)iBlack / (s.width() - 1));
			pValues->dBMaximum = mapto((double)iWhite / (s.width() - 1));
			if (bGamma)
			{
				pValues->dBGamma = (double)(iGamma - iBlack) / (iWhite - iBlack);
				if ((pValues->dBGamma > .49) && (pValues->dBGamma < .51)) pValues->dBGamma = .5;
				pValues->dBGamma = log(0.5) / log(pValues->dBGamma);
			}
			break;
		}
	}
	update();
	pValues->UpdateValues();
	return 0;
}

double TMOGUIToneSlider::mapfrom(double v)
{
	if (bLog) return (log(v) - pValues->dMinimum) / (pValues->dMaximum - pValues->dMinimum);
	else return (v - pValues->dMinimum) / (pValues->dMaximum - pValues->dMinimum);
}

double TMOGUIToneSlider::mapto(double v)
{
	if (bLog) return exp(pValues->dMinimum + v * (pValues->dMaximum - pValues->dMinimum));
	else return pValues->dMinimum + v * (pValues->dMaximum - pValues->dMinimum);
}
