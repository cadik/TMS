// TMOGUICenterView.cpp: implementation of the TMOGUICenterView class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUICenterView.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUICenterView::TMOGUICenterView(QWidget * parent, const char * name, WFlags f)
	:QScrollView(parent, name, f)
{
	bZooming = false;
}

TMOGUICenterView::~TMOGUICenterView()
{

}

void TMOGUICenterView::Zoom(bool bZoom)
{
	bZooming = bZoom;
}

void TMOGUICenterView::viewportPaintEvent ( QPaintEvent * pe )
{
	if (!bZooming) 
	{
		//QScrollView::paintEvent (pe);
	}
	else 
	{
		int a = 0;
	}
}