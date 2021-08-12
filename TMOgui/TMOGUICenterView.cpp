// TMOGUICenterView.cpp: implementation of the TMOGUICenterView class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUICenterView.h"
//Added by qt3to4:
#include <QPaintEvent>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUICenterView::TMOGUICenterView(QWidget * parent, const char * name)
    :QScrollArea(parent)
{
	bZooming = false;
    setObjectName(name);
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
        QScrollArea::paintEvent (pe);
	}
	else 
	{
		int a = 0;
	}
}
