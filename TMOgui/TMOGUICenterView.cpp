// TMOGUICenterView.cpp: implementation of the TMOGUICenterView class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUICenterView.h"
//Added by qt3to4:
#include <QPaintEvent>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUICenterView::TMOGUICenterView(QWidget * parent, const char * name, Qt::WFlags f)
	:Q3ScrollView(parent, name, f)
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