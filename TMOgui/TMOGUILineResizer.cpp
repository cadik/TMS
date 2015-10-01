
#include "TMOGUILineResizer.h"
#include <qcursor.h>

TMOGUILineResizer::TMOGUILineResizer( QWidget* parent, const char* name, WFlags fl )
    : QFrame( parent, name, fl )
{
	pParent = parent;
    this->setCursor( QCursor( 6 ) );
    this->setMouseTracking( TRUE );
    this->setAcceptDrops( TRUE );
    this->setFrameShape( QFrame::VLine );
	this->setFrameShadow( QFrame::Raised );
	this->setLineWidth(2);
	pressed = false;
	minWidth = 0;
}

TMOGUILineResizer::~TMOGUILineResizer()
{
}

void TMOGUILineResizer::SetMinWidth(int min)
{
	minWidth = min;
}

void TMOGUILineResizer::mouseMoveEvent( QMouseEvent * e )
{
	if(pressed && ((this->x() + e->x()) >= minWidth))
	{
		this->move(e->x() + this->x(), this->y());
		emit resizeInfoElem(e->x());
	}
}

void TMOGUILineResizer::mousePressEvent ( QMouseEvent * e )
{
	pressed = true;
}

void TMOGUILineResizer::mouseReleaseEvent ( QMouseEvent * e )
{
	pressed = false;
}