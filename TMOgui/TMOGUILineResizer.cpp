
#include "TMOGUILineResizer.h"
#include <qcursor.h>
//Added by qt3to4:
#include <QMouseEvent>
#include <Q3Frame>

TMOGUILineResizer::TMOGUILineResizer( QWidget* parent, const char* name, Qt::WFlags fl )
    : Q3Frame( parent, name, fl )
{
	pParent = parent;
    this->setCursor( QCursor( 6 ) );
    this->setMouseTracking( TRUE );
    this->setAcceptDrops( TRUE );
    this->setFrameShape( Q3Frame::VLine );
	this->setFrameShadow( Q3Frame::Raised );
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