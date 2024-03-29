
#include "TMOGUILineResizer.h"
#include <qcursor.h>
//Added by qt3to4:
#include <QMouseEvent>
#include <QFrame>

TMOGUILineResizer::TMOGUILineResizer(QWidget *parent, const char *name /*, Qt::WindowFlags fl */)
	: QFrame(parent)
{
	pParent = parent;
	// this->setCursor( QCursor( 6 ) );
	this->setMouseTracking(true);
	this->setAcceptDrops(true);
	this->setFrameShape(QFrame::VLine);
	this->setFrameShadow(QFrame::Raised);
	this->setLineWidth(2);
	this->setObjectName(name);
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

void TMOGUILineResizer::mouseMoveEvent(QMouseEvent *e)
{
	if (pressed && ((this->x() + e->x()) >= minWidth))
	{
		this->move(e->x() + this->x(), this->y());
		emit resizeInfoElem(e->x());
	}
}

void TMOGUILineResizer::mousePressEvent(QMouseEvent *e)
{
	pressed = true;
}

void TMOGUILineResizer::mouseReleaseEvent(QMouseEvent *e)
{
	pressed = false;
}
