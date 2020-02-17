// TMOGUIInfoTool.cpp: implementation of the TMOGUIInfoTool class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIBitmap.h"
#include "ui_resources6.h"
#include "TMOGUIInfoTool.h"
#include <qpainter.h>
#include <qcursor.h>
#include <qlineedit.h>
#include <qscrollbar.h>
#include <qradiobutton.h>
//Added by qt3to4:
#include <QMouseEvent>
#include <QWidget>


TMOGUIInfoTool::TMOGUIInfoTool(QObject * parent, const char * name) : QObject(parent)
{
    toolDialog = new QDialog((QWidget*) parent);
	drawToolIcon = false;
	enableTool = false;
	Reset();
	contextDialogEnable = false;	
    toolContext = new Ui::TMOGUITool();
    toolContext->setupUi(toolDialog); //TODO Qdialog?
	connect( toolContext->scrollBar, SIGNAL( valueChanged(int) ), this, SLOT( changeSize(int) ) );
	connect( toolContext->editSize, SIGNAL( textChanged(const QString &) ), this, SLOT( changeTextSize(const QString &) ) );
	connect( toolContext->radioCircle, SIGNAL( toggled(bool) ), this, SLOT( changeCircle(bool) ) );
	connect( toolContext->radioSquare, SIGNAL( toggled(bool) ), this, SLOT( changeSquare(bool) ) );
    toolDialog->hide();
}

TMOGUIInfoTool::~TMOGUIInfoTool()
{
}

void TMOGUIInfoTool::DrawTool(QPainter & paint)
{
	if(drawToolIcon && enableTool)
	{
		int halfSize = toolSize * parentBitmap->GetRatioD() / 2;
		paint.setPen( Qt::white );
		if(circleShape)
		{
			paint.drawEllipse(x - halfSize, y - halfSize, halfSize + halfSize, halfSize + halfSize);
			paint.setPen( Qt::black );
			paint.drawEllipse(x - halfSize - 1, y - halfSize - 1, halfSize + halfSize + 2, halfSize + halfSize + 2);
		}
		else
		{
			paint.drawRect(x - halfSize, y - halfSize, halfSize + halfSize, halfSize + halfSize);
			paint.setPen( Qt::black );
			paint.drawRect(x - halfSize - 1, y - halfSize - 1, halfSize + halfSize + 2, halfSize + halfSize + 2);
		}
	}
}

void TMOGUIInfoTool::Reset()
{
	toolSize = 30;
	x = y = 0;
	circleShape = true;
}

void TMOGUIInfoTool::MouseAction(QMouseEvent * e)
{
	if(enableTool)
	{
		drawToolIcon = true;
		x = e->x();
		y = e->y();
		emit toolApllied(x, y, toolSize / 2, circleShape);
		parentBitmap->update();
	}
}

void TMOGUIInfoTool::LeaveImage()
{
	if(!contextDialogEnable && enableTool)
	{
		drawToolIcon = false;
		emit toolCancelled();
		parentBitmap->update();
	}
}

void TMOGUIInfoTool::SetEnabled(bool on) 
{ 
	enableTool = on;
	if(on)
        parentBitmap->setCursor(QCursor(Qt::CursorShape::CrossCursor));
	else
        parentBitmap->setCursor(QCursor(Qt::CursorShape::ArrowCursor));
}

void TMOGUIInfoTool::CreateContextMenu()
{
	if(enableTool)
	{
		QString s;
		contextDialogEnable = true;
		toolContext->scrollBar->setValue(toolSize);
		toolContext->editSize->setText(s.setNum(toolSize));
		if(circleShape)
			toolContext->radioCircle->setChecked(true);
		else
			toolContext->radioSquare->setChecked(true);
        // TODO toolContext->move(QCursor::pos());
        // TODO if (toolContext->exec() == QDialog::Rejected) contextDialogEnable = false;
	}
}

void TMOGUIInfoTool::changeTextSize(const QString & s)
{
	bool ok;
	int size = s.toInt(&ok);
	if(!ok || 1 > size || size > 150) return;
	toolSize = size;
	if(parentBitmap)parentBitmap->update();
}

void TMOGUIInfoTool::changeSize(int value)
{
	toolSize = value;
	if(parentBitmap)parentBitmap->update();
}

void TMOGUIInfoTool::changeCircle(bool on)
{
	circleShape = on;
	if(parentBitmap)parentBitmap->update();
}

void TMOGUIInfoTool::changeSquare(bool on)
{
	circleShape = !on;
	parentBitmap->update();
}
