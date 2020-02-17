// TMOGUIZoomTool.cpp: implementation of the TMOGUIZoomTool class.
//
//////////////////////////////////////////////////////////////////////
	
#include "TMOGUIZoomTool.h"
#include "TMOGUIImage.h"
#include "TMOGUIResource.h"
#include <qtoolbutton.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qcombobox.h>
#include <QMenu>
#include <qmenubar.h>
#include <QMdiArea>
#include <qpushbutton.h>




//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIZoomTool::TMOGUIZoomTool(QWidget * parent, const char * name):QToolBar(parent)
{
	pParent = parent;
    setObjectName(name);
	pImage = 0;
	Create();
}

TMOGUIZoomTool::~TMOGUIZoomTool()
{
}

int TMOGUIZoomTool::Create()
{
    this->setWindowTitle( "Zoom manipulations" );
    zoomOutBtn = new QToolButton(this);
    zoomOutBtn->setIcon(QIcon(":/icons/IconMinus.png"));
    zoomOutBtn->setText("Zoom Out");
    zoomOutBtn->setParent( pParent);
    connect(zoomOutBtn, SIGNAL(triggered(QAction *)), this, SLOT( zoomOut() )); // TODO check
	zoomOutBtn->setFixedSize(25, 30);

	zoomEdit = new QLineEdit("100%", this);
	zoomEdit->setFixedWidth(50);
	zoomEdit->setAlignment(Qt::AlignCenter);
	connect(zoomEdit, SIGNAL(returnPressed()), this, SLOT(zoom()));

	zoomChoice = new QToolButton(this);
    QMenu * pZoom = new QMenu();
	
    /* FIXME
	pZoom->insertItem( "&10%", this, SLOT(zoomWithValue(int)), 0, 10);
	pZoom->insertItem( "&25%", this, SLOT(zoomWithValue(int)), 0, 25);
	pZoom->insertItem( "&50%", this, SLOT(zoomWithValue(int)), 0, 50);
	pZoom->insertItem( "&75%", this, SLOT(zoomWithValue(int)), 0, 75);
	pZoom->insertItem( "&100%", this, SLOT(zoomWithValue(int)), 0, 100);
	pZoom->insertItem( "&150%", this, SLOT(zoomWithValue(int)), 0, 150);
	pZoom->insertItem( "&250%", this, SLOT(zoomWithValue(int)), 0, 250);
	pZoom->insertItem( "&500%", this, SLOT(zoomWithValue(int)), 0, 500);
	pZoom->insertSeparator();
	pZoom->insertItem( "&Fit to Screen", pParent, SLOT(fitToScreen()));
	pZoom->insertItem( "&Fit Width", pParent, SLOT(fitToWidth()));
	pZoom->insertItem( "&Fit Height", pParent, SLOT(fitToHeight()));
    */
    zoomChoice->setMenu(pZoom);
    //zoomChoice->setPopupDelay(0);
	zoomChoice->setFixedSize(16, 30);

    zoomInBtn = new QToolButton( this);
    zoomInBtn->setIcon(QIcon(":/icons/IconPlus.png"));
    zoomInBtn->setText("Zoom In");
    zoomInBtn->setParent(pParent);
    connect(zoomInBtn, SIGNAL(triggered(QAction *)), this, SLOT( zoomIn() )); // TODO check
	zoomInBtn->setFixedSize(25, 30);
	this->setDisabled(true);
	return 0;
}

void TMOGUIZoomTool::SetState()
{
	QString s;
	zoomEdit->setText(s.setNum(pImage->pImage->GetRatio())+'%');
}

int TMOGUIZoomTool::SetWindows(QMdiArea* w)
{
    if (w && !w->subWindowList().isEmpty())
		this->setDisabled(false);
	else
		this->setDisabled(true);
	return 0;
}

void TMOGUIZoomTool::windowChanged(TMOGUIImage* pWindow)
{
	QString s;
	pImage = pWindow;
	if(pImage)
		zoomEdit->setText(s.setNum(pImage->pImage->GetRatio())+'%');
}

void TMOGUIZoomTool::zoomWithValue(int value)
{
	pImage->pImage->Zoom((double)value / 100);	
	pImage->pImage->SetSize();	
	SetState();
	pImage->SetImageZoomLabel();
}

void TMOGUIZoomTool::zoom()
{	
	if(!pImage) return;
	QString s;
	QString str = zoomEdit->text();
	if(str.contains('%') != 0)	
		str = str.section('%', 0);
	bool ok;
	int value = str.toInt(&ok);
	if(!ok)
	{
		zoomEdit->setText(s.setNum(pImage->pImage->GetRatio())+'%');
		return;
	}	
	if(10 > value)
	{
		zoomEdit->setText("10%");
		value = 10;
	}
	if(value > 1000)
	{
		value = 1000;
		zoomEdit->setText("1000%");
	}	
	pImage->pImage->Zoom((double)value / 100);
	zoomEdit->setText(s.setNum(pImage->pImage->GetRatio())+'%');
	pImage->pImage->SetSize();
	pImage->SetImageZoomLabel();
}
