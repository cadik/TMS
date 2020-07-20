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
    zoomOutBtn->setIcon(QIcon(":/resources/icons/IconMinus.png"));
    zoomOutBtn->setText("Zoom Out");
    zoomOutBtn->setParent( pParent);

    connect(zoomOutBtn, SIGNAL(clicked()), pParent, SLOT( zoomOut() ));
    //zoomOutBtn->setFixedSize(25, 30);

	zoomEdit = new QLineEdit("100%", this);
    zoomEdit->setFixedWidth(100);
	zoomEdit->setAlignment(Qt::AlignCenter);
    //connect(zoomEdit, SIGNAL(returnPressed()), this, SLOT(zoom()));
    connect(zoomEdit, SIGNAL(editingFinished()), this, SLOT(zoom()));

    zoomChoice = new QToolButton(this);
    QMenu * pZoom = new QMenu();
	
    signalMapper = new QSignalMapper(this);
    connect(signalMapper, SIGNAL(mapped(int)), SLOT(zoomWithValue(int)));
    QAction* action = pZoom->addAction( "&10%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 10);
    action = pZoom->addAction( "&25%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 25);
    action = pZoom->addAction( "&50%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 50);
    action = pZoom->addAction( "&75%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 75);
    action = pZoom->addAction( "&100%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 100);
    action = pZoom->addAction( "&150%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 150);
    action = pZoom->addAction( "&250%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 250);
    action = pZoom->addAction( "&500%", signalMapper, SLOT(map()));
    signalMapper->setMapping(action, 500);
    action = pZoom->addSeparator();
    action = pZoom->addAction( "&Fit to Screen", pParent, SLOT(fitToScreen()));
    action = pZoom->addAction( "&Fit Width", pParent, SLOT(fitToWidth()));
    action = pZoom->addAction( "&Fit Height", pParent, SLOT(fitToHeight()));

    zoomChoice->setMenu(pZoom);
    //zoomChoice->setPopupMode(QToolButton::MenuButtonPopup);
    zoomChoice->setArrowType(Qt::NoArrow);
    zoomChoice->setToolButtonStyle(Qt::ToolButtonStyle::ToolButtonIconOnly);
    zoomChoice->setIcon(QIcon(":/resources/icons/IconArrow.png"));
    connect(zoomChoice, SIGNAL(clicked()), zoomChoice, SLOT(showMenu()));
    //zoomChoice->setPopupDelay(0);
    zoomChoice->setFixedSize(16, 30);


    zoomInBtn = new QToolButton(this);
    zoomInBtn->setIcon(QIcon(":/resources/icons/IconPlus.png"));
    zoomInBtn->setText("Zoom In");
    zoomInBtn->setParent(pParent);
    connect(zoomInBtn, SIGNAL(clicked()), pParent, SLOT( zoomIn() ));

    //zoomInBtn->setFixedSize(25, 30);
    zoomInBtnAct = this->insertWidget(nullptr, zoomInBtn);
    zoomChoiceAct = this->insertWidget(zoomInBtnAct, zoomChoice);
    zoomEditAct = this->insertWidget(zoomChoiceAct, zoomEdit);
    zoomOutBtnAct = this->insertWidget(zoomEditAct, zoomOutBtn);

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

    str = str.remove('%');

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


