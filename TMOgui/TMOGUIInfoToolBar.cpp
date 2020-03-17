// TMOGUIInfoToolBar.cpp: implementation of the TMOGUIInfoToolBar class.
//
//////////////////////////////////////////////////////////////////////
	
#include "TMOGUIInfoToolBar.h"
#include "TMOGUIResource.h"
#include <qtoolbutton.h>
#include <qlabel.h>
#include <QMdiArea>
#include <QMenu>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIInfoToolBar::TMOGUIInfoToolBar(QWidget * parent, const char * name):QToolBar( parent)
{
	pParent = parent;
	Create();
}

TMOGUIInfoToolBar::~TMOGUIInfoToolBar()
{
}

int TMOGUIInfoToolBar::Create()
{
    this->setWindowTitle( "Local tools" );
    toolSettingBtn = new QToolButton(this);
    toolSettingBtn->setFixedSize(11, 30);
    toolSettingBtn->setIcon(QIcon(":/resources/icons/IconArrow.png"));
    toolSettingBtn->setText("Local info selection tool setting");
    connect(toolSettingBtn, SIGNAL(clicked()), pParent, SLOT(showToolSetting()));


    toolBtn = new QToolButton(this);
    toolBtn->setIcon(QIcon(":/resources/icons/IconTool.png"));
    toolBtn->setText("Local info selection tool");
    toolBtn->setCheckable(true);
    connect(toolBtn, SIGNAL(toggled(bool)), pParent, SLOT(activateInfoTool(bool)));


    toolBtnAct = this->insertWidget(nullptr, toolBtn);
    toolSettingBtnAct = this->insertWidget(toolBtnAct, toolSettingBtn);
	


	this->setDisabled(true);
	return 0;
}

int TMOGUIInfoToolBar::SetWindows(QMdiArea* w)
{
    if (w && !w->subWindowList().isEmpty())
		this->setDisabled(false);
	else
	{
        toolBtn->setChecked(false);
		this->setDisabled(true);
	}	
	return 0;
}

bool TMOGUIInfoToolBar::IsActivated() const
{ 
    return toolBtn->isChecked();
}

