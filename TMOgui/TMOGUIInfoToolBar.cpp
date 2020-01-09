// TMOGUIInfoToolBar.cpp: implementation of the TMOGUIInfoToolBar class.
//
//////////////////////////////////////////////////////////////////////
	
#include "TMOGUIInfoToolBar.h"
#include "TMOGUIResource.h"
#include <qtoolbutton.h>
#include <qlabel.h>
#include <qworkspace.h>
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
    toolBtn = new QToolButton(this);//, "Local info selection tool");
    toolBtn->setIcon(*TMOResource::pResource->IconTool->pixmap());
    toolBtn->setText("Local info selection tool");
    toolBtn->setCheckable(true);
	connect (toolBtn, SIGNAL(toggled(bool)), pParent, SLOT(activateInfoTool(bool)));
	
    toolSettingBtn = new QToolButton(this);//, "Lis tool setting");
	toolSettingBtn->setFixedSize(11, 30);
    toolSettingBtn->setIcon(*TMOResource::pResource->IconArrow->pixmap());
    toolSettingBtn->setText("Local info selection tool setting");
	connect (toolSettingBtn, SIGNAL(clicked()), pParent, SLOT(showToolSetting()));

	this->setDisabled(true);
	return 0;
}

int TMOGUIInfoToolBar::SetWindows(QWorkspace* w)
{
	if (w && !w->windowList().isEmpty())
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
