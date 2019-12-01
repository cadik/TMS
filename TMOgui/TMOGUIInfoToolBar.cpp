// TMOGUIInfoToolBar.cpp: implementation of the TMOGUIInfoToolBar class.
//
//////////////////////////////////////////////////////////////////////
	
#include "TMOGUIInfoToolBar.h"
#include "TMOGUIResource.h"
#include <qtoolbutton.h>
#include <qlabel.h>
#include <qworkspace.h>
#include <q3popupmenu.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIInfoToolBar::TMOGUIInfoToolBar(QWidget * parent, const char * name):Q3ToolBar((Q3MainWindow*)parent, name)
{
	pParent = parent;
	Create();
}

TMOGUIInfoToolBar::~TMOGUIInfoToolBar()
{
}

int TMOGUIInfoToolBar::Create()
{
	this->setLabel( "Local tools" );
	toolBtn = new QToolButton(this, "Local info selection tool");
	toolBtn->setIconSet(*TMOResource::pResource->IconTool->pixmap());
	toolBtn->setTextLabel("Local info selection tool");
	toolBtn->setToggleButton(true);
	connect (toolBtn, SIGNAL(toggled(bool)), pParent, SLOT(activateInfoTool(bool)));
	
	toolSettingBtn = new QToolButton(this, "Lis tool setting");
	toolSettingBtn->setFixedSize(11, 30);
	toolSettingBtn->setIconSet(*TMOResource::pResource->IconArrow->pixmap());
	toolSettingBtn->setTextLabel("Local info selection tool setting");	
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
		toolBtn->setOn(false);
		this->setDisabled(true);
	}	
	return 0;
}

bool TMOGUIInfoToolBar::IsActivated() const
{ 
	return toolBtn->isOn(); 
}