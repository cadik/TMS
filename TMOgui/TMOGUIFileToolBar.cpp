// TMOGUIFileToolBar.cpp: implementation of the TMOGUIFileToolBar class.
//
//////////////////////////////////////////////////////////////////////
	
#include "TMOGUIFileToolBar.h"
#include "TMOGUIResource.h"
#include <qtoolbutton.h>
#include <qlabel.h>
#include <qworkspace.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIFileToolBar::TMOGUIFileToolBar(QWidget * parent, const char * name):Q3ToolBar((Q3MainWindow*)parent, name)
{
	pParent = parent;
	Create();
}

TMOGUIFileToolBar::~TMOGUIFileToolBar()
{
}

int TMOGUIFileToolBar::Create()
{
	this->setLabel( "File Operations" );

	fileNewBtn = new QToolButton(*TMOResource::pResource->IconNew->pixmap(), "New file", QString::null, pParent, SLOT( newFile() ), this );
	fileOpenBtn = new QToolButton(*TMOResource::pResource->IconOpen->pixmap(), "Open file", QString::null, pParent, SLOT( openFile() ), this );
	fileSaveBtn = new QToolButton(*TMOResource::pResource->IconSave->pixmap(), "Save file", QString::null, pParent, SLOT( saveFile() ), this );
	fileSaveAllBtn = new QToolButton(*TMOResource::pResource->IconSaveAll->pixmap(), "Save all files", QString::null, pParent, SLOT( saveallFile() ), this );
	filePrintBtn = new QToolButton(*TMOResource::pResource->IconPrint->pixmap(), "Print file", QString::null, pParent, SLOT( printFile() ), this );
	
	fileSaveAllBtn->setDisabled(true);
	fileSaveBtn->setDisabled(true);
	filePrintBtn->setDisabled(true);
	return 0;
}

int TMOGUIFileToolBar::SetWindows(QWorkspace* w)
{
	if (w && !w->windowList().isEmpty())
	{
		fileSaveAllBtn->setDisabled(false);
		fileSaveBtn->setDisabled(false);
		filePrintBtn->setDisabled(false);
	}
	else
	{
		fileSaveAllBtn->setDisabled(true);
		fileSaveBtn->setDisabled(true);
		filePrintBtn->setDisabled(true);
	}
	return 0;
}

