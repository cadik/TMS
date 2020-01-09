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

TMOGUIFileToolBar::TMOGUIFileToolBar(QWidget * parent, const char * name):QToolBar(name, parent)
{
	pParent = parent;
	Create();
}

TMOGUIFileToolBar::~TMOGUIFileToolBar()
{
}

int TMOGUIFileToolBar::Create()
{
    this->setWindowTitle( "File Operations" );

    fileNewBtn = TMOGUIFileToolBar::addButton(*TMOResource::pResource->IconNew->pixmap(), "New file", QString::null, pParent, SLOT( newFile() ), this );
    fileOpenBtn = TMOGUIFileToolBar::addButton(*TMOResource::pResource->IconOpen->pixmap(), "Open file", QString::null, pParent, SLOT( openFile() ), this );
    fileSaveBtn = TMOGUIFileToolBar::addButton(*TMOResource::pResource->IconSave->pixmap(), "Save file", QString::null, pParent, SLOT( saveFile() ), this );
    fileSaveAllBtn = TMOGUIFileToolBar::addButton(*TMOResource::pResource->IconSaveAll->pixmap(), "Save all files", QString::null, pParent, SLOT( saveallFile() ), this );
    filePrintBtn = TMOGUIFileToolBar::addButton(*TMOResource::pResource->IconPrint->pixmap(), "Print file", QString::null, pParent, SLOT( printFile() ), this );
	
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

QToolButton *TMOGUIFileToolBar::addButton(const QIcon &s, const QString &textLabel,
                                         const QString &grouptext, QObject *receiver,
                                         const char *slot, QWidget *parent, const char *name)
{
    QToolButton* temp = new QToolButton(parent);
    temp->setIcon(s);
    temp->setToolTip(textLabel); // TODO check
    temp->setStatusTip(grouptext);
    connect (temp, SIGNAL(clicked()), receiver, slot); // TODO check
    return temp;
}

