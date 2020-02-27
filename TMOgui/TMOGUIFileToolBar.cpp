// TMOGUIFileToolBar.cpp: implementation of the TMOGUIFileToolBar class.
//
//////////////////////////////////////////////////////////////////////
	
#include "TMOGUIFileToolBar.h"
#include "TMOGUIResource.h"
#include <QToolButton>
#include <QLabel>
#include <QMdiArea>


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

    fileNewBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconNew.png"), "New file", QString(), pParent, SLOT( newFile() ), this, nullptr );
    fileOpenBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconOpen.png"), "Open file", QString(), pParent, SLOT( openFile() ), this, fileNewBtn );
    fileSaveBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconSave.png"), "Save file", QString(), pParent, SLOT( saveFile() ), this, fileOpenBtn );
    fileSaveAllBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconSaveAll.png"), "Save all files", QString(), pParent, SLOT( saveallFile() ), this, fileSaveBtn );
    filePrintBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconPrint.png"), "Print file", QString(), pParent, SLOT( printFile() ), this, fileSaveAllBtn );

	fileSaveAllBtn->setDisabled(true);
	fileSaveBtn->setDisabled(true);
	filePrintBtn->setDisabled(true);
	return 0;
}

int TMOGUIFileToolBar::SetWindows(QMdiArea* w)
{
    if (w && !w->subWindowList().isEmpty())  // windowList
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

QAction *TMOGUIFileToolBar::addButton(const QIcon &s, const QString &textLabel,
                                         const QString &grouptext, QObject *receiver,
                                         const char *slot, QWidget *parent, QAction *before)
{
    QToolButton* temp = new QToolButton(parent);
    temp->setToolButtonStyle(Qt::ToolButtonStyle::ToolButtonIconOnly); // TODO FIXME
    temp->setIcon(s);
    temp->setToolTip(textLabel); // TODO check
    temp->setStatusTip(grouptext);
    QAction *act = this->insertWidget(before, temp);
    connect(temp, SIGNAL(clicked()), receiver, slot); // TODO check
    return act;
}

