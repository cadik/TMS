// TMOGUIFileToolBar.cpp: implementation of the TMOGUIFileToolBar class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIFileToolBar.h"
#include "TMOGUIResource.h"
#include "switch.h"
#include <QToolButton>
#include <QLabel>
#include <QMdiArea>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIFileToolBar::TMOGUIFileToolBar(QWidget *parent, const char *name) : QToolBar(name, parent)
{
    pParent = parent;
    Create();
}

TMOGUIFileToolBar::~TMOGUIFileToolBar()
{
}

int TMOGUIFileToolBar::Create()
{
    this->setWindowTitle("File Operations");

    fileNewBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconNew.png"), "New file", QString(), pParent, SLOT(newFile()), this, fileOpenBtn);
    fileOpenBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconOpen.png"), "Open file", QString(), pParent, SLOT(openFile()), this, fileSaveBtn);
    fileSaveBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconSave.png"), "Save file", QString(), pParent, SLOT(saveFile()), this, fileSaveAllBtn);
    fileSaveAllBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconSaveAll.png"), "Save all files", QString(), pParent, SLOT(saveallFile()), this, filePrintBtn);
    filePrintBtn = TMOGUIFileToolBar::addButton(QIcon(":/resources/icons/IconPrint.png"), "Print file", QString(), pParent, SLOT(printFile()), this, fileWorkspace);
    switchWorkspace = new Switch("Advanced");
    fileWorkspace = this->insertWidget(nullptr, switchWorkspace);
    connect(switchWorkspace, SIGNAL(stateChanged(int)), pParent, SLOT(changeWorkspace(int))); // FIXME action, position, initial state

    fileSaveAllBtn->setDisabled(true);
    fileSaveBtn->setDisabled(true);
    filePrintBtn->setDisabled(true);
    return 0;
}

int TMOGUIFileToolBar::SetWindows(QMdiArea *w)
{
    if (w && !w->subWindowList().isEmpty()) // windowList
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
    QToolButton *temp = new QToolButton(parent);
    temp->setToolButtonStyle(Qt::ToolButtonStyle::ToolButtonIconOnly);
    temp->setIcon(s);
    temp->setToolTip(textLabel);
    temp->setStatusTip(grouptext);
    QAction *act = this->insertWidget(before, temp);
    connect(temp, SIGNAL(clicked()), receiver, slot);
    return act;
}
