// TMOGUIMenu.cpp: implementation of the TMOGUIMenu class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIMenu.h"
#include "TMOGUIImage.h"
#include "TMOGUIResource.h"
#include "TMOGUIStyle.h"
#include <QMenu>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <qfile.h>
#include <QTextStream>
#include <qlabel.h>
#include <qmessagebox.h>
#include <qpushbutton.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIMenu::TMOGUIMenu(QWidget *parent, const char *name) : QMenuBar(parent)
{

    pParent = parent;
    pRecent = nullptr;
    this->setObjectName(name);

    qDeleteAll(listRecent); // listRecent.setAutoDelete(true);
    Create();
}

TMOGUIMenu::~TMOGUIMenu()
{
    SaveRecent();
}

int TMOGUIMenu::Create()
{
    pRecent = new QMenu("Recent", this);
    pComponent = new QMenu("Component", this);
    pFile = new QMenu("File", this);
    pEdit = new QMenu("Edit", this);
    pView = new QMenu("View", this);
    pCommand = new QMenu("Command", this);
    pWindows = new QMenu("Window", this);
    pHelpIt = new QMenu("Help", this);

    QPalette p = palette();
    p.setColor(QPalette::Base, Qt::white);

    setAutoFillBackground(true);
    setPalette(p);
    pFile->setPalette(p);
    pEdit->setPalette(p);
    pView->setPalette(p);

    this->addMenu(pFile);
    this->addMenu(pEdit);
    this->addMenu(pView);
    this->addMenu(pCommand);
    this->addMenu(pWindows);
    this->addMenu(pHelpIt);

    LoadRecent();

    pFileAct.insert(1, pFile->addAction(QIcon(":/resources/icons/IconNew.png"), "&New", pParent, SLOT(newFile()), Qt::CTRL + Qt::Key_N));    //, 1);
    pFileAct.insert(2, pFile->addAction(QIcon(":/resources/icons/IconOpen.png"), "&Open", pParent, SLOT(openFile()), Qt::CTRL + Qt::Key_O)); //, 2 );
    pFileAct.insert(3, pFile->addAction("&Close", pParent, SLOT(closeFile()), 0));                                                           //, 3);
    pFile->addSeparator();
    pFileAct.insert(4, pFile->addAction(QIcon(":/resources/icons/IconSave.png"), "&Save", pParent, SLOT(saveFile()), Qt::CTRL + Qt::Key_S)); //, 4 );
    pFileAct.insert(5, pFile->addAction(QIcon(":/resources/icons/IconSave.png"), "Save &As...", pParent, SLOT(saveasFile()), 0));            //, 5 );
    pFileAct.insert(6, pFile->addAction(QIcon(":/resources/icons/IconSaveAll.png"), "Save A&ll", pParent, SLOT(saveallFile()), 0));          //, 6 );
    pFile->addSeparator();
    pFileAct.insert(7, pFile->addAction("Page Set&up...", pParent, SLOT(pageFile()), 0));                                                          //, 7 );
    pFileAct.insert(8, pFile->addAction(QIcon(":/resources/icons/IconPrint.png"), "&Print...", pParent, SLOT(printFile()), Qt::CTRL + Qt::Key_P)); //, 8 );
    pFileAct.insert(9, pFile->addSeparator());
    pFileAct.insert(10, pFile->addAction("E&xit", pParent, SLOT(exitFile()), 0)); //, 10);
    //addAction("&File", this, SLOT(pFile));//, 1);
    pFileAct.insert(9, pFile->insertMenu(pFileAct.take(9), pRecent));

    pEditAct.insert(1, pEdit->addAction(QIcon(":/resources/icons/IconUndo.png"), "&Undo", pParent, SLOT(undoEdit()), Qt::CTRL + Qt::Key_Z)); //, 1);
    //pEditAct.insert(2, pEdit->addAction( "&Redo", pParent, SLOT(redoEdit()), Qt::CTRL+Qt::Key_Y));//, 2);
    pEdit->addSeparator();
    //addAction("&Edit", this ,SLOT(pEdit));//, 2);
    Disable(2, 1);
    //	Disable(2, 2);

    pViewAct.insert(1, pView->addAction("&Output log", pParent, SLOT(viewInfo()), Qt::ALT + Qt::Key_0));     //, 1);
    pViewAct.insert(2, pView->addAction("&Tools", pParent, SLOT(viewRight()), Qt::ALT + Qt::Key_1));         //, 2);
    pViewAct.insert(3, pView->addAction("&Histogram", pParent, SLOT(viewHistogram()), Qt::ALT + Qt::Key_2)); //, 3);
    pView->addSeparator();
    pViewAct.insert(4, pView->addAction("Zoom &In", pParent, SLOT(zoomIn()), Qt::CTRL + Qt::Key_Plus));                                                //, 4);
    pViewAct.insert(5, pView->addAction("Zoom &Out", pParent, SLOT(zoomOut()), Qt::CTRL + Qt::Key_Minus));                                             //, 5);
    pViewAct.insert(6, pView->addAction(QIcon(":/resources/icons/IconFit.png"), "Fit To &Screen", pParent, SLOT(fitToScreen()), Qt::ALT + Qt::Key_3)); //, 6);
    pViewAct.insert(7, pView->addAction(QIcon(":/resources/icons/IconWidth.png"), "Fit &Width", pParent, SLOT(fitToWidth()), Qt::ALT + Qt::Key_4));    //, 7);
    pViewAct.insert(8, pView->addAction(QIcon(":/resources/icons/IconHeight.png"), "Fit &Height", pParent, SLOT(fitToHeight()), Qt::ALT + Qt::Key_5)); //, 8);
    //addAction("&View", this, SLOT(pView));//, 3);
    pViewAct.value(1)->setCheckable(true);
    pViewAct.value(2)->setCheckable(true);
    //pViewAct.value(3)->setCheckable(true);

    pComponentAct.insert(0, pComponent->addAction("&Red", pParent, SLOT(extractRed())));     //, 0);
    pComponentAct.insert(1, pComponent->addAction("&Green", pParent, SLOT(extractGreen()))); //, 1);
    pComponentAct.insert(2, pComponent->addAction("&Blue", pParent, SLOT(extractBlue())));   //, 2);

    pCommandAct.insert(1, pCommand->addAction(QIcon(":/resources/icons/IconDuplicate.png"), "&Duplicate Image", pParent, SLOT(duplicateCommand()), Qt::CTRL + Qt::Key_D)); //, 1);
    pCommandAct.insert(2, pCommand->addAction(QIcon(":/resources/icons/IconSize.png"), "Change Image &Size", pParent, SLOT(sizeCommand()), 0));                            //, 2);
    pCommand->addSeparator();
    pCommandAct.insert(3, pCommand->addAction("&Extract Luminance", pParent, SLOT(extractLumCommand()), Qt::CTRL + Qt::Key_L)); //, 3);
    pCommandAct.insert(4, pCommand->insertMenu(pCommandAct.take(5), pComponent));
    pCommandAct.insert(5, pCommand->addAction("&Merge components", pParent, SLOT(mergeCommand()), 0)); //, 5);
    pCommand->addSeparator();
    pCommandAct.insert(6, pCommand->addAction("Arithmetical &Operation", pParent, SLOT(operationCommand()), 0)); //, 6);
    //addAction("&Command", pCommand, SLOT(pCommand));//, 5);

    //pHelpItAct.insert(1, pHelpIt->addAction( "&Help", pParent, SLOT(showHelp()), Qt::Key_F1));//, 0);//, 1);
    //pHelpIt->addSeparator();
    pHelpItAct.insert(3, pHelpIt->addAction("&About", this, SLOT(about()), 0)); //, 3);
                                                                                //pHelpIt->addSeparator();
                                                                                //pHelpItAct.insert(5, pHelpIt->addAction( QIcon(":/resources/icons/IconThis.png"), "&What's This?", pParent,  SLOT(whatsThis()), Qt::SHIFT+Qt::Key_F1));//, 5);
                                                                                //addAction("&Help", pHelpIt, SLOT(pHelpIt));//, 5);

    pImage = 0;
    SetWindows(0);
    /*
    connect (pWindows, SIGNAL(triggered(QAction*)), pParent, SLOT(activateWindow(int)));
    connect (pRecent, SIGNAL(triggered(QAction*)), pParent, SLOT(openFile(int)));
    connect (pComponent, SIGNAL(triggered(QAction*)), pParent, SLOT(extractComCommand(int)));
    */

    return 0;
}

int TMOGUIMenu::SetWindows(QMdiArea *w)
{
    QList<QMdiSubWindow *> wl; //QWidgetList wl;
    QMdiSubWindow *widget;
    QString s;
    int number = 64;

    if (w)
        wl = w->subWindowList();
    //removeItem(4);
    pWindowsAct.clear();
    pWindows->clear();

    if (w)
    {
        pWindowsAct.insert(1, pWindows->addAction(QIcon(":/resources/icons/IconCascade.png"), "&Cascade", w, &QMdiArea::cascadeSubWindows, 0)); //, 1);
        pWindowsAct.insert(2, pWindows->addAction(QIcon(":/resources/icons/IconTile.png"), "&Tile", w, &QMdiArea::tileSubWindows, 0));          //, 2);
    }
    pWindowsAct.insert(3, pWindows->addAction("Close &All", pParent, SLOT(closeallWindow()), 0)); //, 3);
    pWindows->addSeparator();

    int i = 4;
    foreach (widget, wl)
    {
        TMOGUIImage *pImg = (TMOGUIImage *)widget->widget();
        if (pImg)
        {
            s = pImage->objectName();
        }
        else
        {
            s = "unknown";
            //s = widget->widget()->objectName();
        }
        s = TMOGUIImage::GetName(s);
        if (pImg->bPreview)
            s = "Preview - " + s;

        QAction *act = pWindows->addAction(s, pParent, [=]()
                                           { emit activateWindowAction(widget->windowTitle()); });
        act->setCheckable(true);
        pWindowsAct.insert(i, act); //, 3);

        number++;
        if (pImage)
        {
            QString h = pImage->windowTitle();
            //if(pImage->bPreview) h = "Preview - " + h;
            if (s == h)
                act->setChecked(true);
        }

        i++;
    }

    if (number == 64)
    {
        Disable(1, 3);
        Disable(1, 4);
        Disable(1, 5);
        Disable(1, 6);
        Disable(1, 7);
        Disable(1, 8);
        if (listRecent.count() == 0)
            Disable(1, 9);
        Disable(2, 1);
        //Disable(2,2);
        //Disable(2,3);
        Disable(3, 3);
        Disable(3, 4);
        Disable(3, 5);
        Disable(3, 6);
        Disable(3, 7);
        Disable(3, 8);
        Disable(4, 1);
        Disable(4, 2);
        Disable(4, 3);
        Disable(4, 4);
        Disable(4, 5);
        Disable(4, 6);
        //		Disable(5,1);
        //		Disable(5,2);
        //		Disable(5,3);
    }
    else
    {
        Enable(1, 3);
        Enable(1, 4);
        Enable(1, 5);
        Enable(1, 6);
        Enable(1, 7);
        Enable(1, 8);
        if (listRecent.count())
            Enable(1, 9);
        Enable(3, 3);
        Enable(3, 4);
        Enable(3, 5);
        Enable(3, 6);
        Enable(3, 7);
        Enable(3, 8);
        Enable(4, 1);
        Enable(4, 2);
        Enable(4, 3);
        Enable(4, 4);
        Enable(4, 5);
        Enable(4, 6);
        Enable(5, 1);
        Enable(5, 2);
        Enable(5, 3);
    }
    this->addMenu(pWindows);
    this->addMenu(pHelpIt);
    //this->addAction("&Window", pWindows, SLOT(pWindows));//addAction("&Window", pWindows, 4, 4);
    return 0;
}

void TMOGUIMenu::windowChanged(TMOGUIImage *pImg)
{
    if (!pImg || pImg->windowTitle().isNull())
        return;
    pImage = pImg;
    int index = 4;
    while (pWindowsAct.contains(index))
    {
        pWindowsAct.value(index)->setChecked(false);
        //pWindows->setItemChecked(id, false);
        QString s = pWindowsAct.value(index)->text();
        //QString s = pWindows->text(id);
        QString h = pImage->windowTitle();
        //if(pImage->bPreview) h = "Preview - " + h;
        if (s == h)
            pWindowsAct.value(index)->setChecked(true); //setItemChecked(id, true);

        index++;
    }
    if (!pImg->bPreview)
        pViewAct.value(2)->setChecked(pImage->pToolsButton->isChecked());
    //pView->setItemChecked(pView->idAt(2), pImage->pToolsButton->isOn());
}

int TMOGUIMenu::Enable(int menu, int item)
{
    switch (menu)
    {
    case 1:
        pFileAct.value(item)->setEnabled(true);
        break;
    case 2:
        pEditAct.value(item)->setEnabled(true);
        break;
    case 3:
        pViewAct.value(item)->setEnabled(true);
        break;
    case 4:
        pCommandAct.value(item)->setEnabled(true);
        break;
    case 5:
        pWindowsAct.value(item)->setEnabled(true);
        break;
    default:
        return 1;
    }
    return 0;
}

int TMOGUIMenu::SetHidden(int menu, int item, bool hidden)
{
    switch (menu)
    {
    case 1:
        pFileAct.value(item)->setVisible(!hidden);
        break;
    case 2:
        pEditAct.value(item)->setVisible(!hidden);
        break;
    case 3:
        pViewAct.value(item)->setVisible(!hidden);
        break;
    case 4:
        pCommandAct.value(item)->setVisible(!hidden);
        break;
    case 5:
        pWindowsAct.value(item)->setVisible(!hidden);
        break;
    default:
        return 1;
    }
    return 0;
}

int TMOGUIMenu::GetChecked(int menu, int item)
{
    int res = -1;
    switch (menu)
    {
    case 1:
        res = pFileAct.value(item)->isChecked();
        break;
    case 2:
        res = pEditAct.value(item)->isChecked();
        break;
    case 3:
        res = pViewAct.value(item)->isChecked();
        break;
    case 4:
        res = pCommandAct.value(item)->isChecked();
        break;
    case 5:
        res = pWindowsAct.value(item)->isChecked();
        break;
    default:
        return -1;
    }
    return res;
}

int TMOGUIMenu::SetChecked(int menu, int item, bool checked)
{
    switch (menu)
    {
    case 1:
        pFileAct.value(item)->setChecked(checked);
        break;
    case 2:
        pEditAct.value(item)->setChecked(checked);
        break;
    case 3:
        pViewAct.value(item)->setChecked(checked);
        break;
    case 4:
        pCommandAct.value(item)->setChecked(checked);
        break;
    case 5:
        pWindowsAct.value(item)->setChecked(checked);
        break;
    default:
        return 1;
    }
    return 0;
}

int TMOGUIMenu::Disable(int menu, int item)
{
    switch (menu)
    {
    case 1:
        pFileAct.value(item)->setEnabled(false);
        break;
    case 2:
        pEditAct.value(item)->setEnabled(false);
        break;
    case 3:
        pViewAct.value(item)->setEnabled(false);
        break;
    case 4:
        pCommandAct.value(item)->setEnabled(false);
        break;
    case 5:
        pWindowsAct.value(item)->setEnabled(false);
        break;
    default:
        return 1;
    }
    return 0;
}

int TMOGUIMenu::LoadRecent()
{
    QFile f("recent.dat");

    if (f.open(QIODevice::ReadOnly))
    {
        QTextStream t(&f);
        QString s;
        listRecent.clear();
        int n = 1;
        while (!t.atEnd()) // .eof()
        {
            s = t.readLine();
            listRecent.append(new QString(s));
        }
        f.close();
        SetRecent();
        return 0;
    }
    return -1;
}

int TMOGUIMenu::SaveRecent()
{
    QFile f("recent.dat");

    if (f.open(QIODevice::WriteOnly))
    {
        QTextStream t(&f);

        for (QString *temp : listRecent)
        {
            t << *temp << '\n';
        }
        f.close();
        return 0;
    }
    return -1;
}

int TMOGUIMenu::AddRecent(QString &s)
{

    for (QString *temp : listRecent)
    {
        if (*temp == s)
        {
            listRecent.removeAt(listRecent.indexOf(temp));
            break;
        }
    }
    if (listRecent.count() > 4)
        listRecent.removeLast();
    listRecent.prepend(new QString(s));

    return 0;
}

int TMOGUIMenu::SetRecent()
{

    int number = 128;

    pRecent->clear();
    for (QString *s : listRecent)
    {
        QAction *action = pRecent->addAction(*s, pParent, [=]()
                                             { emit openFile(*s); });

        number++;
    }

    return 0;
}

QString TMOGUIMenu::GetRecent(int ID)
{

    int number = 128;

    for (QString *s : listRecent)
    {
        if (number == ID)
            return *s;
        number++;
    }
    return QString("");
}

void TMOGUIMenu::about()
{
    QMessageBox::about(this, "Tone Mapping Studio",
                       "<p><b>(c) 2004-2024 <br />Ondrej Hajdok, Martin Cadik, Jan Jedlicka</b></p>"
                       "<p>"
                       "<a href=\"mailto:cadikm@centrum.cz\">cadikm@centrum.cz</a></p>");
}
