// TMOGUIMenu.cpp: implementation of the TMOGUIMenu class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIMenu.h"
#include "TMOGUIImage.h"
#include "TMOGUIResource.h"
#include <QMenu>
#include <qworkspace.h>
#include <qfile.h>
#include <q3textstream.h>
#include <qlabel.h>
#include <qmessagebox.h>
#include <qpushbutton.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIMenu::TMOGUIMenu(QWidget * parent, const char * name):QMenuBar(parent)
{
	pParent = parent;
	pRecent = 0;

    qDeleteAll(listRecent); // TODO listRecent.setAutoDelete(true);
	Create();
}

TMOGUIMenu::~TMOGUIMenu()
{
	SaveRecent();
}

int TMOGUIMenu::Create()
{
    pRecent = new QMenu(this);//, "Recent");
    pComponent = new QMenu(this);//, "Component");
    pFile = new QMenu(this);//, "File");
    pEdit = new QMenu(this);//, "Edit");
    pView = new QMenu(this);//, "View");
    pCommand = new QMenu(this);//, "Command");
    pWindows = new QMenu(this);//, "Window");
    pHelpIt = new QMenu(this);//, "Help");

	LoadRecent();

    pFile->addAction( *TMOResource::pResource->IconNew->pixmap(), "&New",  pParent, SLOT(newFile()), Qt::CTRL+Qt::Key_N );//, 1);
    pFile->addAction( *TMOResource::pResource->IconOpen->pixmap(), "&Open", pParent, SLOT(openFile()), Qt::CTRL+Qt::Key_O);//, 2 );
    pFile->addAction( "&Close", pParent, SLOT(closeFile()), 0);//, 3);
    pFile->addSeparator();
    pFile->addAction( *TMOResource::pResource->IconSave->pixmap(), "&Save",  pParent, SLOT(saveFile()), Qt::CTRL+Qt::Key_S);//, 4 );
    pFile->addAction( *TMOResource::pResource->IconSave->pixmap(), "Save &As...", pParent, SLOT(saveasFile()), 0);//, 5 );
    pFile->addAction( *TMOResource::pResource->IconSaveAll->pixmap(), "Save A&ll", pParent, SLOT(saveallFile()), 0);//, 6 );
    pFile->addSeparator();
    pFile->addAction( "Page Set&up...",  pParent, SLOT(pageFile()), 0);//, 7 );
    pFile->addAction( *TMOResource::pResource->IconPrint->pixmap(), "&Print...", pParent, SLOT(printFile()), Qt::CTRL+Qt::Key_P);//, 8 );
    pFile->addSeparator();
    pFile->addAction( "Recent &Files", pFile, SLOT(pRecent), 0);//, 9 );
    pFile->addSeparator();
    pFile->addAction( "E&xit",  pParent, SLOT(exitFile()), 0);//, 10);
    addAction("&File", this, SLOT(pFile));//, 1);

    pEdit->addAction( *TMOResource::pResource->IconUndo->pixmap(), "&Undo", pParent, SLOT(undoEdit()), Qt::CTRL+Qt::Key_Z);//, 1);
//	pEdit->addAction( "&Redo", pParent, SLOT(redoEdit()), CTRL+Key_Y, 2);
    pEdit->addSeparator();
    addAction("&Edit", this ,SLOT(pEdit));//, 2);
	Disable(2, 1);
//	Disable(2, 2);

    pView->addAction( "&Info", pParent, SLOT(viewInfo()), Qt::ALT+Qt::Key_0);//, 1);
    pView->addAction( "&Tone Mapping", pParent, SLOT(viewRight()), Qt::ALT+Qt::Key_1);//, 2);
    pView->addAction( "&Histogram", pParent, SLOT(viewHistogram()), Qt::ALT+Qt::Key_2);//, 3);
    pView->addSeparator();
    pView->addAction( "Zoom &In", pParent, SLOT(zoomIn()), Qt::CTRL+Qt::Key_Plus);//, 4);
    pView->addAction( "Zoom &Out", pParent, SLOT(zoomOut()), Qt::CTRL+Qt::Key_Minus);//, 5);
    pView->addAction( "Fit To &Screen", pParent, SLOT(fitToScreen()), Qt::ALT+Qt::Key_3);//, 6);
    pView->addAction( "Fit &Width", pParent, SLOT(fitToWidth()), Qt::ALT+Qt::Key_4);//, 7);
    pView->addAction( "Fit &Height", pParent, SLOT(fitToHeight()), Qt::ALT+Qt::Key_5);//, 8);
    addAction("&View", this, SLOT(pView));//, 3);
    // TODO pView->setCheckable(true);

    pComponent->addAction("&Red");//, 0);
    pComponent->addAction("&Green");//, 1);
    pComponent->addAction("&Blue");//, 2);

    pCommand->addAction( *TMOResource::pResource->IconDuplicate->pixmap(), "&Duplicate Image", pParent, SLOT(duplicateCommand()), Qt::CTRL+Qt::Key_D);//, 1);
    pCommand->addAction( *TMOResource::pResource->IconSize->pixmap(), "Change Image &Size", pParent, SLOT(sizeCommand()), 0);//, 2);
    pCommand->addSeparator();
    pCommand->addAction( "&Extract Luminance", pParent, SLOT(extractLumCommand()), Qt::CTRL+Qt::Key_L);//, 3);
    pCommand->addAction( "Extract &Component", pComponent, SLOT(pComponent));//, 4);
    pCommand->addAction( "&Merge components", pParent, SLOT(mergeCommand()), 0);//, 5);
    pCommand->addSeparator();
    pCommand->addAction( "Arithmetical &Operation", pParent, SLOT(operationCommand()), 0);//, 6);
    addAction("&Command", pCommand, SLOT(pCommand));//, 5);

    pHelpIt->addAction( "&Help", pParent, SLOT(showHelp()), Qt::Key_F1);//, 0);//, 1);
    pHelpIt->addSeparator();
    pHelpIt->addAction( "&About", this, SLOT(about()), 0);//, 3);
    pHelpIt->addSeparator();
    pHelpIt->addAction( *TMOResource::pResource->IconThis->pixmap(), "&What's This?", pParent,  SLOT(whatsThis()), Qt::SHIFT+Qt::Key_F1);//, 5);
    addAction("&Help", pHelpIt, SLOT(pHelpIt));//, 5);

	pImage = 0;
	SetWindows(0);
	connect (pWindows, SIGNAL(activated(int)), pParent, SLOT(activateWindow(int)));
	connect (pRecent, SIGNAL(activated(int)), pParent, SLOT(openFile(int)));
	connect (pComponent, SIGNAL(activated(int)), pParent, SLOT(extractComCommand(int)));
	
	return 0;
}

int TMOGUIMenu::SetWindows(QWorkspace* w)
{
	QWidgetList wl; 
	QWidget* widget;
	QString s;
	int number = 64;

	if (w) wl = w->windowList();
    // TODO removeItem(4);
	pWindows->clear();

	if(w)
	{
     pWindows->addAction( *TMOResource::pResource->IconCascade->pixmap(), "&Cascade", w, SLOT(cascade()), 0);//, 1);
     pWindows->addAction( *TMOResource::pResource->IconTile->pixmap(), "&Tile", w, SLOT(tile()), 0);//, 2);
	}
    pWindows->addAction( "Close &All", pParent, SLOT(closeallWindow()), 0);//, 3);
    pWindows->addSeparator();

    // TODO check
    foreach(widget, wl){
        s = widget->accessibleName();
        s = TMOGUIImage::GetName(s);
        QAction* act = pWindows->addAction(s);

        QString h = pImage->windowTitle();
        if(s == h)
            act->setChecked(true);
    }

	if (number == 64)
	{		
		Disable(1,3);
		Disable(1,4);
		Disable(1,5);
		Disable(1,6);
		Disable(1,7);
		Disable(1,8);
		if (listRecent.count()==0) Disable(1,9);
		Disable(2,1);
		Disable(2,2);
		Disable(2,3);
		Disable(3,3);
		Disable(3,4);
		Disable(3,5);
		Disable(3,6);
		Disable(3,7);
		Disable(3,8);
		Disable(4,1);
		Disable(4,2);
		Disable(4,3);
		Disable(4,4);
		Disable(4,5);
		Disable(4,6);
		Disable(5,1);
		Disable(5,2);
		Disable(5,3);

	}
	else
	{		
		Enable(1,3);
		Enable(1,4);
		Enable(1,5);
		Enable(1,6);
		Enable(1,7);
		Enable(1,8);
		if (listRecent.count()) Enable(1,9);
		Enable(3,3);
		Enable(3,4);
		Enable(3,5);
		Enable(3,6);
		Enable(3,7);
		Enable(3,8);
		Enable(4,1);
		Enable(4,2);
		Enable(4,3);
		Enable(4,4);
		Enable(4,5);
		Enable(4,6);
		Enable(5,1);
		Enable(5,2);
		Enable(5,3);
	}
    addAction("&Window", pWindows, SLOT(pWindows));//addAction("&Window", pWindows, 4, 4);
	return 0;
}

void TMOGUIMenu::windowChanged(TMOGUIImage* pImg)
{
	if(!pImg)
		return;
	pImage = pImg;
	int index = 4, id;
    /* TODO! while ( (id = pWindows->idAt(index++)) != -1 )
	{
		pWindows->setItemChecked(id, false);
		QString s = pWindows->text(id);
		if(s == pImage->caption())
			pWindows->setItemChecked(id, true);

	}
    pView->setItemChecked(pView->idAt(2), pImage->pToolsButton->isOn());*/
}

int TMOGUIMenu::Enable(int menu, int item)
{
    QMenu* pMenu;
	switch (menu)
	{
	case 1:
		pMenu = pFile;
		break;
	case 2:
		pMenu = pEdit;
		break;
	case 3:
		pMenu = pView;
		break;
	case 4:
		pMenu = pCommand;
		break;
	case 5:
		pMenu = pWindows;
		break;
	default:
		return 1;
	}
    // TODO! pMenu->setItemEnabled(item, true);
	return 0;
}

int TMOGUIMenu::Disable(int menu, int item)
{
    QMenu* pMenu;
	switch (menu)
	{
	case 1:
		pMenu = pFile;
		break;
	case 2:
		pMenu = pEdit;
		break;
	case 3:
		pMenu = pView;
		break;
	case 4:
		pMenu = pCommand;
		break;
	case 5:
		pMenu = pWindows;
		break;
	default:
		return 1;
	}
    //TODO! pMenu->setItemEnabled(item, false);
	return 0;
}

int TMOGUIMenu::LoadRecent()
{
	QFile f("recent.dat");

	if (f.open(QIODevice::ReadOnly))
	{
        Q3TextStream t( &f );        
        QString s;
		listRecent.clear();
        int n = 1;
        while ( !t.eof() )			
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
        Q3TextStream t( &f );        
 // TODO foreach
        for (QString* temp : listRecent)
		{
            t<<*temp<<'\n';
        }
        f.close();
		return 0;
	}
	return -1;
}

int TMOGUIMenu::AddRecent(QString& s)
{

// TODO foreach
    for (QString* temp : listRecent)
	{
		if (*temp == s) 
		{
            listRecent.removeAt(listRecent.indexOf(temp));
			break;
		}
	}
	if (listRecent.count() > 4) listRecent.removeLast();
	listRecent.prepend(new QString(s));

	return 0;
}

int TMOGUIMenu::SetRecent()
{

	int number = 128;

	pRecent->clear();
    for (QString* s : listRecent)
	{
        pRecent->addAction( *s);
        number++;
	}
	
	return 0;
}

QString TMOGUIMenu::GetRecent(int ID)
{

	int number = 128;

    for (QString* s : listRecent)
	{
		if (number == ID) return *s;
        number++;
	}
	return QString("");
}

void TMOGUIMenu::about()
{
    QMessageBox::about( this, "Tone Mapping Studio",
		"<p><b>(c) 2004-2015 <br />Ondrej Hajdok, Martin Cadik</b></p>" 
		"<p>"
		"<a href=\"mailto:cadikm@centrum.cz\">cadikm@centrum.cz</a></p>"
	);
}
