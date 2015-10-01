// TMOGUIMenu.cpp: implementation of the TMOGUIMenu class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIMenu.h"
#include "TMOGUIImage.h"
#include "TMOGUIResource.h"
#include <qpopupmenu.h>
#include <qworkspace.h>
#include <qfile.h>
#include <qtextstream.h>
#include <qlabel.h>
#include <qmessagebox.h>
#include <qpushbutton.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIMenu::TMOGUIMenu(QWidget * parent, const char * name):QMenuBar(parent, name)
{
	pParent = parent;
	pRecent = 0;
	listRecent.setAutoDelete(true);
	Create();
}

TMOGUIMenu::~TMOGUIMenu()
{
	SaveRecent();
}

int TMOGUIMenu::Create()
{
	pRecent = new QPopupMenu(this, "Recent");
	pComponent = new QPopupMenu(this, "Component");
	pFile = new QPopupMenu(this, "File");
	pEdit = new QPopupMenu(this, "Edit");
	pView = new QPopupMenu(this, "View");
	pCommand = new QPopupMenu(this, "Command");
	pWindows = new QPopupMenu(this, "Window");
	pHelpIt = new QPopupMenu(this, "Help");

	LoadRecent();

	pFile->insertItem( *TMOResource::pResource->IconNew->pixmap(), "&New",  pParent, SLOT(newFile()), CTRL+Key_N, 1 );
    pFile->insertItem( *TMOResource::pResource->IconOpen->pixmap(), "&Open", pParent, SLOT(openFile()), CTRL+Key_O, 2 );
	pFile->insertItem( "&Close", pParent, SLOT(closeFile()), 0, 3);
	pFile->insertSeparator();
	pFile->insertItem( *TMOResource::pResource->IconSave->pixmap(), "&Save",  pParent, SLOT(saveFile()), CTRL+Key_S, 4 );
    pFile->insertItem( *TMOResource::pResource->IconSave->pixmap(), "Save &As...", pParent, SLOT(saveasFile()), 0, 5 );
	pFile->insertItem( *TMOResource::pResource->IconSaveAll->pixmap(), "Save A&ll", pParent, SLOT(saveallFile()), 0, 6 );
	pFile->insertSeparator();
	pFile->insertItem( "Page Set&up...",  pParent, SLOT(pageFile()), 0, 7 );
    pFile->insertItem( *TMOResource::pResource->IconPrint->pixmap(), "&Print...", pParent, SLOT(printFile()), CTRL+Key_P, 8 );
	pFile->insertSeparator();
	pFile->insertItem( "Recent &Files", pRecent, 9 );
	pFile->insertSeparator();
	pFile->insertItem( "E&xit",  pParent, SLOT(exitFile()), 0, 10);	
	insertItem("&File", pFile, 1);

	pEdit->insertItem( *TMOResource::pResource->IconUndo->pixmap(), "&Undo", pParent, SLOT(undoEdit()), CTRL+Key_Z, 1);
//	pEdit->insertItem( "&Redo", pParent, SLOT(redoEdit()), CTRL+Key_Y, 2);
	pEdit->insertSeparator();
	insertItem("&Edit", pEdit, 2);
	Disable(2, 1);
//	Disable(2, 2);

	pView->insertItem( "&Info", pParent, SLOT(viewInfo()), ALT+Key_0, 1);
	pView->insertItem( "&Tone Mapping", pParent, SLOT(viewRight()), ALT+Key_1, 2);
	pView->insertItem( "&Histogram", pParent, SLOT(viewHistogram()), ALT+Key_2, 3);
	pView->insertSeparator();
	pView->insertItem( "Zoom &In", pParent, SLOT(zoomIn()), CTRL+Key_Plus, 4);
	pView->insertItem( "Zoom &Out", pParent, SLOT(zoomOut()), CTRL+Key_Minus, 5);
	pView->insertItem( "Fit To &Screen", pParent, SLOT(fitToScreen()), ALT+Key_3, 6);
	pView->insertItem( "Fit &Width", pParent, SLOT(fitToWidth()), ALT+Key_4, 7);
	pView->insertItem( "Fit &Height", pParent, SLOT(fitToHeight()), ALT+Key_5, 8);
	insertItem("&View", pView, 3);
	pView->setCheckable(true);

	pComponent->insertItem("&Red", 0);
	pComponent->insertItem("&Green", 1);
	pComponent->insertItem("&Blue", 2);

	pCommand->insertItem( *TMOResource::pResource->IconDuplicate->pixmap(), "&Duplicate Image", pParent, SLOT(duplicateCommand()), CTRL+Key_D, 1);
	pCommand->insertItem( *TMOResource::pResource->IconSize->pixmap(), "Change Image &Size", pParent, SLOT(sizeCommand()), 0, 2);
	pCommand->insertSeparator();
	pCommand->insertItem( "&Extract Luminance", pParent, SLOT(extractLumCommand()), CTRL+Key_L, 3);
	pCommand->insertItem( "Extract &Component", pComponent, 4);
	pCommand->insertItem( "&Merge components", pParent, SLOT(mergeCommand()), 0, 5);
	pCommand->insertSeparator();
	pCommand->insertItem( "Arithmetical &Operation", pParent, SLOT(operationCommand()), 0, 6);
	insertItem("&Command", pCommand, 5);

	pHelpIt->insertItem( "&Help", pParent, SLOT(showHelp()), Key_F1, 0, 1);
	pHelpIt->insertSeparator();
	pHelpIt->insertItem( "&About", this, SLOT(about()), 0, 3);
	pHelpIt->insertSeparator();
	pHelpIt->insertItem( *TMOResource::pResource->IconThis->pixmap(), "&What's This?", pParent,  SLOT(whatsThis()), SHIFT+Key_F1, 5);
	insertItem("&Help", pHelpIt, 5);	

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
	removeItem(4);
	pWindows->clear();

	if(w)
	{
	 pWindows->insertItem( *TMOResource::pResource->IconCascade->pixmap(), "&Cascade", w, SLOT(cascade()), 0, 1);
	 pWindows->insertItem( *TMOResource::pResource->IconTile->pixmap(), "&Tile", w, SLOT(tile()), 0, 2);
	}
	pWindows->insertItem( "Close &All", pParent, SLOT(closeallWindow()), 0, 3);
	pWindows->insertSeparator();

	for (widget = wl.first(); widget; widget = wl.next())
	{
		s = widget->name();
		s = TMOGUIImage::GetName(s);
		int id = pWindows->insertItem( s, number++);
		QString h = pImage->caption();
		if(s == h)
			pWindows->setItemChecked(id, true);
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
	insertItem("&Window", pWindows, 4, 4);
	return 0;
}

void TMOGUIMenu::windowChanged(TMOGUIImage* pImg)
{
	if(!pImg)
		return;
	pImage = pImg;
	int index = 4, id;
	while ( (id = pWindows->idAt(index++)) != -1 ) 
	{
		pWindows->setItemChecked(id, false);
		QString s = pWindows->text(id);
		if(s == pImage->caption())
			pWindows->setItemChecked(id, true);

	}
	pView->setItemChecked(pView->idAt(2), pImage->pToolsButton->isOn());
}

int TMOGUIMenu::Enable(int menu, int item)
{
	QPopupMenu* pMenu;
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
	pMenu->setItemEnabled(item, true);
	return 0;
}

int TMOGUIMenu::Disable(int menu, int item)
{
	QPopupMenu* pMenu;
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
	pMenu->setItemEnabled(item, false);
	return 0;
}

int TMOGUIMenu::LoadRecent()
{
	QFile f("recent.dat");

	if (f.open(IO_ReadOnly))
	{
        QTextStream t( &f );        
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
	QString *temp;
	if (f.open(IO_WriteOnly))
	{
        QTextStream t( &f );        
 
        for ( temp = listRecent.first(); temp; temp = listRecent.next() )			
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
	QString *temp;

	for (temp = listRecent.first(); temp; temp = listRecent.next())
	{
		if (*temp == s) 
		{
			listRecent.remove(temp);
			break;
		}
	}
	if (listRecent.count() > 4) listRecent.removeLast();
	listRecent.prepend(new QString(s));

	return 0;
}

int TMOGUIMenu::SetRecent()
{
	QString* s;
	int number = 128;

	pRecent->clear();
	for (s = listRecent.first(); s ; s = listRecent.next())
	{
		pRecent->insertItem( *s, number++);
	}
	
	return 0;
}

QString TMOGUIMenu::GetRecent(int ID)
{
	QString* s;
	int number = 128;

	for (s = listRecent.first(); s ; s = listRecent.next(), number++)
	{
		if (number == ID) return *s;
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
