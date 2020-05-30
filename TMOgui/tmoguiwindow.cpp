#include <qapplication.h>
#include <qmessagebox.h>
//#include <qassistantclient.h>
#include <qdir.h>
#include <qpushbutton.h>
//#include <q3hbox.h>
#include <QDockWidget>
#include <qsplitter.h>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QFileDialog>
#include <qcombobox.h>
#include <qcheckbox.h>
#include <qimage.h>
#include <qradiobutton.h>
#include <QLabel>
#include <qlineedit.h>
#include <QtPrintSupport/QPrinter>
#include <QtPrintSupport/QPrintDialog>
#include <qpainter.h>
#include <qregexp.h>
#include <qmatrix.h>
#include <QToolBar>
#include <qaction.h>
#include <QTextStream>
//Added by qt3to4:
#include <QList>
#include <QPixmap>
#include <QFrame>
#include "../tmolib/TMO.h"
#include "tmoguiwindow.h"
#include "TMOGUIToneMapping.h"
#include "TMOGUIStatistics.h"
#include "TMOGUIOutput.h"
#include "TMOGUIAdjust.h"
#include "TMOGUIToneSlider.h"
#include "TMOGUIResource.h"
#include "TMOGUIHisto.h"
#include "TMOGUIInfoTool.h"
#include "TMOGUITransformation.h"



#include "TMOGUISaveDialog.h"
#include <qmap.h>

//#include <iostream>

//Ui::TMOGUIResource* TMOResource::pResource = 0;

TMOGUIWindow::TMOGUIWindow(QWidget* parent, const char* name, Qt::WindowFlags f )
    : QMainWindow( parent, f )
{
    setup();
    this->setObjectName(name);
}
void TMOGUIWindow::setup()
{
	pMargins[0] = pMargins[1] = pMargins[2] = pMargins[3] = 10;

    // TMOResource::pResource = new TMOGUIResource(this, "Resources");
    // TMOResource::pResource->setupUi(this); // = new TMOGUIResource(this, "Resources");
    // TODO TMOResource::pResource->hide();

    setWindowTitle("TMOGUI");
    setWindowIcon(QIcon(QString::fromUtf8(":/resources/icons/IconMain.png"))); //*TMOResource::pResource->IconMain->pixmap());
	sPrevFileName = "";
	Create();
}

TMOGUIWindow::~TMOGUIWindow()
{
	SavePosition();
    delete iTool;
}

int TMOGUIWindow::Create()
{
    pMenu = new TMOGUIMenu(this, "Menu");
    pStatus = new TMOGUIStatus(this, "Status");
    pProgress = new TMOGUIProgressBar(pStatus, "Progress");
    pProgress->SetLabel("");
    pRightSplitter = new QSplitter(this);//, "RightSplitter");
    pRightSplitter->setHandleWidth(3);
    pRightSplitter->setStretchFactor(1,3);
    pRightSplitter->setChildrenCollapsible(false);
    pRightSplitter->setFrameStyle( QFrame::Sunken | QFrame::Panel );
    pSplitter = new QSplitter(pRightSplitter);//, "BottomSplitter");
    pSplitter->setOrientation(Qt::Vertical);
    pSplitter->setChildrenCollapsible(false);
    pWorkspace = new QMdiArea(pSplitter);//, "Workspace");
        pWorkspace->setViewMode(QMdiArea::ViewMode::SubWindowView);
        pWorkspace->setTabsMovable(true);
    pInfo = TMOGUIOutput::pInfo = new TMOGUIInfo(pSplitter, "Info");
    pRight = new TMOGUIRightBar(pRightSplitter, "RightBar");

    pInfo->bVisible = true;
    setCentralWidget(pRightSplitter);
	
    qDeleteAll(listImage);
    LoadPosition();

    //QToolBar *pToolBar = addToolBar("main toolbar");
    pFileTool = new TMOGUIFileToolBar(this);
    pTools = new TMOGUIZoomTool(this);
    pInfoTool = new TMOGUIInfoToolBar(this);
    pDialog = nullptr;
    iTool = new TMOGUIInfoTool(this);

    // assistant = new QAssistantClient( "" );
    //TODO
    // assistant = new QAssistantClient( QString(), 0, 0 );
    // setAssistantArguments();

    this->setMenuBar(pMenu);
    this->setStatusBar(pStatus);
    this->setCentralWidget(pRightSplitter);

    this->addToolBar(Qt::TopToolBarArea, pFileTool);
    this->addToolBarBreak();
    this->insertToolBar(pFileTool, pTools);
    this->addToolBarBreak();
    this->insertToolBar(pTools, pInfoTool);


    //connect( assistant, SIGNAL(error(const QString&)), this, SLOT(showAssistantErrors(const QString&)) );
    connect(pRight->GetMapping()->pOk, SIGNAL(clicked()), this, SLOT(transform()));
    connect(pMenu, &TMOGUIMenu::openFile, this, QOverload<QString>::of(&TMOGUIWindow::openFile));
    connect(pMenu, &TMOGUIMenu::activateWindowAction, this, QOverload<int>::of(&TMOGUIWindow::activateWindow));
    connect(pWorkspace, &QMdiArea::subWindowActivated, this, &TMOGUIWindow::windowChanged);
    connect(this, &TMOGUIWindow::imageSelected, pInfo->pStats, &TMOGUIStatistics::windowChanged);
    connect(this, &TMOGUIWindow::imageSelected, pRight->pFilters, &TMOGUIFilters::windowChanged);
    connect(this, &TMOGUIWindow::imageSelected, pTools, &TMOGUIZoomTool::windowChanged);
    connect(this, &TMOGUIWindow::imageSelected, pMenu, &TMOGUIMenu::windowChanged);
    connect(pRight, SIGNAL(closeBar()), this, SLOT(viewRight()));
    connect(pInfo, SIGNAL(closeBar()), this, SLOT(viewInfo()));
    connect(pWorkspace, &QMdiArea::close, this, &TMOGUIWindow::closeActiveWindow);
	return 0;
}

void TMOGUIWindow::showHelp()
{
    // TODO assistant->showPage( "../doc/index.html" );
}

void TMOGUIWindow::setAssistantArguments()
{
    QStringList cmdLst;
	cmdLst << "-profile" << "../doc/help.adp";
    // TODO assistant->setArguments( cmdLst );
}

void TMOGUIWindow::showAssistantErrors( const QString &err )
{
    QMessageBox::critical( this, "Assistant Error", err );
}

void TMOGUIWindow::exitFile()
{
	closeallWindow();
	qApp->quit();
}

void TMOGUIWindow::openFile(QString fileName)
{
	QString number, name;
    TMOGUIImage* newfile = nullptr;

	int iCount = 0;
	
        if ( fileName.isEmpty() ) return;

    for (TMOGUIImage* temp : listImage)
	{
        if (!temp->pImage || !temp->imageName || temp->imageName->isEmpty()){
            listImage.removeOne(temp);
            continue;
        }
        name = *temp->imageName;
		if (name.contains(fileName)) iCount++;
	}
	
	if (iCount)
	{
		number.setNum(iCount);
        QString newfilename = fileName + " [" + number + "]";
        newfile = new TMOGUIImage(pProgress, pWorkspace, newfilename.toStdString().c_str());
	}
	else
        newfile = new TMOGUIImage(pProgress, pWorkspace, fileName.toStdString().c_str());

    QMdiSubWindow* subw = pWorkspace->addSubWindow(newfile);
    //connect(subw, &QMdiSubWindow::windowStateChanged, this, &TMOGUIWindow::subWindowStateChanged);


    if (newfile->Open(fileName.toStdString().c_str())) delete newfile;
	else
	{
        listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));		
		WindowChangedToolActivated(newfile);
		pMenu->AddRecent(fileName);
		pMenu->SetRecent();
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

void TMOGUIWindow::openFile()
{

    openFile( QFileDialog::getOpenFileName
        ( this,  QString(), QString(),
          "HDR images (*.raw *.hdrraw *.tif *.tiff *.pic *.hdr *.exr *.jpg *.jpeg);; *.raw;; *.hdrraw;; *.tif;; *.tiff;; *.pic;; *.hdr;; *.exr;; *.jpg;; *.jpeg;; *.png;; *.ppm;; All files (*.*)",
          new QString("HDR images (*.raw *.hdrraw *.tif *.tiff *.pic *.hdr *.exr *.jpg *.jpeg)") ) // TODO Check
        	);
}

void TMOGUIWindow::saveasFile()
{
	QString fileName;
	TMOGUIImage* pImage;
	TMOImage OutImage;
	QString s;
		
    fileName = QString(pWorkspace->activeSubWindow()->widget()->objectName());
	pImage = FindImage(fileName);
	if (!pImage) return;

	int iIdx = 0, iFound;
    if ((iFound = fileName.indexOf("[")) > 0)
	{
		QString sTemp = fileName.mid(iFound + 1);
        sTemp = sTemp.left(sTemp.indexOf("]"));
		iIdx = sTemp.toInt();
		fileName = fileName.left(iFound);
	}
	
    if ((iFound = fileName.lastIndexOf(".", -1)) > 0)
	{
		QString sTemp;
		if (iIdx) fileName = fileName.left(iFound) + "[" + sTemp.setNum(iIdx) + "].tif";
		else fileName = fileName.left(iFound) + ".tif";
	}
	
	filterMap filter;
	filter["Tiff image (*.tif *.tiff)"] = fileType(TMO_TIFF_8,"tif");
	filter["Tiff image HDR(*.tif *.tiff)"] = fileType(TMO_TIFF_32,"tif");
	filter["raw HDR(*.raw)"] = fileType(TMO_RAW_32,"raw");
	filter["*.exr"] = fileType(TMO_EXR_16,"exr");
	filter["pfm HDR (*.pfm)"] = fileType(TMO_PFM_32,"pfm");
	filter["HDR (*.hdr)"] = fileType(TMO_HDR_32,"hdr");
	filter["jpeg (*.jpeg *.jpg)"] = fileType(TMO_JPEG_8,"jpeg");
	filter["png (*.png)"] = fileType(TMO_PNG_8,"png");
	filter["ppm (*.ppm)"] = fileType(TMO_PPM_8,"ppm");
	TMOGUISaveDialog * fd = new TMOGUISaveDialog(fileName,&filter,this,"",true);
    fd->setOption(QFileDialog::DontUseNativeDialog, true);

    //fileName = fd->getSaveFileName(this, fileName, fileName ,*fd->filtersString, &defaultFilter);
    if(fd->exec() == QDialog::Accepted)
    {
        fileName = *fd->file;
         OutImage.New(*pImage->GetImage());
         pImage->pImage->GetImage(&OutImage);
         if(OutImage.SaveAs(fileName.toStdString().c_str(), fd->getSelectedFileType()))
          QMessageBox::critical( nullptr, "Tone Mapping Studio", QString("Failed to save file : \n\n") + fileName + "\n(probably unknown format)");
             else
         {
          pImage->GetImage()->WriteLine(L"");
          s = QString("File saved to : ") + fileName;
          pImage->GetImage()->WriteLine(pImage->GetString(s.unicode()));
         }

	}

}

void TMOGUIWindow::saveFile()
{
	QString fileName;
	TMOGUIImage *pImage;
	TMOImage OutImage;
	int iFound;

    fileName = pWorkspace->activeSubWindow()->widget()->objectName();
	pImage = FindImage(fileName);
	if (!pImage) return;

	int iIdx = 0;
    if ((iFound = fileName.indexOf("[")) > 0)
	{
		QString sTemp = fileName.mid(iFound + 1);
        sTemp = sTemp.left(sTemp.indexOf("]"));
		iIdx = sTemp.toInt();
		fileName = fileName.left(iFound);
	}
	
    if ((iFound = fileName.lastIndexOf(".", -1)) > 0)
	{
		QString sTemp;
		if (iIdx) fileName = fileName.left(iFound) + "[" + sTemp.setNum(iIdx) + "].tif";
		else fileName = fileName.left(iFound) + ".tif";
	}

	OutImage.New(*pImage->GetImage());
	pImage->pImage->GetImage(&OutImage);

    OutImage.SaveAs(fileName.toStdString().c_str());
    if (fileName == pWorkspace->activeSubWindow()->widget()->objectName())
	{
		pImage->GetImage()->WriteLine(L"");
		pImage->GetImage()->WriteLine(L"Image saved to the original location.");
	}
	else
	{
		pImage->GetImage()->WriteLine(L"");
		QString s = QString("File saved to : ") + fileName;
		pImage->GetImage()->WriteLine(pImage->GetString(s.unicode()));
	}
}

void TMOGUIWindow::openFile(int ID)
{
	int iCount = 0;
    TMOGUIImage* newfile = nullptr;
	QString number,fileName = pMenu->GetRecent(ID), name;


    if ( fileName.isEmpty() ) return;

    for (TMOGUIImage* temp : listImage)
    {
        if (!temp->pImage || !temp->imageName || temp->imageName->isEmpty()){
            listImage.removeOne(temp);
            continue;
        }
        name = *temp->imageName;
		if (name.contains(fileName)) iCount++;
	}
	
	if (iCount)
	{
		number.setNum(iCount);
        newfile = new TMOGUIImage(pProgress, pWorkspace, (fileName + " [" + number + "]").toStdString().c_str());
    }
    else
        newfile = new TMOGUIImage(pProgress, pWorkspace, fileName.toStdString().c_str());

    pWorkspace->addSubWindow(newfile);

    if (newfile->Open(fileName.toStdString().c_str())) delete newfile;
	else
	{
        listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));
		WindowChangedToolActivated(newfile);
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

void TMOGUIWindow::closeActiveWindow()
{
    QString s;
    s = pWorkspace->activeSubWindow()->widget()->objectName();
    listImage.removeOne(FindImage(s));
    pWorkspace->close();
}

void TMOGUIWindow::closeFile()
{
	QString s;
	TMOGUIImage* im;

    s = pWorkspace->activeSubWindow()->widget()->objectName();
	im = FindImage(s);
	if (im->Terminate()) return;
    listImage.removeOne(im);
	
	pTools->SetWindows(pWorkspace);
	pInfoTool->SetWindows(pWorkspace);
	pFileTool->SetWindows(pWorkspace);
	pMenu->SetWindows(pWorkspace);
}

void TMOGUIWindow::closeallWindow()
{
	listImage.clear();

    pWorkspace->closeAllSubWindows();
	pMenu->SetWindows(pWorkspace);
	pTools->SetWindows(pWorkspace);
    pInfo->SetOutput(nullptr);
}

void TMOGUIWindow::saveallFile()
{
	QString fileName;

	TMOImage OutImage;
	int iFound;

    refreshWindowsList();
    for (TMOGUIImage *pImage : listImage)
	{
        fileName = *pImage->imageName;
		if (!pImage) return;

        if ((iFound = fileName.indexOf("[")) > 0)
		{
			fileName = fileName.left(iFound);
		}
		
        if ((iFound = fileName.lastIndexOf(".", -1)) > 0)
		{
			fileName = fileName.left(iFound) + ".tif";
		}

		OutImage.New(*pImage->GetImage());
		pImage->pImage->GetImage(&OutImage);

        OutImage.SaveAs(fileName.toStdString().c_str());
        if (fileName == pWorkspace->activeSubWindow()->widget()->objectName())
		{
			pImage->GetImage()->WriteLine(L"");
			pImage->GetImage()->WriteLine(L"Image saved to the original location.");
		}
		else
		{
			pImage->GetImage()->WriteLine(L"");
			QString s = QString("File saved to : ") + fileName;
			pImage->GetImage()->WriteLine(pImage->GetString(s.unicode()));
		}
	}
}

TMOGUIImage* TMOGUIWindow::FindImage(QString name)
{


    for (TMOGUIImage* retval : listImage)
	{
        if (retval->imageName == name) return retval;
	}
    return nullptr;
}

void TMOGUIWindow::activateWindow(int id)
{
    QList<QMdiSubWindow *> wl;

    QString sName;
	TMOGUIImage *pImage;
	
    //?int number = 64;
    int number = 0;

    refreshWindowsList();
    wl = pWorkspace->subWindowList();

    for (QMdiSubWindow* widget : wl)
	{
		if (number++ == id)
		{
			widget->setFocus();			
            sName = widget->widget()->objectName(); //name()
			pImage = FindImage(sName);
            if (pImage->CanUndo()) pMenu->Enable(2, 1);
            else pMenu->Disable(2, 1);
			pInfo->SetOutput(pImage->pOutput);
			pMenu->windowChanged(pImage);
			break;
		}
	}
}

void TMOGUIWindow::transform()
{
    TMOImage *pSrc = nullptr, *pDst = nullptr;
    TMOGUIImage *pImage = nullptr;
	QWidget *pWidget;
	TMOGUITransformation *pTransform;
	TMO* pTMO;
	QString sName;

	pTMO = pRight->GetTMO();
	if (!pTMO) return;
    if (!pWorkspace->activeSubWindow()) return;
    pWidget = pWorkspace->activeSubWindow()->widget();
    if (!pWidget) return;
    sName = pWidget->objectName(); //name()
	pImage = FindImage(sName);
	if (!pImage) return;

	pTransform = pImage->Transform();
	if (pTransform)
	{
		pSrc = pImage->GetImage();
        pTMO->SetImage(*pSrc);
        //pTMO->SetSource(*pSrc);
		pTransform->SetTMO(pTMO);
	}
}

void TMOGUIWindow::finishTransform()
{
    TMOGUIImage *pImage = nullptr;
	QString sName;

    sName = pWorkspace->activeSubWindow()->widget()->objectName();//name()
	pImage = FindImage(sName);
	if (pImage->CanUndo()) pMenu->Enable(2, 1);
	else pMenu->Disable(2, 1);
}

void TMOGUIWindow::undoEdit()
{
    TMOGUIImage *pImage = nullptr;
	QString sName;

    if(!pWorkspace || !(pWorkspace->activeSubWindow()->widget()))return;
    sName = pWorkspace->activeSubWindow()->widget()->objectName();//name()
	pImage = FindImage(sName);
	if (!pImage) return;

	pImage->Swap();
}

void TMOGUIWindow::windowChanged(QMdiSubWindow* pWidget)
{
	if (!pWidget) 
	{
        refreshWindowsList();
        emit imageSelected(nullptr);
	}
	else
	{

        QString sName = pWidget->widget()->objectName();//name()
		TMOGUIImage* pImage = FindImage(sName);
		if (!pImage) return;
		pInfo->SetOutput(pImage->pOutput);
		emit imageSelected(pImage);
		if (pImage->IsModified()) pMenu->Enable(2, 1);
		else pMenu->Disable(2,1);
		WindowChangedToolActivated(pImage);
	}
}

void TMOGUIWindow::duplicateCommand()
{
    QWidget *pWindow = pWorkspace->activeSubWindow()->widget();
	if (!pWindow) return;
    QString sName = pWindow->objectName();//name()
	TMOGUIImage* pImage = FindImage(sName);
	if (!pImage) return;

	TMOGUIImage *newfile = GetNewImage(sName);

	if (newfile)
	{
		newfile->New(pImage);
		listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));
		WindowChangedToolActivated(newfile);
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

TMOGUIImage* TMOGUIWindow::GetNewImage(const QString &sName)
{
	int iCount = 0;
    TMOGUIImage* newfile = nullptr;
	QString number,fileName = sName, name;
	
    if ( fileName.isEmpty() ) return nullptr;
	
    QRegExp r = QRegExp(" \\[[0123456789]+\\]", Qt::CaseSensitive, QRegExp::RegExp);
    int iFind = fileName.indexOf(r);
	fileName = fileName.left(iFind);
    
    for (TMOGUIImage* temp : listImage)
	{
        if (!temp->pImage || !temp->imageName || temp->imageName->isEmpty()){
            listImage.removeOne(temp);
            continue;
        }
        name = *temp->imageName;

        if (name.contains(fileName)) iCount++; // CHECK != NULL
	}

	if (iCount)
	{
		number.setNum(iCount);
        newfile = new TMOGUIImage(pProgress, pWorkspace, (fileName + " [" + number + "]").toStdString().c_str());

    }
	else
        newfile = new TMOGUIImage(pProgress, pWorkspace, fileName.toStdString().c_str());

    pWorkspace->addSubWindow(newfile);

	return newfile;
}

void TMOGUIWindow::sizeCommand()
{
	TMOGUIImage* pImage = GetActiveImage();
    QDialog* qDialog = new QDialog(this);
	QString s;
	int iWidth, iHeight;
	if (!pImage) return;
	
    pDialog = new Ui::TMOGUIImageSize();
    pDialog->setupUi(qDialog);
	iFlags = iFlags & ~1;
	pDialog->ratioLabel->setHidden(true);
	dRatio = (double)pImage->GetImage()->GetWidth() / pImage->GetImage()->GetHeight();
	pDialog->LineEdit1->setText(s.setNum(pImage->GetImage()->GetWidth()));
	pDialog->LineEdit2->setText(s.setNum(pImage->GetImage()->GetHeight()));
    connect (pDialog->LineEdit1, &QLineEdit::textChanged, this, &TMOGUIWindow::ImageSizeWidth);
    connect (pDialog->LineEdit2, &QLineEdit::textChanged, this, &TMOGUIWindow::ImageSizeHeight);
    connect (pDialog->CheckBox1, &QCheckBox::toggled, this, &TMOGUIWindow::ImageSizeConstrain);
    connect (this, &TMOGUIWindow::signalImageSizeWidth, pDialog->LineEdit1, &QLineEdit::setText);
    connect (this, &TMOGUIWindow::signalImageSizeHeight, pDialog->LineEdit2, &QLineEdit::setText);
    if (qDialog->exec() == QDialog::Rejected) return;
	iWidth = pDialog->LineEdit1->text().toInt();
	iHeight = pDialog->LineEdit2->text().toInt();
	pImage->SetImageSize(iWidth, iHeight);	
	delete pDialog;
    pDialog = nullptr;
}

TMOGUIImage* TMOGUIWindow::GetActiveImage()
{
    QWidget *pWindow = pWorkspace->activeSubWindow()->widget();
    if (!pWindow) return nullptr;
    QString sName = pWindow->objectName();
	return FindImage(sName);
}

void TMOGUIWindow::ImageSizeConstrain(bool bChecked)
{
	QString s;
	TMOGUIImage* pImage = GetActiveImage();
	if (!pImage && !pDialog) return;
	if (bChecked) 
	{
		iFlags = iFlags | 7;
		pDialog->LineEdit1->setText(s.setNum(pImage->GetImage()->GetWidth()));
		pDialog->LineEdit2->setText(s.setNum(pImage->GetImage()->GetHeight()));
	}
	else 
		iFlags = iFlags & ~1;

}

void TMOGUIWindow::ImageSizeWidth(const QString &s)
{
	if ((iFlags&1) && (iFlags&4))
	{
		int i = s.toInt();
		iFlags = iFlags & ~2;
		QString outstring;

		i /= dRatio;
		outstring.setNum(i);
		emit signalImageSizeHeight(outstring);
		iFlags = iFlags | 2;
	}
}

void TMOGUIWindow::ImageSizeHeight(const QString &s)
{
	if ((iFlags&1) && (iFlags&2))
	{
		int i = s.toInt();
		iFlags = iFlags & ~4;
		QString outstring;

		i *= dRatio;
		outstring.setNum(i);
		emit signalImageSizeWidth(outstring);
		iFlags = iFlags | 4;
	}
}

void TMOGUIWindow::extractLumCommand()
{
	TMOGUIImage* pImage = GetActiveImage();
	if (!pImage) return;
    TMOGUIImage *newfile = GetNewImage(*pImage->imageName);

	if (newfile)
	{
		newfile->Extract(pImage, 3);
		listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));
		WindowChangedToolActivated(newfile);
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

void TMOGUIWindow::extractRed(){
    extractComCommand(0);
}
void TMOGUIWindow::extractGreen(){
    extractComCommand(1);
}
void TMOGUIWindow::extractBlue(){
    extractComCommand(2);
}

void TMOGUIWindow::extractComCommand(int iComponent)
{
	TMOGUIImage* pImage = GetActiveImage();
	if (!pImage) return;
    TMOGUIImage *newfile = GetNewImage(*pImage->imageName);

	if (newfile)
	{
		newfile->Extract(pImage, iComponent);
		listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));
		WindowChangedToolActivated(newfile);
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

void TMOGUIWindow::mergeCommand()
{
    Ui::TMOGUIMergeComponents *pDialog;
    QDialog *qDialog = new QDialog(this);

    pDialog = new Ui::TMOGUIMergeComponents();
    pDialog->setupUi(qDialog);
    connect (pDialog->ComboBox1, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIWindow::MergeComponentsRed);
    connect (pDialog->ComboBox2, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIWindow::MergeComponentsGreen);
    connect (pDialog->ComboBox3, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIWindow::MergeComponentsBlue);
    connect (this, &TMOGUIWindow::signalMergeComRed, pDialog->PixmapLabel1, &QLabel::setPixmap);
    connect (this, &TMOGUIWindow::signalMergeComGreen, pDialog->PixmapLabel2, &QLabel::setPixmap);
    connect (this, &TMOGUIWindow::signalMergeComBlue, pDialog->PixmapLabel3, &QLabel::setPixmap);
		
	QString s;
	int i = 0;
    for (TMOGUIImage *temp: listImage)
	{
        if (!temp->pImage || temp->imageName == nullptr){
            listImage.removeOne(temp);
            continue;
        }
        s = TMOGUIImage::GetName(*temp->imageName);
        pDialog->ComboBox1->insertItem(i, s);
        pDialog->ComboBox2->insertItem(i, s);
        pDialog->ComboBox3->insertItem(i, s);
        i++;
	}

	MergeComponentsRed(0);
	MergeComponentsGreen(0);
	MergeComponentsBlue(0);

	s = "";
    if (pImages[0]) s.append(TMOGUIImage::GetName(*pImages[0]->imageName) + "_");
    if (pImages[1]) s.append(TMOGUIImage::GetName(*pImages[1]->imageName) + "_");
    if (pImages[2]) s.append(TMOGUIImage::GetName(*pImages[2]->imageName) + "_");
	s = s.left(s.length());

    // TODO if (pDialog->exec() == QDialog::Rejected) return;

    TMOGUIImage *newfile = new TMOGUIImage(pProgress, pWorkspace, s.toStdString().c_str());

    pWorkspace->addSubWindow(newfile);

	if (newfile)
	{
		newfile->MergeComponents(pImages[0], pImages[1], pImages[2]);	
		listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));
		WindowChangedToolActivated(newfile);
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

void TMOGUIWindow::OperationFirst(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	
    QList<TMOGUIImage*>::const_iterator listImageIter = listImage.begin();
    temp = (*listImageIter); // TODO listImage.next()
    for (int i = iImage; temp && i; i--)
    {
        listImageIter++;
        temp = (*listImageIter);
    }


	pImages[0] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
    QImage tempImage(160, 160, QImage::Format::Format_RGB32);
    QMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->transformed(QTransform(m));
	emit signalMergeComRed(aPixmap);
}

void TMOGUIWindow::OperationSecond(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	
    QList<TMOGUIImage*>::const_iterator listImageIter = listImage.begin();
    temp = (*listImageIter); // TODO listImage.next()
    for (int i = iImage; temp && i; i--)
    {
        listImageIter++;
        temp = (*listImageIter);
    }

	pImages[1] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
    QImage tempImage(160, 160, QImage::Format::Format_RGB32);
    QMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->transformed(QTransform(m));
	emit signalMergeComGreen(aPixmap);
}

void TMOGUIWindow::MergeComponentsRed(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	QRgb* pRgb;
	
    QList<TMOGUIImage*>::const_iterator listImageIter = listImage.begin();
    temp = (*listImageIter); // TODO listImage.next()
    for (int i = iImage; temp && i; i--)
    {
        listImageIter++;
        temp = (*listImageIter);
    }

	pImages[0] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
    QImage tempImage(160, 160, QImage::Format::Format_RGB32);
    QMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->transformed(QTransform(m));
	
    tempImage = aPixmap.toImage();
	for (int jx = 0; jx < aPixmap.height(); jx++)
	{
		pRgb = (QRgb*)tempImage.scanLine(jx);
		for (int ix = 0; ix < aPixmap.width(); ix++, pRgb++)
		{
			*pRgb = qRgb( qRed(*pRgb), 0 , 0);
		}
	}
	aPixmap.convertFromImage(tempImage, Qt::ColorOnly);
	
	emit signalMergeComRed(aPixmap);
}

void TMOGUIWindow::MergeComponentsGreen(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	QRgb* pRgb;
	
    QList<TMOGUIImage*>::const_iterator listImageIter = listImage.begin();
    temp = (*listImageIter); // TODO listImage.next()
    for (int i = iImage; temp && i; i--)
    {
        listImageIter++;
        temp = (*listImageIter);
    }

	pImages[1] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
    QImage tempImage(160, 160, QImage::Format::Format_RGB32);
    QMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->transformed(QTransform(m));
	
    tempImage = aPixmap.toImage();
	for (int jx = 0; jx < aPixmap.height(); jx++)
	{
		pRgb = (QRgb*)tempImage.scanLine(jx);
		for (int ix = 0; ix < aPixmap.width(); ix++, pRgb++)
		{
			*pRgb = qRgb( 0, qGreen(*pRgb), 0);
		}
	}
	aPixmap.convertFromImage(tempImage, Qt::ColorOnly);
	emit signalMergeComGreen(aPixmap);
}

void TMOGUIWindow::MergeComponentsBlue(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	QRgb* pRgb;
	
    QList<TMOGUIImage*>::const_iterator listImageIter = listImage.begin();
    temp = (*listImageIter); // TODO listImage.next()
    for (int i = iImage; temp && i; i--)
    {
        listImageIter++;
        temp = (*listImageIter);
    }

	pImages[2] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
    QImage tempImage(160, 160, QImage::Format::Format_RGB32);
    QMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->transformed(QTransform(m));
	
    tempImage = aPixmap.toImage();
	for (int jx = 0; jx < aPixmap.height(); jx++)
	{
		pRgb = (QRgb*)tempImage.scanLine(jx);
		for (int ix = 0; ix < aPixmap.width(); ix++, pRgb++)
		{
			*pRgb = qRgb( 0, 0, qBlue(*pRgb));
		}
	}
	aPixmap.convertFromImage(tempImage, Qt::ColorOnly);
	emit signalMergeComBlue(aPixmap);
}

void TMOGUIWindow::operationCommand()
{
    Ui::TMOGUIOperation *pDialog;
    QDialog* qDialog = new QDialog();

    pDialog = new Ui::TMOGUIOperation();
    pDialog->setupUi(qDialog);
    connect (pDialog->ComboBox1, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIWindow::OperationFirst);
    connect (pDialog->ComboBox2, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIWindow::OperationSecond);
    connect (pDialog->ComboBox3, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIWindow::ImageOperation);
    connect (this, &TMOGUIWindow::signalMergeComRed, pDialog->PixmapLabel1, &QLabel::setPixmap);
    connect (this, &TMOGUIWindow::signalMergeComGreen, pDialog->PixmapLabel2, &QLabel::setPixmap);

		
	QString s;
	int i = 0;
    for (TMOGUIImage *temp : listImage)
	{
        if (!temp->pImage || !temp->imageName || temp->imageName->isEmpty()){
            listImage.removeOne(temp);
            continue;
        }
        s = TMOGUIImage::GetName(*temp->imageName);
        pDialog->ComboBox1->insertItem(i, s);
        pDialog->ComboBox2->insertItem(i, s);
        i++;
	}

    pDialog->ComboBox3->insertItem(0,"Addition");
    pDialog->ComboBox3->insertItem(1,"Subtraction");
    pDialog->ComboBox3->insertItem(2,"Multiplication");
    pDialog->ComboBox3->insertItem(3,"Division");
    pDialog->ComboBox3->insertItem(4,"Involution");
    pDialog->ComboBox3->insertItem(5,"Maximum");
    pDialog->ComboBox3->insertItem(6,"Minimum");
    pDialog->ComboBox3->insertItem(7,"Average");
	
	OperationFirst(0);
	OperationSecond(0);
	iFlags = 0;
	
	s = "";
    if (pImages[0]) s.append(TMOGUIImage::GetName(*pImages[0]->imageName) + "_");
    if (pImages[1]) s.append(TMOGUIImage::GetName(*pImages[1]->imageName));
	s = s.left(s.length());

    // TODO if (pDialog->exec() == QDialog::Rejected) return;

    TMOGUIImage *newfile = new TMOGUIImage(pProgress, pWorkspace, s.toStdString().c_str());

    pWorkspace->addSubWindow(newfile);

	if (newfile)
	{
		newfile->ImageOperation(pImages[0], pImages[1], iFlags);	
		listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));		
		WindowChangedToolActivated(newfile);
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

void TMOGUIWindow::ImageOperation(int iOperation)
{
	iFlags = iOperation;
}

void TMOGUIWindow::newFile()
{
    Ui::TMOGUINewFile *pDialog;
    QDialog *qDialog = new QDialog(this);
	QString s;
	
    pDialog = new Ui::TMOGUINewFile();
    pDialog->setupUi(qDialog); // = new TMOGUINewFile(this, "NewFileDialog", true);

    connect (pDialog->RadioButton1, &QRadioButton::toggled, this, &TMOGUIWindow::NewImageConstant);
    connect (pDialog->LineEdit1, &QLineEdit::textChanged, this, &TMOGUIWindow::SetFileName);
    connect (pDialog->LineEdit2, &QLineEdit::textChanged, this, &TMOGUIWindow::SetWidth);
    connect (pDialog->LineEdit3, &QLineEdit::textChanged, this, &TMOGUIWindow::SetHeight);
    connect (pDialog->LineEdit5, &QLineEdit::textChanged, this, &TMOGUIWindow::SetRMin);
    connect (pDialog->LineEdit6, &QLineEdit::textChanged, this, &TMOGUIWindow::SetGMin);
    connect (pDialog->LineEdit7, &QLineEdit::textChanged, this, &TMOGUIWindow::SetBMin);
    connect (pDialog->LineEdit8, &QLineEdit::textChanged, this, &TMOGUIWindow::SetRMax);
    connect (pDialog->LineEdit9, &QLineEdit::textChanged, this, &TMOGUIWindow::SetGMax);
    connect (pDialog->LineEdit10, &QLineEdit::textChanged, this, &TMOGUIWindow::SetBMax);
    connect (pDialog->ComboBox1, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIWindow::SetOperation);
		
    pDialog->ComboBox1->insertItem(0,"Top Left");
    pDialog->ComboBox1->insertItem(1,"Top Center");
    pDialog->ComboBox1->insertItem(2,"Top Right");
    pDialog->ComboBox1->insertItem(3,"Middle Left");
    pDialog->ComboBox1->insertItem(4,"Middle Center");
    pDialog->ComboBox1->insertItem(5,"Middle Right");
    pDialog->ComboBox1->insertItem(6,"Bottom Left");
    pDialog->ComboBox1->insertItem(7,"Bottom Center");
    pDialog->ComboBox1->insertItem(8,"Bottom Right");

	pDialog->LineEdit1->setText("untitled.tif");
	pDialog->LineEdit2->setText("256");
	pDialog->LineEdit3->setText("256");

	pDialog->LineEdit5->setText("0.0");
	pDialog->LineEdit6->setText("0.0");
	pDialog->LineEdit7->setText("0.0");

	pDialog->LineEdit8->setText("1.0");
	pDialog->LineEdit9->setText("1.0");
	pDialog->LineEdit10->setText("1.0");

	pColors[0] = pColors[1] = pColors[2] = 0.0;
	pColors[3] = pColors[4] = pColors[5] = 1.0;
	iFlags = 0;
	iOperation = 0;
	
    // TODO if (pDialog->exec() == QDialog::Rejected) return;

	TMOGUIImage *newfile = GetNewImage(sFileName);

	if (iFlags) iOperation = 9;

	if (newfile)
	{
		newfile->New(iWidth, iHeight, pColors, iOperation);	
		listImage.append(newfile);
		newfile->pImage->iTool = iTool;
		emit imageSelected(newfile);
		connect (newfile, SIGNAL(closeFile()), this, SLOT(closeFile()));
		connect (newfile, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
		connect (newfile->pToolsButton, SIGNAL(clicked()), this, SLOT(viewHistogram()));
		WindowChangedToolActivated(newfile);
		pMenu->SetWindows(pWorkspace);
		pTools->SetWindows(pWorkspace);
		pInfoTool->SetWindows(pWorkspace);
		pFileTool->SetWindows(pWorkspace);
		pInfo->SetOutput(newfile->pOutput);
		pProgress->hide();
	}
}

void TMOGUIWindow::NewImageConstant(bool bSelected)
{
	if (bSelected) iFlags = 1;
	else iFlags = 0;
	bZoom = bSelected;
}

void TMOGUIWindow::SetWidth(const QString& s)
{
	iWidth = s.toInt();
}

void TMOGUIWindow::SetHeight(const QString& s)
{
	iHeight = s.toInt();
}

void TMOGUIWindow::SetRMin(const QString& s)
{
	pColors[0] = s.toDouble();
}

void TMOGUIWindow::SetGMin(const QString& s)
{
	pColors[1] = s.toDouble();
}

void TMOGUIWindow::SetBMin(const QString& s)
{
	pColors[2] = s.toDouble();
}

void TMOGUIWindow::SetRMax(const QString& s)
{
	pColors[3] = s.toDouble();
}

void TMOGUIWindow::SetGMax(const QString& s)
{
	pColors[4] = s.toDouble();
}

void TMOGUIWindow::SetBMax(const QString& s)
{
	pColors[5] = s.toDouble();
}

void TMOGUIWindow::SetFileName(const QString& s)
{
	sFileName = s;
}


void TMOGUIWindow::SetOperation(int iOp)
{
	iOperation = iOp;
}

void TMOGUIWindow::pageFile()
{
    Ui::TMOGUIPageSetup *pDialog;
    QDialog *qDialog = new QDialog();
	QString s;

    pDialog = new Ui::TMOGUIPageSetup();
    pDialog->setupUi(qDialog);


    connect (pDialog->CheckBox1, &QCheckBox::toggled, this, &TMOGUIWindow::NewImageConstant);
    connect (pDialog->LineEdit1, &QLineEdit::textChanged, this, &TMOGUIWindow::SetRMin);
    connect (pDialog->LineEdit2, &QLineEdit::textChanged, this, &TMOGUIWindow::SetGMin);
    connect (pDialog->LineEdit3, &QLineEdit::textChanged, this, &TMOGUIWindow::SetBMin);
    connect (pDialog->LineEdit4, &QLineEdit::textChanged, this, &TMOGUIWindow::SetRMax);
		
	pDialog->LineEdit1->setText(s.setNum(pMargins[0]));
	pDialog->LineEdit2->setText(s.setNum(pMargins[1]));
	pDialog->LineEdit3->setText(s.setNum(pMargins[2]));
	pDialog->LineEdit4->setText(s.setNum(pMargins[3]));

	pColors[0] = pColors[1] = pColors[2] = pColors[3] = 10;
	iFlags = 0;
		
    // TODO if (pDialog->exec() == QDialog::Rejected) return;

	pMargins[0] = pColors[0];
	pMargins[1] = pColors[1];
	pMargins[2] = pColors[2];
	pMargins[3] = pColors[3];
}

void TMOGUIWindow::printFile()
{
	LoadPosition();
    /*
	QPrinter printer;
	QString fileName;
	TMOGUIImage* pImage;
	TMOImage OutImage;
	QString s;
		
    fileName = QString(pWorkspace->activeSubWindow()->widget()->objectName());
	pImage = FindImage(fileName);
	if (!pImage) return;

	printer.setFullPage(true);

    QPrintDialog dialog(printer, this);
    if (dialog.exec())
	{
        QPainter paint;
        if ( !paint.begin( &printer ) )
            return;
        pImage->pImage->DrawIt(&paint, pMargins, bZoom);
    }
    */
}

int TMOGUIWindow::SavePosition()
{
	QFile f("position.dat");
	QString temp;
    QList<int> vl;
    QList<int>::Iterator i;

	if (f.open(QIODevice::WriteOnly))
	{
        QTextStream t( &f );

		temp.setNum(x());
		t << "X = " << temp + "\n";
		temp.setNum(y());
		t << "Y = " << temp + "\n";
		temp.setNum(width());
		t << "WIDTH = " << temp + "\n";
		temp.setNum(height());
		t << "HEIGHT = " << temp + "\n";
		if (isMaximized()) t << "MAXIMIZED = ON\n";
		else t << "MAXIMIZED = OFF\n";
		if (pInfo->bVisible) t << "BOTTOM = ON\n";
		else t << "BOTTOM = OFF\n";
		vl = pSplitter->sizes();
		i = vl.begin();
		temp.setNum(*i);
		i++;
		t << "BOTTOMSPLITTER = " << temp;
		temp.setNum(*i);
		t << ", " << temp << "\n";
		if (pRight->bVisible) t << "RIGHT = ON\n";
		else t << "RIGHT = OFF\n";
        // FIXME vl = pRightSplitter->sizes();
		i = vl.begin();
		temp.setNum(*i);
		i++;
		t << "RIGHTSPLITTER = " << temp;
		temp.setNum(*i);
		t << ", " << temp << "\n";
        f.close();
		return 0;
	}
	return -1;
}

int TMOGUIWindow::LoadPosition()
{
	QFile f("position.dat");
	QString temp;
	int x = 0, y = 0, bottom = 0, right = 0, bottom1 = 0, right1 = 0;
	bool bBottom = true, bRight = true, bMaximized = false;
    QList<int> vl;

	if (f.open(QIODevice::ReadOnly))
	{
        QTextStream t( &f );
        QString s;
        while ( !t.atEnd() )
		{
            s = t.readLine();
            if (s.indexOf("WIDTH =") == 0)
			{
				iWidth = s.mid(8).toInt();
				continue;
			}
            if (s.indexOf("HEIGHT =") == 0)
			{
				iHeight = s.mid(9).toInt();
				continue;
			}
            if (s.indexOf("X =") == 0)
			{
				x = s.mid(4).toInt();
				continue;
			}
            if (s.indexOf("Y =") == 0)
			{
				y = s.mid(4).toInt();
				continue;
			}
            if (s.indexOf("RIGHT = OFF") == 0)
			{
				bRight = false;
				continue;
			}
            if (s.indexOf("BOTTOM = OFF") == 0)
			{
				bBottom = false;
				continue;
			}
            if (s.indexOf("MAXIMIZED = ON") == 0)
			{
				bMaximized = true;
				continue;
			}
            if (s.indexOf("BOTTOMSPLITTER =") == 0)
			{
                bottom = s.mid(17, s.indexOf(", ") - 17).toInt();
                bottom1 = s.mid(s.indexOf(", ") + 2).toInt();
				continue;
			}
            if (s.indexOf("RIGHTSPLITTER =") == 0)
			{
                right = s.mid(16, s.indexOf(", ") - 16).toInt();
                right1 = s.mid(s.indexOf(", ") + 2).toInt();
				continue;
			}
        }
		f.close();

        if (bMaximized){
            showMaximized();
        }else
			if (iWidth && iHeight) resize(iWidth, iHeight);
		if (x && y) move(x, y);
		pRight->bVisible = bRight;
		if (bRight) 
		{
			pRight->show();
            pMenu->SetChecked(3, 2, true);

		}
		else pRight->hide();
		pInfo->bVisible = bBottom;
		if (bBottom)
		{
			pInfo->show();
            pMenu->SetChecked(3, 1, true);
		}
		else pInfo->hide();
		if (right) 
		{
			vl.clear();
			vl.append(right);
			vl.append(right1);
			pRightSplitter->setSizes(vl);
		}
		if (bottom) 
		{
			vl.clear();
			vl.append(bottom);
			vl.append(bottom1);
			pSplitter->setSizes(vl);
		}

		return 0;
	}
	return -1;
}

void TMOGUIWindow::zoom()
{
	TMOGUIImage* pImage = GetActiveImage();
	QString s;
	int iWidth, iHeight;
	if (!pImage) return;
}

void ref_test(int& a)
{
 
}

int ref_testa()
{
 return 1;
}

void TMOGUIWindow::fitToScreen()
{
	TMOGUIImage* pImage = GetActiveImage();
	if(!pImage) return;
    pImage->fitToScreen(pWorkspace->activeSubWindow()->size()); //TODO check
	pImage->parentWidget()->move(0,0);
	pTools->SetState();
}

void TMOGUIWindow::zoomIn()
{
	TMOGUIImage* pImage = GetActiveImage();
	if(!pImage) return;
	pImage->zoomIn();
	pTools->SetState();
}

void TMOGUIWindow::zoomOut()
{
	TMOGUIImage* pImage = GetActiveImage();
	if(!pImage) return;
	pImage->zoomOut();
	pTools->SetState();
}

void TMOGUIWindow::fitToWidth()
{
	TMOGUIImage* pImage = GetActiveImage();
	if(!pImage) return;
    pImage->fitToWidth(pWorkspace->activeSubWindow()->size());
	pImage->parentWidget()->move(0,0);
	pTools->SetState();
}

void TMOGUIWindow::fitToHeight()
{
	TMOGUIImage* pImage = GetActiveImage();
	if(!pImage) return;
    pImage->fitToHeight(pWorkspace->activeSubWindow()->size());
	pImage->parentWidget()->move(0,0);
	pTools->SetState();
}

void TMOGUIWindow::activateInfoTool(bool on)
{
	TMOGUIImage* pImage = GetActiveImage();
	if(!pImage) return;
	WindowChangedToolActivated(pImage);
	iTool->SetEnabled(on);
}

void TMOGUIWindow::showToolSetting()
{
    iTool->toolDialog->show();
    //iTool->toolContext->move(this->x() + pInfoTool->x() + 4, this->y() + 78);
}

void TMOGUIWindow::WindowChangedToolActivated(TMOGUIImage * pImage)
{
	if(pInfoTool->IsActivated())
	{
		if(sPrevFileName != "")
		{
			TMOGUIImage* pImagePrev = FindImage(sPrevFileName);
			if(pImagePrev)
				pImagePrev->pImage->DeactivateTool();
		}
		pImage->pImage->ActivateTool(iTool);
        sPrevFileName = *pImage->imageName;
	}
}

void TMOGUIWindow::viewInfo()
{	
	pInfo->bVisible = !pInfo->bVisible;
	if(pInfo->bVisible)
		pInfo->show();
	else
		pInfo->hide();
    pMenu->SetChecked(3, 1, !pMenu->GetChecked(3, 1));
}


void TMOGUIWindow::viewRight()
{
	pRight->bVisible = !pRight->bVisible;
	if(pRight->bVisible)
		pRight->show();
	else
		pRight->hide();
    pMenu->SetChecked(3, 2, !pMenu->GetChecked(3, 2));//->setItemChecked(2, !pMenu->pView->isItemChecked(2));
}

void TMOGUIWindow::viewHistogram()
{
    QWidget *pWindow = pWorkspace->activeSubWindow()->widget();
	if (!pWindow) return;
    QString sName = pWindow->objectName();
	TMOGUIImage* pImage = FindImage(sName);
	if (!pImage) return;
	pImage->showtools();
    pMenu->SetChecked(3, 3, !pMenu->GetChecked(3, 3));
}

void TMOGUIWindow::refreshWindowsList(){
    for (TMOGUIImage* temp : listImage)
    {
        if (!temp->pImage || !temp->imageName || temp->imageName->isEmpty()){
            listImage.removeOne(temp);
            continue;
        }
    }
    pMenu->SetWindows(pWorkspace);
}

