#include <qapplication.h>
#include <qmessagebox.h>
#include <qassistantclient.h>
#include <qdir.h>
#include <qpushbutton.h>
#include <qhbox.h>
#include <qsplitter.h>
#include <qworkspace.h>
#include <qworkspace.h>
#include <qfiledialog.h>
#include <qcombobox.h>
#include <qcheckbox.h>
#include <qimage.h>
#include <qradiobutton.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qprinter.h>
#include <qpainter.h>
#include <qregexp.h>
#include <qwmatrix.h>
#include <qtoolbar.h>
#include <qaction.h>
#include <qtextstream.h>
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
#include "resources.h"
#include "resources1.h"
#include "resources2.h"
#include "resources3.h"
#include "resources4.h"
#include "resources6.h"

#include "TMOGUISaveDialog.h"
#include <qmap.h>

//#include <iostream>

TMOGUIResource* TMOResource::pResource = 0;

TMOGUIWindow::TMOGUIWindow( QWidget* parent, const char* name, WFlags f )
	: QMainWindow( parent, name, f )
{
	pMargins[0] = pMargins[1] = pMargins[2] = pMargins[3] = 10;
	TMOResource::pResource = new TMOGUIResource(this, "Resources");
	TMOResource::pResource->hide();
	setCaption("TMOGUI");
	setIcon(*TMOResource::pResource->IconMain->pixmap());
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
	pRightSplitter = new QSplitter(this, "RightSplitter");
	pRightSplitter->setFrameStyle( QFrame::Sunken | QFrame::Panel );
	pSplitter = new QSplitter(pRightSplitter, "BottomSplitter");
	pSplitter->setOrientation(Vertical);
	pWorkspace = new QWorkspace(pSplitter, "Workspace");
	pInfo = TMOGUIOutput::pInfo = new TMOGUIInfo(pSplitter, "Info");
	pRight = new TMOGUIRightBar(pRightSplitter, "RightBar");
	pInfo->bVisible = true;
	setCentralWidget(pRightSplitter);
	
	listImage.setAutoDelete(true);
	LoadPosition();
	pFileTool = new TMOGUIFileToolBar(this); 
	pTools = new TMOGUIZoomTool(this);
	pInfoTool = new TMOGUIInfoToolBar(this);	
	pDialog = 0;
	iTool = new TMOGUIInfoTool(this);
	//assistant = new QAssistantClient( "" );
	assistant = new QAssistantClient( QString(), 0, 0 );
    setAssistantArguments();

    connect( assistant, SIGNAL(error(const QString&)), this, SLOT(showAssistantErrors(const QString&)) );	
	connect(pRight->GetMapping()->pOk, SIGNAL(clicked()), this, SLOT(transform()));
	connect(pWorkspace, SIGNAL(windowActivated(QWidget*)), this, SLOT(windowChanged(QWidget*)));
	connect(this, SIGNAL(imageSelected(TMOGUIImage*)), pInfo->pStats, SLOT(windowChanged(TMOGUIImage*)));
	connect(this, SIGNAL(imageSelected(TMOGUIImage*)), pRight->pFilters, SLOT(windowChanged(TMOGUIImage*)));
	connect(this, SIGNAL(imageSelected(TMOGUIImage*)), pTools, SLOT(windowChanged(TMOGUIImage*)));
	connect(this, SIGNAL(imageSelected(TMOGUIImage*)), pMenu, SLOT(windowChanged(TMOGUIImage*)));
	connect(pRight, SIGNAL(closeBar()), this, SLOT(viewRight()));
	connect(pInfo, SIGNAL(closeBar()), this, SLOT(viewInfo()));
	return 0;
}

void TMOGUIWindow::showHelp()
{
    assistant->showPage( "../doc/index.html" );
}

void TMOGUIWindow::setAssistantArguments()
{
    QStringList cmdLst;
	cmdLst << "-profile" << "../doc/help.adp";
    assistant->setArguments( cmdLst );
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
	TMOGUIImage* newfile = 0, *temp;
	int iCount = 0;
	
        if ( fileName.isEmpty() ) return;

	for (temp = listImage.first(); temp; temp = listImage.next())
	{
		name = temp->name();
		if (name.contains(fileName)) iCount++;
	}
	
	if (iCount)
	{
		number.setNum(iCount);
		newfile = new TMOGUIImage(pProgress, pWorkspace, fileName + " [" + number + "]");
	}
	else
		newfile = new TMOGUIImage(pProgress, pWorkspace, fileName);

	if (newfile->Open(fileName)) delete newfile;
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
		( QString::null, "HDR images (*.raw *.hdrraw *.tif *.tiff *.pic *.hdr *.exr *.jpg *.jpeg);; *.raw;; *.hdrraw;; *.tif;; *.tiff;; *.pic;; *.hdr;; *.exr;; *.jpg;; *.jpeg;; All files (*.*)", this )
        	);
}

void TMOGUIWindow::saveasFile()
{
	QString fileName;
	TMOGUIImage* pImage;
	TMOImage OutImage;
	QString s;
		
	fileName = QString(pWorkspace->activeWindow()->name());
	pImage = FindImage(fileName);
	if (!pImage) return;

	int iIdx = 0, iFound;
	if ((iFound = fileName.find("[")) > 0) 
	{
		QString sTemp = fileName.mid(iFound + 1);
		sTemp = sTemp.left(sTemp.find("]"));
		iIdx = sTemp.toInt();
		fileName = fileName.left(iFound);
	}
	
	if ((iFound = fileName.findRev(".", -1)) > 0) 
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
	TMOGUISaveDialog * fd = new TMOGUISaveDialog(fileName,&filter,this,"",true);
	fd->setSelectedFilter(4);
	if(fd->exec() == QDialog::Accepted)
        {
	 fileName = fd->selectedFile();
	 OutImage.New(*pImage->GetImage());
	 pImage->pImage->GetImage(&OutImage);
	 if(OutImage.SaveAs(fileName, fd->getSelectedFileType()))
	  QMessageBox::critical( 0, "Tone Mapping Studio", QString("Failed to save file : \n\n") + fileName + "\n(probably unknown format)");
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

	fileName = pWorkspace->activeWindow()->name();
	pImage = FindImage(fileName);
	if (!pImage) return;

	int iIdx = 0;
	if ((iFound = fileName.find("[")) > 0) 
	{
		QString sTemp = fileName.mid(iFound + 1);
		sTemp = sTemp.left(sTemp.find("]"));
		iIdx = sTemp.toInt();
		fileName = fileName.left(iFound);
	}
	
	if ((iFound = fileName.findRev(".", -1)) > 0) 
	{
		QString sTemp;
		if (iIdx) fileName = fileName.left(iFound) + "[" + sTemp.setNum(iIdx) + "].tif";
		else fileName = fileName.left(iFound) + ".tif";
	}

	OutImage.New(*pImage->GetImage());
	pImage->pImage->GetImage(&OutImage);

	OutImage.SaveAs(fileName);
	if (fileName == pWorkspace->activeWindow()->name())
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
	TMOGUIImage* newfile = 0, *temp;
	QString number,fileName = pMenu->GetRecent(ID), name;
	
    if ( fileName.isEmpty() ) return;

	for (temp = listImage.first(); temp; temp = listImage.next())
	{
		name = temp->name();
		if (name.contains(fileName)) iCount++;
	}
	
	if (iCount)
	{
		number.setNum(iCount);
		newfile = new TMOGUIImage(pProgress, pWorkspace, fileName + " [" + number + "]");
	}
	else
		newfile = new TMOGUIImage(pProgress, pWorkspace, fileName);

	if (newfile->Open(fileName)) delete newfile;
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

void TMOGUIWindow::closeFile()
{
	QString s;
	TMOGUIImage* im;

	s = pWorkspace->activeWindow()->name();
	im = FindImage(s);
	if (im->Terminate()) return;
	listImage.remove(im);
	
	pTools->SetWindows(pWorkspace);
	pInfoTool->SetWindows(pWorkspace);
	pFileTool->SetWindows(pWorkspace);
	pMenu->SetWindows(pWorkspace);
}

void TMOGUIWindow::closeallWindow()
{
	listImage.clear();

	pMenu->SetWindows(pWorkspace);
	pTools->SetWindows(pWorkspace);
	pInfo->SetOutput(0);
}

void TMOGUIWindow::saveallFile()
{
	QString fileName;
	TMOGUIImage *pImage;
	TMOImage OutImage;
	int iFound;

	for (pImage = listImage.first(); pImage; pImage = listImage.next())
	{
		fileName = pImage->name();
		if (!pImage) return;

		if ((iFound = fileName.find("[")) > 0) 
		{
			fileName = fileName.left(iFound);
		}
		
		if ((iFound = fileName.findRev(".", -1)) > 0) 
		{
			fileName = fileName.left(iFound) + ".tif";
		}

		OutImage.New(*pImage->GetImage());
		pImage->pImage->GetImage(&OutImage);

		OutImage.SaveAs(fileName);
		if (fileName == pWorkspace->activeWindow()->name())
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
	TMOGUIImage* retval = 0;

	for (retval = listImage.first(); retval; retval = listImage.next())
	{
		if (retval->name() == name) break;
	}
	return retval;
}

void TMOGUIWindow::activateWindow(int id)
{
	QWidgetList wl; 
	QWidget* widget;
	QString sName;
	TMOGUIImage *pImage;
	
	int number = 64;

	wl = pWorkspace->windowList();

	for (widget = wl.first(); widget; widget = wl.next())
	{
		if (number++ == id)
		{
			widget->setFocus();			
			sName = widget->name();
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
	TMOImage *pSrc = 0, *pDst = 0;
	TMOGUIImage *pImage = 0;
	QWidget *pWidget;
	TMOGUITransformation *pTransform;
	TMO* pTMO;
	QString sName;

	pTMO = pRight->GetTMO();
	if (!pTMO) return;
	pWidget = pWorkspace->activeWindow();
	if (!pWidget) return;
	sName = pWidget->name();
	pImage = FindImage(sName);
	if (!pImage) return;

	pTransform = pImage->Transform();
	if (pTransform)
	{
		pSrc = pImage->GetImage();
		pTMO->SetSource(*pSrc);
		pTransform->SetTMO(pTMO);
	}
}

void TMOGUIWindow::finishTransform()
{
	TMOGUIImage *pImage = 0;
	QString sName;

	sName = pWorkspace->activeWindow()->name();
	pImage = FindImage(sName);
	if (pImage->CanUndo()) pMenu->Enable(2, 1);
	else pMenu->Disable(2, 1);
}

void TMOGUIWindow::undoEdit()
{
	TMOGUIImage *pImage = 0;
	QString sName;

	if(!pWorkspace || !(pWorkspace->activeWindow()))return;
	sName = pWorkspace->activeWindow()->name();
	pImage = FindImage(sName);
	if (!pImage) return;

	pImage->Swap();
}

void TMOGUIWindow::windowChanged(QWidget* pWidget)
{
	if (!pWidget) 
	{
		emit imageSelected(0);
	}
	else
	{
		QString sName = pWidget->name();
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
	QWidget *pWindow = pWorkspace->activeWindow();
	if (!pWindow) return;
	QString sName = pWindow->name();
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
	TMOGUIImage* newfile = 0, *temp;
	QString number,fileName = sName, name;
	
    if ( fileName.isEmpty() ) return 0;
	
	QRegExp r = QRegExp(" \\[[0123456789]+\\]", TRUE, FALSE);
	int iFind = fileName.find(r);
	fileName = fileName.left(iFind);
    
	for (temp = listImage.first(); temp; temp = listImage.next())
	{
		name = temp->name();

		if (name.contains(fileName)>=0) iCount++;
	}

	if (iCount)
	{
		number.setNum(iCount);
		newfile = new TMOGUIImage(pProgress, pWorkspace, fileName + " [" + number + "]");
	}
	else
		newfile = new TMOGUIImage(pProgress, pWorkspace, fileName);

	return newfile;
}

void TMOGUIWindow::sizeCommand()
{
	TMOGUIImage* pImage = GetActiveImage();
	QString s;
	int iWidth, iHeight;
	if (!pImage) return;
	
	pDialog = new TMOGUIImageSize(this, "ImageSizeDialog", true);
	iFlags = iFlags & ~1;
	pDialog->ratioLabel->setHidden(true);
	dRatio = (double)pImage->GetImage()->GetWidth() / pImage->GetImage()->GetHeight();
	pDialog->LineEdit1->setText(s.setNum(pImage->GetImage()->GetWidth()));
	pDialog->LineEdit2->setText(s.setNum(pImage->GetImage()->GetHeight()));
	connect (pDialog->LineEdit1, SIGNAL(textChanged(const QString&)), this, SLOT(ImageSizeWidth(const QString &)));
	connect (pDialog->LineEdit2, SIGNAL(textChanged(const QString&)), this, SLOT(ImageSizeHeight(const QString &)));
	connect (pDialog->CheckBox1, SIGNAL(toggled(bool)), this, SLOT(ImageSizeConstrain(bool)));
	connect (this, SIGNAL(signalImageSizeWidth(const QString&)), pDialog->LineEdit1, SLOT(setText(const QString &)));
	connect (this, SIGNAL(signalImageSizeHeight(const QString&)), pDialog->LineEdit2, SLOT(setText(const QString &)));
	if (pDialog->exec() == QDialog::Rejected) return;
	iWidth = pDialog->LineEdit1->text().toInt();
	iHeight = pDialog->LineEdit2->text().toInt();
	pImage->SetImageSize(iWidth, iHeight);	
	delete pDialog;
	pDialog = 0;
}

TMOGUIImage* TMOGUIWindow::GetActiveImage()
{
	QWidget *pWindow = pWorkspace->activeWindow();
	if (!pWindow) return 0;
	QString sName = pWindow->name();
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
	TMOGUIImage *newfile = GetNewImage(pImage->name());

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

void TMOGUIWindow::extractComCommand(int iComponent)
{
	TMOGUIImage* pImage = GetActiveImage();
	if (!pImage) return;
	TMOGUIImage *newfile = GetNewImage(pImage->name());

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
	TMOGUIMergeComponents *pDialog;
	TMOGUIImage *temp;
	
	pDialog = new TMOGUIMergeComponents(this, "MergeComponentsDialog", true);
	connect (pDialog->ComboBox1, SIGNAL(activated(int)), this, SLOT(MergeComponentsRed(int)));
	connect (pDialog->ComboBox2, SIGNAL(activated(int)), this, SLOT(MergeComponentsGreen(int)));
	connect (pDialog->ComboBox3, SIGNAL(activated(int)), this, SLOT(MergeComponentsBlue(int)));
	connect (this, SIGNAL(signalMergeComRed(const QPixmap &)), pDialog->PixmapLabel1, SLOT(setPixmap(const QPixmap &)));
	connect (this, SIGNAL(signalMergeComGreen(const QPixmap &)), pDialog->PixmapLabel2, SLOT(setPixmap(const QPixmap &)));
	connect (this, SIGNAL(signalMergeComBlue(const QPixmap &)), pDialog->PixmapLabel3, SLOT(setPixmap(const QPixmap &)));
		
	QString s;
	int i = 0;
	for (temp = listImage.first(); temp; temp = listImage.next(), i++)
	{
		s = TMOGUIImage::GetName(temp->name());
		pDialog->ComboBox1->insertItem(s, i);
		pDialog->ComboBox2->insertItem(s, i);
		pDialog->ComboBox3->insertItem(s, i);
	}

	MergeComponentsRed(0);
	MergeComponentsGreen(0);
	MergeComponentsBlue(0);

	s = "";
	if (pImages[0]) s.append(TMOGUIImage::GetName(pImages[0]->name()) + "_");
	if (pImages[1]) s.append(TMOGUIImage::GetName(pImages[1]->name()) + "_");
	if (pImages[2]) s.append(TMOGUIImage::GetName(pImages[2]->name()) + "_");
	s = s.left(s.length());

	if (pDialog->exec() == QDialog::Rejected) return;

	TMOGUIImage *newfile = new TMOGUIImage(pProgress, pWorkspace, s);

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
	
	temp = listImage.first();
	for (int i = iImage; temp && i; i--) temp = listImage.next();

	pImages[0] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
	QImage tempImage(160, 160, 32);
    QWMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->xForm(m);
	emit signalMergeComRed(aPixmap);
}

void TMOGUIWindow::OperationSecond(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	
	temp = listImage.first();
	for (int i = iImage; temp && i; i--) temp = listImage.next();

	pImages[1] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
	QImage tempImage(160, 160, 32);
    QWMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->xForm(m);
	emit signalMergeComGreen(aPixmap);
}

void TMOGUIWindow::MergeComponentsRed(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	QRgb* pRgb;
	
	temp = listImage.first();
	for (int i = iImage; temp && i; i--) temp = listImage.next();

	pImages[0] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
	QImage tempImage(160, 160, 32);
    QWMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->xForm(m);
	
	tempImage = aPixmap.convertToImage();
	for (int jx = 0; jx < aPixmap.height(); jx++)
	{
		pRgb = (QRgb*)tempImage.scanLine(jx);
		for (int ix = 0; ix < aPixmap.width(); ix++, pRgb++)
		{
			*pRgb = qRgb( qRed(*pRgb), 0 , 0);
		}
	}
	aPixmap.convertFromImage(tempImage, ColorOnly);
	
	emit signalMergeComRed(aPixmap);
}

void TMOGUIWindow::MergeComponentsGreen(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	QRgb* pRgb;
	
	temp = listImage.first();
	for (int i = iImage; temp && i; i--) temp = listImage.next();

	pImages[1] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
	QImage tempImage(160, 160, 32);
    QWMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->xForm(m);
	
	tempImage = aPixmap.convertToImage();
	for (int jx = 0; jx < aPixmap.height(); jx++)
	{
		pRgb = (QRgb*)tempImage.scanLine(jx);
		for (int ix = 0; ix < aPixmap.width(); ix++, pRgb++)
		{
			*pRgb = qRgb( 0, qGreen(*pRgb), 0);
		}
	}
	aPixmap.convertFromImage(tempImage, ColorOnly);
	emit signalMergeComGreen(aPixmap);
}

void TMOGUIWindow::MergeComponentsBlue(int iImage)
{
	TMOGUIImage *temp;
	QPixmap aPixmap;
	QRgb* pRgb;
	
	temp = listImage.first();
	for (int i = iImage; temp && i; i--) temp = listImage.next();

	pImages[2] = temp;
	QPixmap *tempPixmap = temp->pImage->pSrcPixmap;
	QImage tempImage(160, 160, 32);
    QWMatrix m;
	double aspect;
	if (tempPixmap->height() > tempPixmap->width())
		aspect = 160.0 / tempPixmap->height();
	else
		aspect = 160.0 / tempPixmap->width();

    m.scale(aspect, aspect);
    aPixmap = tempPixmap->xForm(m);
	
	tempImage = aPixmap.convertToImage();
	for (int jx = 0; jx < aPixmap.height(); jx++)
	{
		pRgb = (QRgb*)tempImage.scanLine(jx);
		for (int ix = 0; ix < aPixmap.width(); ix++, pRgb++)
		{
			*pRgb = qRgb( 0, 0, qBlue(*pRgb));
		}
	}
	aPixmap.convertFromImage(tempImage, ColorOnly);
	emit signalMergeComBlue(aPixmap);
}

void TMOGUIWindow::operationCommand()
{
	TMOGUIOperation *pDialog;
	TMOGUIImage *temp;
		
	pDialog = new TMOGUIOperation(this, "OperationsDialog", true);
	connect (pDialog->ComboBox1, SIGNAL(activated(int)), this, SLOT(OperationFirst(int)));
	connect (pDialog->ComboBox2, SIGNAL(activated(int)), this, SLOT(OperationSecond(int)));
	connect (pDialog->ComboBox3, SIGNAL(activated(int)), this, SLOT(ImageOperation(int)));
	connect (this, SIGNAL(signalMergeComRed(const QPixmap &)), pDialog->PixmapLabel1, SLOT(setPixmap(const QPixmap &)));
	connect (this, SIGNAL(signalMergeComGreen(const QPixmap &)), pDialog->PixmapLabel2, SLOT(setPixmap(const QPixmap &)));
		
	QString s;
	int i = 0;
	for (temp = listImage.first(); temp; temp = listImage.next(), i++)
	{
		s = TMOGUIImage::GetName(temp->name());
		pDialog->ComboBox1->insertItem(s, i);
		pDialog->ComboBox2->insertItem(s, i);
	}

	pDialog->ComboBox3->insertItem("Addition", 0);
	pDialog->ComboBox3->insertItem("Subtraction", 1);
	pDialog->ComboBox3->insertItem("Multiplication", 2);
	pDialog->ComboBox3->insertItem("Division", 3);
	pDialog->ComboBox3->insertItem("Involution", 4);
	pDialog->ComboBox3->insertItem("Maximum", 5);
	pDialog->ComboBox3->insertItem("Minimum", 6);
	pDialog->ComboBox3->insertItem("Average", 7);
	
	OperationFirst(0);
	OperationSecond(0);
	iFlags = 0;
	
	s = "";
	if (pImages[0]) s.append(TMOGUIImage::GetName(pImages[0]->name()) + "_");
	if (pImages[1]) s.append(TMOGUIImage::GetName(pImages[1]->name()));
	s = s.left(s.length());

	if (pDialog->exec() == QDialog::Rejected) return;

	TMOGUIImage *newfile = new TMOGUIImage(pProgress, pWorkspace, s);

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
	TMOGUINewFile *pDialog;
	QString s;
	
	pDialog = new TMOGUINewFile(this, "NewFileDialog", true);

	connect (pDialog->RadioButton1, SIGNAL(toggled(bool)), this, SLOT(NewImageConstant(bool)));
	connect (pDialog->LineEdit1, SIGNAL(textChanged(const QString &)), this, SLOT(SetFileName(const QString &)));
	connect (pDialog->LineEdit2, SIGNAL(textChanged(const QString &)), this, SLOT(SetWidth(const QString &)));
	connect (pDialog->LineEdit3, SIGNAL(textChanged(const QString &)), this, SLOT(SetHeight(const QString &)));
	connect (pDialog->LineEdit5, SIGNAL(textChanged(const QString &)), this, SLOT(SetRMin(const QString &)));
	connect (pDialog->LineEdit6, SIGNAL(textChanged(const QString &)), this, SLOT(SetGMin(const QString &)));
	connect (pDialog->LineEdit7, SIGNAL(textChanged(const QString &)), this, SLOT(SetBMin(const QString &)));
	connect (pDialog->LineEdit8, SIGNAL(textChanged(const QString &)), this, SLOT(SetRMax(const QString &)));
	connect (pDialog->LineEdit9, SIGNAL(textChanged(const QString &)), this, SLOT(SetGMax(const QString &)));
	connect (pDialog->LineEdit10, SIGNAL(textChanged(const QString &)), this, SLOT(SetBMax(const QString &)));
	connect (pDialog->ComboBox1, SIGNAL(activated(int)), this, SLOT(SetOperation(int)));
		
	pDialog->ComboBox1->insertItem("Top Left", 0);
	pDialog->ComboBox1->insertItem("Top Center", 1);
	pDialog->ComboBox1->insertItem("Top Right", 2);
	pDialog->ComboBox1->insertItem("Middle Left", 3);
	pDialog->ComboBox1->insertItem("Middle Center", 4);
	pDialog->ComboBox1->insertItem("Middle Right", 5);
	pDialog->ComboBox1->insertItem("Bottom Left", 6);
	pDialog->ComboBox1->insertItem("Bottom Center", 7);
	pDialog->ComboBox1->insertItem("Bottom Right", 8);

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
	
	if (pDialog->exec() == QDialog::Rejected) return;

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
	TMOGUIPageSetup *pDialog;
	QString s;
	
	pDialog = new TMOGUIPageSetup(this, "PageSetupDialog", true);

	connect (pDialog->CheckBox1, SIGNAL(toggled(bool)), this, SLOT(NewImageConstant(bool)));
	connect (pDialog->LineEdit1, SIGNAL(textChanged(const QString &)), this, SLOT(SetRMin(const QString &)));
	connect (pDialog->LineEdit2, SIGNAL(textChanged(const QString &)), this, SLOT(SetGMin(const QString &)));
	connect (pDialog->LineEdit3, SIGNAL(textChanged(const QString &)), this, SLOT(SetBMin(const QString &)));
	connect (pDialog->LineEdit4, SIGNAL(textChanged(const QString &)), this, SLOT(SetRMax(const QString &)));
		
	pDialog->LineEdit1->setText(s.setNum(pMargins[0]));
	pDialog->LineEdit2->setText(s.setNum(pMargins[1]));
	pDialog->LineEdit3->setText(s.setNum(pMargins[2]));
	pDialog->LineEdit4->setText(s.setNum(pMargins[3]));

	pColors[0] = pColors[1] = pColors[2] = pColors[3] = 10;
	iFlags = 0;
		
	if (pDialog->exec() == QDialog::Rejected) return;

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
		
	fileName = QString(pWorkspace->activeWindow()->name());
	pImage = FindImage(fileName);
	if (!pImage) return;

	printer.setFullPage(true);
    if ( printer.setup( this ) ) 
	{
        QPainter paint;
        if ( !paint.begin( &printer ) )
            return;
        pImage->pImage->DrawIt(&paint, pMargins, bZoom);
    }*/
}

int TMOGUIWindow::SavePosition()
{
	QFile f("position.dat");
	QString temp;
	QValueList<int> vl;
	QValueList<int>::Iterator i;

	if (f.open(IO_WriteOnly))
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
		vl = pRightSplitter->sizes();
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
	QValueList<int> vl;

	if (f.open(IO_ReadOnly))
	{
        QTextStream t( &f );        
        QString s;
        while ( !t.eof() )			
		{
            s = t.readLine();
			if (s.find("WIDTH =") == 0) 
			{
				iWidth = s.mid(8).toInt();
				continue;
			}
			if (s.find("HEIGHT =") == 0) 
			{
				iHeight = s.mid(9).toInt();
				continue;
			}
			if (s.find("X =") == 0) 
			{
				x = s.mid(4).toInt();
				continue;
			}
			if (s.find("Y =") == 0) 
			{
				y = s.mid(4).toInt();
				continue;
			}
			if (s.find("RIGHT = OFF") == 0) 
			{
				bRight = false;
				continue;
			}
			if (s.find("BOTTOM = OFF") == 0) 
			{
				bBottom = false;
				continue;
			}
			if (s.find("MAXIMIZED = ON") == 0) 
			{
				bMaximized = true;
				continue;
			}
			if (s.find("BOTTOMSPLITTER =") == 0) 
			{
				bottom = s.mid(17, s.find(", ") - 17).toInt();
				bottom1 = s.mid(s.find(", ") + 2).toInt();
				continue;
			}
			if (s.find("RIGHTSPLITTER =") == 0) 
			{
				right = s.mid(16, s.find(", ") - 16).toInt();
				right1 = s.mid(s.find(", ") + 2).toInt();
				continue;
			}
        }
		f.close();

		if (bMaximized)
			showMaximized();
		else
			if (iWidth && iHeight) resize(iWidth, iHeight);
		if (x && y) move(x, y);
		pRight->bVisible = bRight;
		if (bRight) 
		{
			pRight->show();
			pMenu->pView->setItemChecked(2, true);
		}
		else pRight->hide();
		pInfo->bVisible = bBottom;
		if (bBottom)
		{
			pInfo->show();
			pMenu->pView->setItemChecked(1, true);
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
	pImage->fitToScreen(pWorkspace->size());
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
	pImage->fitToWidth(pWorkspace->size());
	pImage->parentWidget()->move(0,0);
	pTools->SetState();
}

void TMOGUIWindow::fitToHeight()
{
	TMOGUIImage* pImage = GetActiveImage();
	if(!pImage) return;
	pImage->fitToHeight(pWorkspace->size());
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
	iTool->toolContext->show();
	iTool->toolContext->move(this->x() + pInfoTool->x() + 4, this->y() + 78);
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
		sPrevFileName = pImage->name();
	}
}

void TMOGUIWindow::viewInfo()
{	
	pInfo->bVisible = !pInfo->bVisible;
	if(pInfo->bVisible)
		pInfo->show();
	else
		pInfo->hide();
	pMenu->pView->setItemChecked(1, !pMenu->pView->isItemChecked(1));
}


void TMOGUIWindow::viewRight()
{
	pRight->bVisible = !pRight->bVisible;
	if(pRight->bVisible)
		pRight->show();
	else
		pRight->hide();
	pMenu->pView->setItemChecked(2, !pMenu->pView->isItemChecked(2));		
}

void TMOGUIWindow::viewHistogram()
{
	QWidget *pWindow = pWorkspace->activeWindow();
	if (!pWindow) return;
	QString sName = pWindow->name();
	TMOGUIImage* pImage = FindImage(sName);
	if (!pImage) return;
	pImage->showtools();
	pMenu->pView->setItemChecked(3, !pMenu->pView->isItemChecked(3));
}
