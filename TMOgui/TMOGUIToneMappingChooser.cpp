// TMOGUIToneMappingChooser.cpp: implementation of the TMOGUIToneMappingChooser class that displays preview of all TMO.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"
#include "TMOGUIToneMappingChooser.h"
#include "TMOGUIToneMapping.h"
#include "TMOGUIImage.h"
#include "TMOGUITransformation.h"
#include "TMOGUIBitmap.h"
#include <qpushbutton.h>
#include <qlabel.h>
#include <QResizeEvent>
#include <qimage.h>
#include <qpixmap.h>
#include <qtemporarydir.h>
#include <qfile.h>
#include <QFileInfo>
#include <QStandardPaths>

#include "lqstring.h"

TMOGUIToneMappingChooser::TMOGUIToneMappingChooser(QWidget *parent) : QScrollArea(parent)
{
    imageWidgetSize = QSize(300, 250);
    imageFitSize = QSize(290,240);
    listPreview.clear();
    parentWidget = (TMOGUIToneMapping*) parent;
    this->setWidgetResizable(true);
    this->setAlignment(Qt::AlignTop);
    big_box = new QWidget(this->viewport()); //Q3VBox
    bigBoxLayout = new QVBoxLayout(big_box);
    bigBoxLayout->setAlignment(Qt::AlignTop);

    big_box->setLayout(bigBoxLayout);
    //bigBoxLayout->addStretch();
    this->setWidget(big_box);


    QString previewFileName = "preview-src.tif";
    path = new QString(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation));
    //QString absolutePath = QStandardPaths::locate(QStandardPaths::AppDataLocation, previewFileName);
    QDir dir(*path);
    if (!dir.exists())
        dir.mkpath(*path);

    QString absolutePath = dir.absoluteFilePath(previewFileName);
    if (QFile::exists(absolutePath))
        QFile::remove(absolutePath);


    QFile::copy(":/resources/images/preview-source.tif", absolutePath);
    pPreviewImage = new TMOGUIImage(nullptr, this, "Preview-src");
    pPreviewImage->Open(absolutePath.toStdString().c_str());
    pPreviewImage->hide();
}

TMOGUIToneMappingChooser::~TMOGUIToneMappingChooser()
{
    if(pPreviewImage) delete pPreviewImage;
    if(!listPreview.empty()) listPreview.clear();
}

void TMOGUIToneMappingChooser::AddTMOPreview(TMO* pTMO, int indexLib, int indexTMO)
{

    lqstring s;
    QString name;

    TMOGUIImage* pImage = nullptr;

    try
    {
        name = s.setUnicodeCodes(pTMO->GetName(), wcslen(pTMO->GetName()));
    }
    catch(...)
    {
        // bad plugin format // plugin version mismatch
    }

    QWidget *widget = new QWidget(big_box);
    widget->setFixedSize(imageWidgetSize);
    QVBoxLayout *vLayout = new QVBoxLayout(widget);
    //QLabel *lbl = new QLabel(s, widget);
    QPushButton* btn = new QPushButton(s, widget);

    pImage = new TMOGUIImage(nullptr, widget, name.toStdString().c_str());
    if(pImage->New(pPreviewImage)) return;
    //pImage->setFixedHeight(pPreviewImage->GetImage()->GetHeight());
    connect(pImage, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
    connect(btn, &QPushButton::released, [=]{
        parentWidget->chooseTechnique(indexLib, indexTMO);
    });

    //vLayout->addWidget(lbl);
    vLayout->addWidget(pImage);
    vLayout->addWidget(btn);

    pImage->hideAll(true);
    pImage->fitToWidth(imageFitSize);

    // FIXME Son14 not working
    if(name != "Son14") listPreview.insert(pImage, pTMO);
    bigBoxLayout->addWidget(widget, Qt::AlignTop);
}

void TMOGUIToneMappingChooser::display()
{
    TMOImage *pSrc = nullptr;
    TMOGUITransformation *pTransform = nullptr;
    QMap<TMOGUIImage*, TMO*>::Iterator i;
    TMOGUIImage *pImage = nullptr;
    TMO *pTMO = nullptr;

    listPreviewIt = listPreview.begin();
    if(listPreviewIt == listPreview.end()) return;

    pImage = listPreviewIt.key();
    pTMO = listPreviewIt.value();


    if(!openPreview(pImage)){
        pSrc = pImage->GetImage();
        //pImage->setFixedSize(pSrc->GetWidth(), pSrc->GetHeight());
        pImage->fitToWidth(imageFitSize);
        pImage->hideAll(true);
        delete listPreviewIt.value();
        listPreview.remove(pImage);
        display();
        return;
    }

    pTransform = pImage->Transform();
    if (pTransform)
    {
        pSrc = pImage->GetImage();
        pTMO->SetImage(*pSrc);

        //pImage->setFixedSize(pSrc->GetWidth(), pSrc->GetHeight());
        pImage->fitToWidth(imageFitSize);
        pTransform->SetTMO(pTMO);
    }

}
void TMOGUIToneMappingChooser::viewportResizeEvent ( QResizeEvent * e )
{
    QSize s = e->size();
    backWidth = s.width() - 63;
    /*for(TMOGUIImage* pImage : listPreview){
        pImage->fitToWidth(s);
    }*/
}

void TMOGUIToneMappingChooser::finishTransform()
{
    savePreview(listPreviewIt.key());
    if(!listPreview.empty())
    {
        display();
    }
}

int TMOGUIToneMappingChooser::savePreview(TMOGUIImage* pImage)
{
    QString fileName = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/" + pImage->objectName() + ".tif";
    TMOImage OutImage;
    QString s;

    OutImage.New(*pImage->GetImage());
    pImage->pImage->GetImage(&OutImage);
    return OutImage.SaveAs(fileName.toStdString().c_str(), TMO_TIFF_8);
}

int TMOGUIToneMappingChooser::openPreview(TMOGUIImage* pImage)
{
    QString previewPath = pImage->objectName() + ".tif";
    QString absolutePath = QStandardPaths::locate(QStandardPaths::AppDataLocation, previewPath);
    if(!absolutePath.isEmpty()){
        return pImage->Open(absolutePath.toStdString().c_str());
    } else {
        return -1;
    }
}



