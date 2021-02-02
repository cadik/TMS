// TMOGUIToneMappingChooser.cpp: implementation of the TMOGUIToneMappingChooser class that displays preview of all TMO.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"
#include "TMOGUIToneMappingChooser.h"
#include "TMOGUIToneMapping.h"
#include "TMOGUIPreviewImage.h"
#include "TMOGUITransformation.h"
#include "TMOGUIBitmap.h"
#include <qpushbutton.h>
#include <qlabel.h>
#include <QResizeEvent>
#include <qimage.h>
#include <qpixmap.h>
#include <qtemporarydir.h>
#include <qfile.h>
#include <qsizepolicy.h>
#include <QFileInfo>
#include <QStandardPaths>

#include "lqstring.h"

TMOGUIToneMappingChooser::TMOGUIToneMappingChooser(QWidget *parent, TMOGUIImage *pImg) : QScrollArea(parent)
{
    //setAutoFillBackground(true);

    imageWidgetSize = QSize(300, 250);
    imageFitSize = QSize(290,240);
    listPreview.clear();
    parentWidget = (TMOGUIToneMapping*) parent;
    this->setWidgetResizable(true);
    this->setAlignment(Qt::AlignTop);
    //this->setBackgroundRole(QPalette::Light);
    big_box = new QWidget(this->viewport());
    bigBoxLayout = new QVBoxLayout(big_box);
    bigBoxLayout->setAlignment(Qt::AlignTop);

    big_box->setLayout(bigBoxLayout);
    //bigBoxLayout->addStretch();
    this->setWidget(big_box);

    if(!pImg){
        QString previewFileName = "preview-src.tif";
        imgName = previewFileName;
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
    } else {
        pPreviewImage = pImg;
        imgName = pImg->GetImage()->GetFilename();
    }
    pDefaultImage = pPreviewImage;
    //displayAll();
}

TMOGUIToneMappingChooser::~TMOGUIToneMappingChooser()
{
    //if(!listPreview.empty()) listPreview.clear();
    if(!widgetList.empty()) widgetList.clear();
    //if(pPreviewImage) delete pPreviewImage;
}

void TMOGUIToneMappingChooser::AddTMOPreview(TMO* pTMO, int indexLib, int indexTMO)
{

    lqstring s;
    QString name;



    try
    {
        name = s.setUnicodeCodes(pTMO->GetName(), wcslen(pTMO->GetName()));
    }
    catch(...)
    {
        // bad plugin format // plugin version mismatch
    }

    QPushButton* btn = new QPushButton(name, this);
    TMOGUIPreviewImage *widget = new TMOGUIPreviewImage(big_box, name, btn, pPreviewImage, pTMO);

    widget->getImage()->fitToScreen(widget->getImageViewSize());

    // TODO send widget

    connect(widget, SIGNAL(finishTransform()), this, SLOT(finishTransform()));
    connect(btn, &QPushButton::released, [=]{
        parentWidget->chooseTechnique(indexLib, indexTMO);
    });

    widgetList.append(widget);
    bigBoxLayout->addWidget(widget, Qt::AlignTop);


}

void TMOGUIToneMappingChooser::displayAll()
{
    foreach(TMOGUIPreviewImage* w, widgetList){
        w->display();
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
    //SavePreview(listPreviewIt.key());
}

void TMOGUIToneMappingChooser::windowChanged(TMOGUIImage* pImage){

    if(pImage){

        if(!pImage->bPreview){
            pPreviewImage->NewSmall(pImage);
        } else {
            pPreviewImage = pDefaultImage;
        }

    }
    foreach(TMOGUIPreviewImage* w, widgetList){
        w->setImage(pPreviewImage);
    }

    displayAll();
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



