#include "../tmolib/TMO.h"
#include "TMOGUIPreviewImage.h"
#include "TMOGUIImage.h"
#include "TMOGUITransformation.h"
#include "TMOGUIStatus.h"


TMOGUIPreviewImage::TMOGUIPreviewImage(QWidget *parent, QString name, QPushButton *btn, TMOGUIImage *image, TMO* prTMO) : QWidget(parent)
{
    QVBoxLayout *vLayout = new QVBoxLayout(this);

    pTMO = prTMO;
    bTransforming = false;
    bTransformed = false;
    TMOname = new QString(name);
    QString origImageName = *image->imageName;
    imageName = new QString(origImageName.append("-").append(name));
    TMOGUIStatus* pStatus = new TMOGUIStatus(this, "Status");
    pBar = new TMOGUIProgressBar(pStatus, name.toStdString().c_str());
    pImage = new TMOGUIImage(pBar, this, imageName->toStdString().c_str(), true);
    pDefaultImage = pImage;
    if(pImage->NewSmall(image)) return; //TODO check SMALL too slow and size not exact
    pStatus->setHidden(true);
    //pImage->setFixedHeight(image->GetImage()->GetHeight());
    pImage->hideAll(true);
    imageFitSize = pImage->GetViewSize();//QSize(400,240);
    vLayout->addWidget(pImage);

    loadingLbl = new QLabel("Loading...", pImage);
    //loadingLbl->setFixedHeight(80);
    loadingLbl->setAlignment(Qt::AlignCenter);
    loadingLbl->setStyleSheet("background: rgba(255, 255, 255, 0.6)");
    //loadingLbl->setAttribute(Qt::WA_TranslucentBackground);
    loadingLbl->setGeometry(pImage->GetViewGeometry());


    //vLayout->addWidget(loadingLbl);
    vLayout->addWidget(btn);

    connect(pImage, SIGNAL(finishTransform()), this, SLOT(finishImageTransform()));

}

TMOGUIPreviewImage::~TMOGUIPreviewImage()
{
    //TODO if(pTransformation) delete pTransformation;
}

TMOGUIImage* TMOGUIPreviewImage::getImage(){
    return pImage;
}

void TMOGUIPreviewImage::setImage(TMOGUIImage* srcImage){

    bTransformed = false;
    loadingLbl->setHidden(false);
    //pImage->setHidden(true);
    if(pTransformation && bTransforming){
        pTransformation->Cancel();
        bTransforming = false;
    }

    if(srcImage){
        pImage->New(srcImage);
    } else {
        pImage->New(pDefaultImage);
    }

    pImage->hideAll(true);



}

QSize TMOGUIPreviewImage::getImageViewSize()
{
    return pImage->GetViewSize();
}

int TMOGUIPreviewImage::display(){
    // TODO add loading picture
    // TODO save preview
    if(bTransformed) return 0;

    TMOImage *pSrc = nullptr;
    if(pTransformation) pTransformation = nullptr;
    loadingLbl->setHidden(false);
    //pImage->setHidden(true);
    pTransformation = pImage->Transform();

    if (pTransformation && pTMO)
    {
        pSrc = pImage->GetImage();
        if(pSrc){
            pTMO->SetImage(*pSrc);
            pTransformation->setPriority(QThread::LowestPriority);
            //pImage->setFixedSize(pSrc->GetWidth(), pSrc->GetHeight());
            pImage->fitToScreen(getImageViewSize());
            pTransformation->SetTMO(pTMO);
            bTransforming = true;

            return 0;
        }
        return -1;
    }

    return -1;
}

void TMOGUIPreviewImage::finishImageTransform(){
    bTransforming = false;
    bTransformed = true;
    loadingLbl->setHidden(true);
    //pImage->setHidden(false);
    pImage->fitToScreen(getImageViewSize());
    emit finishTransform();
}
