#ifndef TMOGUIPREVIEWIMAGE_H
#define TMOGUIPREVIEWIMAGE_H

#include <QWidget>
#include <QFrame>
#include "TMOGUIImage.h"
//#include "TMOImage.h"
//#include "TMOGUITransformation.h"
#include <QPushButton>

class TMO;
class TMOImage;
class TMOGUITransformation;

class TMOGUIPreviewImage : public QFrame
{
    Q_OBJECT
public:
    TMOGUIPreviewImage(QWidget *parent, QString name, QPushButton *btn, TMOGUIImage* image, TMO* TMO);
    ~TMOGUIPreviewImage();
    int display();
    TMOGUIImage* getImage();
    void setImage(TMOGUIImage* srcImage);
    QSize getImageViewSize();
    QString* TMOname;
    bool bTransforming;
    bool bTransformed = false;
    QPushButton* pBtn;

protected:
    TMOGUIImage* pImage;
    TMOGUIImage* pDefaultImage;
    TMO* pTMO;
    QLabel* loadingLbl;
    TMOGUITransformation* pTransformation;
    QString* imageName;
    TMOGUIProgressBar* pBar;

private:
    QSize imageFitSize;

public slots:
    void finishImageTransform();

signals:
        void finishTransform();

};

#endif // TMOGUIPREVIEWIMAGE_H
