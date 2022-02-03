// TMOGUIProgressBar.h: interface for the TMOGUIProgressBar class.
//
//////////////////////////////////////////////////////////////////////
#ifndef TMOGUI_PROGRESSBAR_H
#define TMOGUI_PROGRESSBAR_H

#include <qmap.h>
//Added by qt3to4:
#include <QLabel>
#include <QWidget>
#include <QProgressBar>
#include <QPushButton>
#include <QHBoxLayout>

class QLabel;
class QStatusBar;
class TMOImage;

class TMOGUIProgressBar : public QWidget
{
    Q_OBJECT
public:
    int SetProgress(int part, int all);
    //virtual int Assign(TMOImage* pImage);
    virtual int SetLabel(QString text);
    virtual bool isVisible();
    TMOGUIProgressBar(QStatusBar *parent = 0, const char *name = 0) : QWidget((QWidget *)parent) //Q3HBox
    {
        QHBoxLayout *layout = new QHBoxLayout(this);
        setLayout(layout);
        iLast = -1;
        bVisible = false;
        bCancel = false;
        pParent = parent;
        // TODO layout - setSpacing(8);
        pProgress = new QProgressBar(this);        //, "ProgressBar");
        pButton = new QPushButton("Cancel", this); //, "CancelButton");
        pLabel = new QLabel("Working...", this);   //, "ProgressLabel");
        // TODO pProgress->setCenterIndicator(true);
        pProgress->setValue(0); //setProgress
        pProgress->setMinimumWidth(pProgress->width() * 1.5);
        pProgress->setFixedHeight(23);
        pButton->setFixedWidth(pButton->width() * .6);
        pButton->setFixedHeight(23);

        layout->addWidget(pProgress);
        layout->addWidget(pLabel);
        layout->addWidget(pButton);
        connect(pButton, SIGNAL(clicked()), this, SLOT(cancel()));
        hide();
    }
    virtual ~TMOGUIProgressBar();
    //static int ProgressBar(TMOImage* pImage, int part, int all);
protected:
    bool bCancel;
    bool bVisible;
    int iLast;
    QProgressBar *pProgress;
    QPushButton *pButton;
    QLabel *pLabel;
    QStatusBar *pParent;
    //static QMap<TMOImage*, TMOGUIProgressBar*> mapLocal;
public slots:
    void cancel();
signals:
    void cancelsignal();
};

#endif // TMOGUI_PROGRESSBAR_H
