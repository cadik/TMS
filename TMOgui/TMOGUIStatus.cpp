// TMOGUIStatus.cpp: implementation of the TMOGUIStatus class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIStatus.h"
#include <qlabel.h>
#include <qapplication.h>
//Added by qt3to4:
#include <QFrame>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIStatus::TMOGUIStatus(QWidget *parent, const char *name)
    : QStatusBar(parent)
{
    pProgress = 0;
    pPosition = 0;
    Create();
}

TMOGUIStatus::~TMOGUIStatus()
{
    if (pProgress)
        delete pProgress;
    pProgress = 0;
}

int TMOGUIStatus::Create()
{
    //	message("Ready",2000); TODO Check status ready message
    this->setMaximumHeight(40);
    pProgress = new TMOGUIProgressBar(this, "ProgressBar");
    addWidget(pProgress, 1);
    pPosition = new QLabel(this); //, "Position");
    pPosition->setFrameStyle(QFrame::Panel);
    pPosition->setText("Tone mapping studio, 2004-2024");
    pPosition->setMaximumHeight(35);
    addPermanentWidget(pPosition, 0);
    return 0;
}
