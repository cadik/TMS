// TMOGUIStatus.cpp: implementation of the TMOGUIStatus class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIStatus.h"
#include <qhbox.h>
#include <qlabel.h>
#include <qapplication.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIStatus::TMOGUIStatus( QWidget * parent, const char * name)  
	: QStatusBar(parent, name)
{
//	pProgress = 0;
	pPosition = 0;
	Create();
}

TMOGUIStatus::~TMOGUIStatus()
{
//	if (pProgress) delete pProgress;
//	pProgress = 0;
}

int TMOGUIStatus::Create()
{
//	message("Ready",2000);
//	pProgress = new TMOGUIProgressBar(this, "ProgressBar");
	pPosition = new QLabel(this, "Position");
	pPosition->setFrameStyle( QFrame::Panel);
	pPosition->setText("Tone mapping studio, 2004-2015");
	addWidget(pPosition,0,true);
	return 0;
}
