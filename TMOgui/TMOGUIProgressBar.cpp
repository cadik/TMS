// TMOGUIProgressBar.cpp: implementation of the TMOGUIProgressBar class.
//
//////////////////////////////////////////////////////////////////////
#include <q3progressbar.h>
#include <qstatusbar.h>
#include <qpushbutton.h>
#include <qapplication.h>
#include <qlabel.h>
#include "../tmolib/TMO.h"
#include "TMOGUIProgressBar.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
//QMap<TMOImage*, TMOGUIProgressBar*> TMOGUIProgressBar::mapLocal;

TMOGUIProgressBar::TMOGUIProgressBar(QStatusBar * parent, const char * name) :
	Q3HBox(parent, name)
{
	iLast = -1;
	bVisible = false;
	bCancel = false;
	pParent = parent;
	setSpacing(8);
	pProgress = new Q3ProgressBar(this, "ProgressBar");
	pButton = new QPushButton("Cancel", this, "CancelButton");
	pLabel = new QLabel("Working...", this, "ProgressLabel");
	pProgress->setCenterIndicator(true);
	pProgress->setProgress(0);
	pProgress->setMinimumWidth(pProgress->width()*1.5);
	pProgress->setFixedHeight(16);
	pButton->setFixedWidth(pButton->width()*.5);
	pButton->setFixedHeight(16);
	connect(pButton, SIGNAL(clicked()), this, SLOT(cancel()));
	hide();
}

TMOGUIProgressBar::~TMOGUIProgressBar()
{
}

int TMOGUIProgressBar::SetProgress(int part, int all)
{
	int iValue;

	if (all) iValue = (part * 100) / all;
	else iValue = part;

	if (!isVisible()) 
	{
		bCancel = false;
		show();
		pParent->addWidget(this);
	}
	bVisible = true;

	if (iValue == iLast) return bCancel;
	else iLast = iValue;

	pProgress->setProgress(iValue);
	qApp->processEvents();
		
	if (part == all) 
	{
		bVisible = false;
		hide();
		pParent->removeWidget(this);
	}
	
	return bCancel;
}

bool TMOGUIProgressBar::isVisible()
{
	return bVisible;
}

int TMOGUIProgressBar::SetLabel(QString text)
{
	pLabel->setText(text);
	return 0;
}

void TMOGUIProgressBar::cancel()
{
	bCancel = true;
	emit cancelsignal();
}

