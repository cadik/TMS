// TMOGUIProgressBar.cpp: implementation of the TMOGUIProgressBar class.
//
//////////////////////////////////////////////////////////////////////
#include <QProgressBar>
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

TMOGUIProgressBar::~TMOGUIProgressBar()
{
}

int TMOGUIProgressBar::SetProgress(int part, int all)
{
	int iValue;

	if (all)
		iValue = (part * 100) / all;
	else
		iValue = part;

	if (!isVisible())
	{
		bCancel = false;
		show();
		pParent->addWidget(this);
	}
	bVisible = true;

	if (iValue == iLast)
		return bCancel;
	else
		iLast = iValue;

	pProgress->setValue(iValue); //setProgress
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
