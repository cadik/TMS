// TMOGUIRightBar.cpp: implementation of the TMOGUIRightBar class.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"
#include "TMOGUIRightBar.h"
#include "TMOGUIToneMapping.h"
#include "TMOGUIFilters.h"
#include "TMOGUIResource.h"
#include <qpushbutton.h>
#include <qtabwidget.h>
#include <q3vbox.h>
#include <qcursor.h>
#include <qlabel.h>
//Added by qt3to4:
#include <QResizeEvent>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIRightBar::TMOGUIRightBar( QWidget* parent, const char* name, Qt::WFlags f )
	: QTabWidget( parent, name, f )
{
	bVisible = false;
	pRightButton = new QPushButton("X", this, "RightPushButton");
	pRightButton->setPixmap(*TMOResource::pResource->IconX->pixmap());
	pRightButton->setFixedSize(16,16);
	pRightButton->setFlat(true);

	//setFrameStyle(Panel|Raised);

	pToneMapping = new TMOGUIToneMapping(this, "ToneMapping");
	pFilters = new TMOGUIFilters(this, "Filters");
	addTab(pToneMapping, "TMO");
	addTab(pFilters, "Filters");

	connect (pRightButton, SIGNAL(clicked()), this, SLOT(hideright()));
}

TMOGUIRightBar::~TMOGUIRightBar()
{

}

TMO* TMOGUIRightBar::GetTMO()
{
	if (!pToneMapping) return 0;

	if (!pToneMapping->pTMO) return 0;

	return pToneMapping->pTMO[pToneMapping->iCurTechnique];
}

void TMOGUIRightBar::resizeEvent ( QResizeEvent * re)
{
	pRightButton->move(re->size().width()-16, 0);
	QTabWidget::resizeEvent(re);
}

void TMOGUIRightBar::hideright()
{
	emit closeBar();
}