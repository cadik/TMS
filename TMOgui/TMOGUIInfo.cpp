// TMOGUIInfo.cpp: implementation of the TMOGUIInfo class.
//
//////////////////////////////////////////////////////////////////////
#include "TMOGUIInfo.h"
#include "TMOGUIOutput.h"
#include "TMOGUIStatistics.h"
#include "TMOGUIResource.h"
#include <qpushbutton.h>
#include <qlabel.h>
//Added by qt3to4:
#include <QResizeEvent>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIInfo::TMOGUIInfo(QWidget* parent, const char * name) :
	QTabWidget(parent, name)
{
	pOutput = 0;
	pLabel = new QLabel("   No Window Selected", this);
		
	pStats = new TMOGUIStatistics(this, "Statistics");
	addTab(pStats, "&Info");

	pClose = new QPushButton(this, "CloseButton");
	pClose->setPixmap(*TMOResource::pResource->IconX->pixmap());
	pClose->setFixedSize(16,16);
	pClose->setFlat(true);

	SetOutput(0);
	
	connect (pClose, SIGNAL(clicked()), this, SLOT(hideinfo()));
}

TMOGUIInfo::~TMOGUIInfo()
{

}

void TMOGUIInfo::resizeEvent ( QResizeEvent * re)
{
	pClose->move(re->size().width()-16, 0);
	QTabWidget::resizeEvent(re);
}

int TMOGUIInfo::SetOutput(QWidget* pOut)
{
	int iCurrentID = currentPageIndex();
	if (pOutput) removePage(pOutput);
	if (pOut) pOutput = pOut;
	else pOutput = pLabel;
	addTab(pOutput, "&Output");
	if (currentPageIndex() != iCurrentID)	showPage(pOutput); 
	return 0;
}

void TMOGUIInfo::hideinfo()
{
	emit closeBar();
}