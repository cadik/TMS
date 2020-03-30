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
#include <qcursor.h>
#include <qlabel.h>
//Added by qt3to4:
#include <QResizeEvent>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIRightBar::TMOGUIRightBar( QWidget* parent, const char* name, Qt::WindowFlags f )
    : QTabWidget( parent)
{
	bVisible = false;
    pRightButton = new QPushButton(this);//, "RightPushButton");
    pRightButton->setIcon(QIcon(":/resources/icons/IconX.png")); //setPixmap
	pRightButton->setFixedSize(16,16);
	pRightButton->setFlat(true);

    this->setCornerWidget(pRightButton);
    this->setTabShape(TabShape::Rounded);


    //this->setFrameStyle(QFrame::Panel|QFrame::Raised);

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
    if (!pToneMapping) return nullptr;

    if (!pToneMapping->pTMO) return nullptr;

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
