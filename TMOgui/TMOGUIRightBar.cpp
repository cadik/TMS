// TMOGUIRightBar.cpp: implementation of the TMOGUIRightBar class.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"
#include "TMOGUIRightBar.h"
#include "TMOGUIToneMapping.h"
#include "TMOGUIStatistics.h"
#include "TMOGUIFilters.h"
#include "TMOGUIResource.h"
#include <qpushbutton.h>
#include <qtabwidget.h>
#include <qtabbar.h>
#include <qcursor.h>
#include <qlabel.h>
//Added by qt3to4:
#include <QResizeEvent>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIRightBar::TMOGUIRightBar(QWidget* parent, const char* name)
    : QTabWidget( parent)
{
	bVisible = false;
    bAdvanced = true;
    pRightButton = new QPushButton(this);//, "RightPushButton");
    pRightButton->setIcon(QIcon(":/resources/icons/IconX.png")); //setPixmap
	pRightButton->setFixedSize(16,16);
	pRightButton->setFlat(true);

    this->setCornerWidget(pRightButton);
    this->setTabShape(TabShape::Rounded);
    this->setMinimumWidth(420);
    setAutoFillBackground(true);


    //this->setFrameStyle(QFrame::Panel|QFrame::Raised);

    pStats = new TMOGUIStatistics(this/*, "Info"*/);
	pToneMapping = new TMOGUIToneMapping(this, "ToneMapping");
    pFiltersTab = new TMOGUIFilters(this, "Filters");

    addTab(pToneMapping, "TMO");
    addTab(pStats, "Info");
    addTab(pFiltersTab, "Filters");

    pTabBar = tabBar();
	connect (pRightButton, SIGNAL(clicked()), this, SLOT(hideright()));
    connect(pToneMapping, SIGNAL(change()), this, SLOT(paramChanged()));
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

void TMOGUIRightBar::changeWorkspace(bool advanced){
    if(advanced && !bAdvanced){
        pTabBar->setHidden(!advanced);
        addTab(pStats, "Info");
        addTab(pFiltersTab, "Filters");
    } else if (!advanced && bAdvanced){
        removeTab(indexOf(pStats));
        removeTab(indexOf(pFiltersTab));
        pTabBar->setHidden(!advanced);
    }
    bAdvanced = advanced;
    pToneMapping->changeWorkspace(bAdvanced);
}

void TMOGUIRightBar::hideright()
{
	emit closeBar();
}

void TMOGUIRightBar::paramChanged(){
    emit change();
}
