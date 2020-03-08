// TMOGUIAdjust.cpp: implementation of the TMOGUIAdjust class.
//
//////////////////////////////////////////////////////////////////////
//#include "../tmolib/TMO.h"
#include "TMOGUIAdjust.h"
#include "TMOGUIHisto.h"
#include "TMOGUIProgressBar.h"
#include "TMOGUIToneSlider.h"
#include "TMOGUIResource.h"
#include <math.h>
//#include <qslider.h>
#include <QSlider>
//#include <qlayout.h>
#include <QLayout>
//#include <qpushbutton.h>
#include <QPushButton>
//#include <qlineedit.h>
#include <QLineEdit>
//#include <qlabel.h>
#include <QLabel>
//#include <qtooltip.h>
#include <QToolTip>
//Added by qt3to4:
#include <QHBoxLayout>
#include <QGridLayout>
#include <QFrame>
#include "gamma.h"

//#include <iostream>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIAdjust::TMOGUIAdjust(QWidget* parent, const char * name):
    QFrame(parent)
{
	iMode = 0;
    pValues = nullptr;
	bLog = true;
	bUpdate = false;

    setFrameStyle( QFrame::Panel | QFrame::Sunken );
    QGridLayout *pLayout = new QGridLayout(this); //, 8, 4);
    pLayout->setContentsMargins(0,0,0,0);
    pLayout->setSpacing(0);
    pLayout->setColumnMinimumWidth(0,30);
    pLayout->setColumnStretch(0,0);
    pLayout->setColumnStretch(1,1);
    pLayout->setColumnStretch(2,0);
    pLayout->setColumnStretch(3,0);
    pLayout->addItem(new QSpacerItem(20,0), 0, 2); //pLayout->addColSpacing(2,20);

	
    QSlider *pSlider = new QSlider(this); //, "HistoSlider");
    pSlider->setFixedSize(30,64);
	pSlider->setRange(-50, 50);
	pSlider->setValue(0);
	pSlider->setOrientation(Qt::Vertical);
    // TODO tickmarks
    pSlider->setTickPosition(QSlider::TicksRight);
    pSlider->setTickInterval(10);
    pLayout->addWidget(pSlider, 0, 0, 4, 0, Qt::AlignLeft); //addMultiCellWidget(pSlider, 0, 3, 0, 0);

    pHisto = new TMOGUIHisto(this, "Histogram");

    pLayout->addWidget(pHisto, 0, 1,  4,  1, Qt::AlignBottom); //addMultiCellWidget(pHisto, 0, 3, 1, 1, Qt::AlignBottom);

    pLayout->addItem(new QSpacerItem(0,2), 4, 0, 4, 0);//addRowSpacing(4, 2);

	pToneSlider = new TMOGUIToneSlider(this, "ToneSlider");
    pLayout->addWidget(pToneSlider, 5, 1);

    pRed = new QPushButton(this); // Red Button
	pRed->setFlat(true);
    pRed->setFixedSize(16,16);
    pRed->setIcon(QIcon(QPixmap(":/resources/icons/IconRed.png")));// setIcon(QIcon(":/resources/icons/IconRed.png"));
    pRed->setToolTip("Select Red Channel");
    pLayout->addWidget(pRed, 0, 2, Qt::AlignCenter);
	
    pGreen = new QPushButton(this); // Green Button
	pGreen->setFlat(true);
    pGreen->setFixedSize(16,16);
    pGreen->setIcon(QIcon(":/resources/icons/IconGreen.png"));//setIcon(QIcon(":/resources/icons/IconGreen.png"));
    pGreen->setToolTip("Select Green Channel");
    pLayout->addWidget(pGreen, 1, 2, Qt::AlignCenter);

    pBlue = new QPushButton(this); // Blue Button
	pBlue->setFlat(true);
    pBlue->setFixedSize(16,16);
    pBlue->setIcon(QIcon(":/resources/icons/IconBlue.png"));
    pBlue->setToolTip("Select Blue Channel");
    pLayout->addWidget(pBlue, 2, 2, Qt::AlignCenter);

    pLum = new QPushButton(this);
	pLum->setFlat(true);
    pLum->setFixedSize(16,16);
    pLum->setIcon(QIcon(":/resources/icons/IconLumDown.png"));
    pLum->setToolTip("Select Luminance");
    pLayout->addWidget(pLum, 3, 2, Qt::AlignCenter);
	
    pLinear = new QPushButton(this);
	pLinear->setFlat(true);
    pLinear->setFixedSize(32,32);
    pLinear->setIcon(QIcon(":/resources/icons/IconLin.png"));
    pLinear->setToolTip("View In Linear Distribution");
    //pLayout->addMultiCellWidget(pLinear, 0, 1, 3, 3, Qt::AlignCenter);
    pLayout->addWidget(pLinear, 0, 3, 2, 1, Qt::AlignCenter);

    pLog = new QPushButton(this);
	pLog->setFlat(true);
    pLog->setFixedSize(32,32);
    pLog->setIcon(QIcon(":/resources/icons/IconLogDown.png"));
    pLog->setToolTip("View In Logarithmic Distribution");
    //pLayout->addMultiCellWidget(pLog, 2, 3, 3, 3, Qt::AlignCenter);
    pLayout->addWidget(pLog, 2, 3, 2, 1, Qt::AlignCenter);
    pLayout->addItem(new QSpacerItem(0,2), 6, 0);//addRowSpacing(6, 2);

    QHBoxLayout* hbox = new QHBoxLayout();
    QLabel* pLabel = new QLabel(this);
    pLabel->setPixmap(QPixmap(":/resources/icons/IconBlack.png"));
    pLabel->setFrameStyle( QFrame::Panel | QFrame::Sunken );
	hbox->addWidget(pLabel);
    pBlack = new QLineEdit(this);
    pBlack->setAlignment(Qt::AlignRight);
    pBlack->setFixedWidth(72);
	hbox->addWidget(pBlack);
	hbox->insertStretch(-1,1);
    pLabel = new QLabel(this);
    pLabel->setPixmap(QPixmap(":/resources/icons/IconGamma.png"));
    pLabel->setFrameStyle( QFrame::Panel | QFrame::Sunken );
	hbox->addWidget(pLabel);
    pGamma = new QLineEdit(this);
    pGamma->setAlignment(Qt::AlignRight);
	pGamma->setFixedWidth(72);
	hbox->addWidget(pGamma);
	hbox->insertStretch(-1,1);
    pLabel = new QLabel(this);
    pLabel->setPixmap(QPixmap(":/resources/icons/IconWhite.png"));
    pLabel->setFrameStyle( QFrame::Panel | QFrame::Sunken );
	hbox->addWidget(pLabel);
    pWhite = new QLineEdit(this);
    pWhite->setAlignment(Qt::AlignRight);
	pWhite->setFixedWidth(72);
	hbox->addWidget(pWhite);
    pLayout->addLayout(hbox, 7, 1);

    setLayout(pLayout);


    connect (pSlider, &QSlider::sliderMoved, this, &TMOGUIAdjust::scalechanged);
	connect (pLinear, SIGNAL(clicked()), pHisto, SLOT(setlinear()));
	connect (pLog, SIGNAL(clicked()), pHisto, SLOT(setlog()));
	connect (pRed, SIGNAL(clicked()), pHisto, SLOT(setr()));
	connect (pGreen, SIGNAL(clicked()), pHisto, SLOT(setg()));
	connect (pBlue, SIGNAL(clicked()), pHisto, SLOT(setb()));
	connect (pLum, SIGNAL(clicked()), pHisto, SLOT(setl()));
	connect (pLinear, SIGNAL(clicked()), pToneSlider, SLOT(setlinear()));
	connect (pLog, SIGNAL(clicked()), pToneSlider, SLOT(setlog()));
	connect (pRed, SIGNAL(clicked()), pToneSlider, SLOT(setr()));
	connect (pGreen, SIGNAL(clicked()), pToneSlider, SLOT(setg()));
	connect (pBlue, SIGNAL(clicked()), pToneSlider, SLOT(setb()));
	connect (pLum, SIGNAL(clicked()), pToneSlider, SLOT(setl()));
	connect (pRed, SIGNAL(clicked()), this, SLOT(setr()));
	connect (pGreen, SIGNAL(clicked()), this, SLOT(setg()));
	connect (pBlue, SIGNAL(clicked()), this, SLOT(setb()));
	connect (pLum, SIGNAL(clicked()), this, SLOT(setl()));
    connect (pBlack, &QLineEdit::textChanged, this, &TMOGUIAdjust::setblack);
    connect (pWhite, &QLineEdit::textChanged, this, &TMOGUIAdjust::setwhite);
    connect (pGamma, &QLineEdit::textChanged, this, &TMOGUIAdjust::setgamma);
	connect (pBlack, SIGNAL(returnPressed()), this, SLOT(updateall()));
	connect (pWhite, SIGNAL(returnPressed()), this, SLOT(updateall()));
	connect (pGamma, SIGNAL(returnPressed()), this, SLOT(updateall()));
	connect (pLinear, SIGNAL(clicked()), this, SLOT(setlinear()));
	connect (pLog, SIGNAL(clicked()), this, SLOT(setlog()));
}

TMOGUIAdjust::~TMOGUIAdjust()
{

}

int TMOGUIAdjust::Create(TMOImage* pSrc, TMOGUIAdjustValues* pVals)
{
	pValues = pVals;
	pHisto->Create(pSrc, pVals);
	pToneSlider->Create(pVals);
	connect (pVals, SIGNAL(valueschanged()), pHisto, SLOT(update()));
	connect (pVals, SIGNAL(valueschanged()), this, SLOT(valueschanged()));
	valueschanged();
	return 0;
}

int TMOGUIAdjust::SetImage(TMOImage *pSrc)
{
	pHisto->SetWindow(pSrc);
	return 0;
}

void TMOGUIAdjust::scalechanged(int i)
{
	pHisto->setscale(i);
}


int TMOGUIAdjust::reset()
{
	pHisto->reset(0);
	pHisto->reset(1);
	pHisto->reset(2);
	pHisto->reset(3);
	return 0;
}

void TMOGUIAdjust::valueschanged()
{
	bUpdate = false;
	switch(iMode)
	{
	case 0:
		pBlack->setText(FormatVal(0.299 * pValues->dRMinimum + 0.587 * pValues->dGMinimum + 0.114 * pValues->dBMinimum));
		pWhite->setText(FormatVal(0.299 * pValues->dRMaximum + 0.587 * pValues->dGMaximum + 0.114 * pValues->dBMaximum));
		pGamma->setText(FormatVal(1/(0.299 * pValues->dRGamma + 0.587 * pValues->dGGamma + 0.114 * pValues->dBGamma))); //pAn
		break;
	case 1:
		pBlack->setText(FormatVal(pValues->dRMinimum));
		pWhite->setText(FormatVal(pValues->dRMaximum));
		pGamma->setText(FormatVal(1/pValues->dRGamma));
		break;
	case 2:
		pBlack->setText(FormatVal(pValues->dGMinimum));
		pWhite->setText(FormatVal(pValues->dGMaximum));
		pGamma->setText(FormatVal(1/pValues->dGGamma));
		break;
	case 3:
		pBlack->setText(FormatVal(pValues->dBMinimum));
		pWhite->setText(FormatVal(pValues->dBMaximum));
		pGamma->setText(FormatVal(1/pValues->dBGamma));
		break;
	}
	bUpdate = true;
}

void TMOGUIAdjust::setr()
{
    pRed->setIcon(QIcon(":/resources/icons/IconRedDown.png"));
    pGreen->setIcon(QIcon(":/resources/icons/IconGreen.png"));
    pBlue->setIcon(QIcon(":/resources/icons/IconBlue.png"));
    pLum->setIcon(QIcon(":/resources/icons/IconLum.png"));
	iMode = 1;
	valueschanged();
}

void TMOGUIAdjust::setg()
{
    pRed->setIcon(QIcon(":/resources/icons/IconRed.png"));
    pGreen->setIcon(QIcon(":/resources/icons/IconGreenDown.png"));
    pBlue->setIcon(QIcon(":/resources/icons/IconBlue.png"));
    pLum->setIcon(QIcon(":/resources/icons/IconLum.png"));
	iMode = 2;
	valueschanged();
}

void TMOGUIAdjust::setb()
{
    pRed->setIcon(QIcon(":/resources/icons/IconRed.png"));
    pGreen->setIcon(QIcon(":/resources/icons/IconGreen.png"));
    pBlue->setIcon(QIcon(":/resources/icons/IconBlueDown.png"));
    pLum->setIcon(QIcon(":/resources/icons/IconLum.png"));
	iMode = 3;
	valueschanged();
}

void TMOGUIAdjust::setl()
{
    pRed->setIcon(QIcon(":/resources/icons/IconRed.png"));
    pGreen->setIcon(QIcon(":/resources/icons/IconGreen.png"));
    pBlue->setIcon(QIcon(":/resources/icons/IconBlue.png"));
    pLum->setIcon(QIcon(":/resources/icons/IconLumDown.png"));
	iMode = 0;
	valueschanged();
}

void TMOGUIAdjust::setblack(const QString& s)
{
	double value = s.toDouble(), min, max;

	if (!bUpdate) return;
	if (bLog) min = exp(pValues->dMinimum);
	else min = pValues->dMinimum;

	switch(iMode)
	{
	case 0:
		max = 0.299 * pValues->dRMaximum + 0.587 * pValues->dGMaximum + 0.114 * pValues->dBMaximum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dRMinimum = pValues->dGMinimum = pValues->dBMinimum = value;
		pToneSlider->setl();
		pHisto->update();
		break;
	case 1:
		max = pValues->dRMaximum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dRMinimum = value;
		pToneSlider->setr();
		pHisto->update();
		break;
	case 2:
		max = pValues->dGMaximum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dGMinimum = value;
		pToneSlider->setg();
		pHisto->update();
		break;
	case 3:
		max = pValues->dBMaximum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dBMinimum = value;
		pToneSlider->setb();
		pHisto->update();
		break;
	}
}


void TMOGUIAdjust::setgamma(const QString& s)
{
	double value = s.toDouble(), min, max;//pAn
	min = GAMMA_MIN;
	max = GAMMA_MAX;
	if (!bUpdate) return;
		
	if(value < GAMMA_MIN)
	{
	 value = GAMMA_MIN;
	}
	if(value > GAMMA_MAX)
	{
	 value = GAMMA_MAX;
	}

	switch(iMode)
	{
	case 0:
		pValues->dRGamma = pValues->dGGamma = pValues->dBGamma = 1/value;
		pToneSlider->setl();
		break;
	case 1:
		pValues->dRGamma = 1/value;
		pToneSlider->setr();
		pHisto->update();
		break;
	case 2:
		pValues->dGGamma = 1/value;
		pToneSlider->setg();
		pHisto->update();
		break;
	case 3:
		pValues->dBGamma = 1/value;
		pToneSlider->setb();
		pHisto->update();
		break;
	}
}

void TMOGUIAdjust::setwhite(const QString& s)
{
	double value = s.toDouble(), min, max;

	if (!bUpdate) return;
	if (bLog) max = exp(pValues->dMaximum);
	else max = pValues->dMaximum;

	switch(iMode)
	{
	case 0:
		min = 0.299 * pValues->dRMinimum + 0.587 * pValues->dGMinimum + 0.114 * pValues->dBMinimum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dRMaximum = pValues->dGMaximum = pValues->dBMaximum = value;
		pToneSlider->setl();
		pHisto->update();
		break;
	case 1:
		min = pValues->dRMinimum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dRMaximum = value;
		pToneSlider->setr();
		pHisto->update();
		break;
	case 2:
		min = pValues->dGMinimum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dGMaximum = value;
		pToneSlider->setg();
		pHisto->update();
		break;
	case 3:
		min = pValues->dBMinimum;
		value = value > min ? value : min;
		value = value < max ? value : max;
		pValues->dBMaximum = value;
		pToneSlider->setb();
		pHisto->update();
		break;
	}
}

void TMOGUIAdjust::setlog()
{
    pLog->setIcon(QIcon(":/resources/icons/IconLogDown.png"));
    pLinear->setIcon(QIcon(":/resources/icons/IconLin.png"));
	bLog = true;
}

void TMOGUIAdjust::setlinear()
{
    pLog->setIcon(QIcon(":/resources/icons/IconLog.png"));
    pLinear->setIcon(QIcon(":/resources/icons/IconLinDown.png"));
	bLog = false;
}

void TMOGUIAdjust::updateall()
{
	double gamma = pGamma->text().toDouble();
	if(gamma < GAMMA_MIN)
	{
	 gamma = GAMMA_MIN;
         pGamma->setText(FormatVal(gamma));
	}
	if(gamma > GAMMA_MAX)
	{
	 gamma = GAMMA_MAX;
         pGamma->setText(FormatVal(gamma));
	}

	pValues->RenderValues();
}

QString TMOGUIAdjust::FormatVal(double val)
{
	QString s;
	if ((val < 0.0001)||(val > 99999.999))
		s.setNum(val,'g',3);
	else if (val <= 0.999999)
		s.setNum(val,'f',6);
	else s.setNum(val,'f',3);
	return s;
}


