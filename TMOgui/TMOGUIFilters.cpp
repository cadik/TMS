// TMOGUIFilters.cpp: implementation of the TMOGUIFilters class.
//
//////////////////////////////////////////////////////////////////////
#include <math.h>
#include "TMOGUIAdjust.h"
#include "TMOGUIFilters.h"
#include "TMOGUIImage.h"
#include <qlabel.h>
#include <qlineedit.h>
#include <qscrollbar.h>
#include <qlayout.h>
#include <qradiobutton.h>
#include <qpushbutton.h>
#include <QButtonGroup>
#include <qcheckbox.h>
#include <QGroupBox>
#include <qpixmap.h>
#include <qslider.h>
#include <qcolordialog.h>
//Added by qt3to4:
#include <QGridLayout>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIFilters::TMOGUIFilters(QWidget *parent, const char *name) : QWidget(parent)
{
	setAutoFillBackground(true);

	pWhite = pRWhite = pGWhite = pBWhite = 0;
	pBlack = pRBlack = pGBlack = pBBlack = 0;
	pLuminance = pRed = pGreen = pBlue = 0;
	pPixmap = 0;
	pImage = 0;
	setObjectName(name);
	Create();
}

TMOGUIFilters::~TMOGUIFilters()
{
	if (pPixmap)
		delete pPixmap;
}

int TMOGUIFilters::Create()
{
	QRadioButton *pRadio;
	QButtonGroup *pButtonGroup;

	pLayout = new QGridLayout(this); //,24,7);

	pLayout->addItem(new QSpacerItem(10, 0), 0, 0);	 //addColSpacing(0,10);
	pLayout->addItem(new QSpacerItem(10, 0), 0, 6);	 //addColSpacing(6,10);
	pLayout->addItem(new QSpacerItem(10, 0), 0, 3);	 //addColSpacing(3,10);
	pLayout->addItem(new QSpacerItem(0, 10), 0, 0);	 //addRowSpacing(0,10);
	pLayout->addItem(new QSpacerItem(0, 20), 1, 0);	 //addRowSpacing(1,20);
	pLayout->addItem(new QSpacerItem(0, 10), 4, 0);	 //addRowSpacing(4,10);
	pLayout->addItem(new QSpacerItem(0, 10), 5, 0);	 //addRowSpacing(5,10);
	pLayout->addItem(new QSpacerItem(0, 10), 7, 0);	 //addRowSpacing(7,10);
	pLayout->addItem(new QSpacerItem(0, 10), 14, 0); //addRowSpacing(14,10);
	pLayout->addItem(new QSpacerItem(0, 10), 18, 0); //addRowSpacing(18,10);
	pLayout->addItem(new QSpacerItem(0, 10), 21, 0); //addRowSpacing(21,10);
	pLayout->addItem(new QSpacerItem(0, 10), 23, 0); //addRowSpacing(23,10);
	pLayout->setColumnStretch(1, 1);
	pLayout->setColumnStretch(4, 1);
	//pLayout->setColumnStretch(0,0); //?
	pLayout->setColumnStretch(6, 0);
	pLayout->setRowStretch(17, 1);

	pButtonGroup = new QButtonGroup(this); // "Threshold", this, "GroupButton");

	//pLayout->addWidget(pButtonGroup, 1, 1, 4, 5, 0);
	//addMultiCellWidget(pButtonGroup,1,4,1,5);

	pRadio = new QRadioButton("Luminance", this); //, "LuminanceButton" );
	pButtonGroup->addButton(pRadio);
	connect(pRadio, SIGNAL(clicked()), this, SLOT(luminance()));
	pLayout->addWidget(pRadio, 2, 2, 1, 3);
	//addMultiCellWidget(pRadio,2,2,2,4);
	pRadio = new QRadioButton("Component", this); //, "ComponentButton" );
	pButtonGroup->addButton(pRadio);
	connect(pRadio, SIGNAL(clicked()), this, SLOT(component()));
	pRadio->setChecked(true);
	pLayout->addWidget(pRadio, 3, 2, 1, 3);
	//addMultiCellWidget(pRadio,3,3,2,4);

	pLuminance = new QLabel("Luminannce", this);
	pLayout->addWidget(pLuminance, 6, 1, 1, 6);
	//addMultiCellWidget(pLuminance,6,6,1,6);
	pBlackScroll = new QScrollBar(Qt::Horizontal, this); //, "BlackScrollBar");
	pBlackScroll->setRange(0, 100);
	pLayout->addWidget(pBlackScroll, 7, 1);
	pBlack = new QLineEdit(this); //, "BlackEdit" );
	pBlack->setFixedWidth(48);
	pLayout->addWidget(pBlack, 7, 2);
	pWhiteScroll = new QScrollBar(Qt::Horizontal, this); //, "WhiteScrollBar");
	pWhiteScroll->setRange(0, 100);
	pLayout->addWidget(pWhiteScroll, 7, 4);
	pWhite = new QLineEdit(this); //, "WhiteEdit" );
	pWhite->setFixedWidth(48);
	pLayout->addWidget(pWhite, 7, 5);

	pRed = new QLabel("Red", this);
	pLayout->addWidget(pRed, 8, 1, 1, 6);
	//addMultiCellWidget(pRed,8,8,1,6);
	pRBlackScroll = new QScrollBar(Qt::Horizontal, this); //, "RBlackScrollBar");
	pRBlackScroll->setRange(0, 100);
	pLayout->addWidget(pRBlackScroll, 9, 1);
	pRBlack = new QLineEdit(this); //, "RBlackEdit" );
	pRBlack->setFixedWidth(48);
	pLayout->addWidget(pRBlack, 9, 2);
	pRWhiteScroll = new QScrollBar(Qt::Horizontal, this); //, "RWhiteScrollBar");
	pRWhiteScroll->setRange(0, 100);
	pLayout->addWidget(pRWhiteScroll, 9, 4);
	pRWhite = new QLineEdit(this); //, "RWhiteEdit" );
	pRWhite->setFixedWidth(48);
	pLayout->addWidget(pRWhite, 9, 5);

	pGreen = new QLabel("Green", this);
	pLayout->addWidget(pGreen, 10, 1, 1, 6);
	//addMultiCellWidget(pGreen,10,10,1,6);
	pGBlackScroll = new QScrollBar(Qt::Horizontal, this); //, "RBlackScrollBar");
	pGBlackScroll->setRange(0, 100);
	pLayout->addWidget(pGBlackScroll, 11, 1);
	pGBlack = new QLineEdit(this); //, "RBlackEdit" );
	pGBlack->setFixedWidth(48);
	pLayout->addWidget(pGBlack, 11, 2);
	pGWhiteScroll = new QScrollBar(Qt::Horizontal, this); //, "RWhiteScrollBar");
	pGWhiteScroll->setRange(0, 100);
	pLayout->addWidget(pGWhiteScroll, 11, 4);
	pGWhite = new QLineEdit(this); //, "RWhiteEdit" );
	pGWhite->setFixedWidth(48);
	pLayout->addWidget(pGWhite, 11, 5);

	pBlue = new QLabel("Blue", this);
	pLayout->addWidget(pBlue, 12, 1, 1, 6);
	//addMultiCellWidget(pBlue,12,12,1,6);
	pBBlackScroll = new QScrollBar(Qt::Horizontal, this); //, "RBlackScrollBar");
	pBBlackScroll->setRange(0, 100);
	pLayout->addWidget(pBBlackScroll, 13, 1);
	pBBlack = new QLineEdit(this); //, "RBlackEdit" );
	pBBlack->setFixedWidth(48);
	pLayout->addWidget(pBBlack, 13, 2);
	pBWhiteScroll = new QScrollBar(Qt::Horizontal, this); //, "RWhiteScrollBar");
	pBWhiteScroll->setRange(0, 100);
	pLayout->addWidget(pBWhiteScroll, 13, 4);
	pBWhite = new QLineEdit(this); //, "RWhiteEdit" );
	pBWhite->setFixedWidth(48);
	pLayout->addWidget(pBWhite, 13, 5);

	pOutside = new QCheckBox("Outside values", this); //, "OutsideCheck");
	pLayout->addWidget(pOutside, 15, 1, 1, 5);
	//addMultiCellWidget(pOutside, 15, 15, 1, 5);
	connect(pOutside, &QCheckBox::toggled, this, &TMOGUIFilters::outsidevalues);
	pConnect = new QCheckBox("Connect values", this); //, "ConnectCheck");
	pLayout->addWidget(pConnect, 16, 1, 1, 5);
	//addMultiCellWidget(pConnect, 16, 16, 1, 5);
	connect(pConnect, &QCheckBox::toggled, this, &TMOGUIFilters::connectvalues);

	QLabel *pLabel = new QLabel("Intensity", this); //, "FilterInt");
	pLayout->addWidget(pLabel, 19, 1, 1, 2);
	//addMultiCellWidget(pLabel, 19, 19, 1, 2);
	pLabel = new QLabel("Color", this); //, "ColorLabel");
	pLayout->addWidget(pLabel, 19, 5, Qt::AlignHCenter);

	pPixmap = new QPixmap(32, 24);
	pPixmap->fill(QColor(255, 0, 0));

	pSlider = new QSlider(Qt::Horizontal, this); //, "Color");
	pSlider->setRange(0, 100);
	pLayout->addWidget(pSlider, 20, 1, 1, 4);
	//addMultiCellWidget(pSlider, 20, 20, 1, 4);

	pColor = new QPushButton(this); //, "ColorBox");
	pColor->setIcon(QIcon(*pPixmap));
	// TODO? pColor->setIconSize();
	pColor->setFixedSize(32, 24);
	pLayout->addWidget(pColor, 20, 5);

	pEnable = new QPushButton("Enable", this); //, "EnableFilters");
	pEnable->setCheckable(true);
	pEnable->setFixedHeight(24);
	pEnable->setFixedWidth(64);
	pLayout->addWidget(pEnable, 22, 2, Qt::AlignRight);

	QPushButton *pReset = new QPushButton("Reset", this); //, "ResetFilters");
	pReset->setFixedHeight(24);
	pReset->setFixedWidth(64);
	pLayout->addWidget(pReset, 22, 4, 1, 5, Qt::AlignLeft);
	//addMultiCellWidget(pReset, 22, 22, 4, 5, Qt::AlignLeft);

	component();

	connect(pBlackScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::blackslider);
	connect(pBlackScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pBlack, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pWhiteScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::whiteslider);
	connect(pWhiteScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pWhite, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pRBlackScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::rblackslider);
	connect(pRBlackScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pRBlack, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pRWhiteScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::rwhiteslider);
	connect(pRWhiteScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pRWhite, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pGBlackScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::gblackslider);
	connect(pGBlackScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pGBlack, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pGWhiteScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::gwhiteslider);
	connect(pGWhiteScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pGWhite, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pBBlackScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::bblackslider);
	connect(pBBlackScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pBBlack, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pBWhiteScroll, &QScrollBar::valueChanged, this, &TMOGUIFilters::bwhiteslider);
	connect(pBWhiteScroll, SIGNAL(sliderReleased()), this, SLOT(color()));
	connect(pBWhite, SIGNAL(returnPressed()), this, SLOT(updatesliders()));
	connect(pEnable, &QPushButton::toggled, this, &TMOGUIFilters::enablefilters);
	connect(pReset, SIGNAL(clicked()), this, SLOT(resetvalues()));
	connect(pSlider, SIGNAL(sliderReleased()), this, SLOT(intensity()));

	connect(pColor, SIGNAL(clicked()), this, SLOT(colorchange()));

	return 0;
}

void TMOGUIFilters::windowChanged(TMOGUIImage *pWindow)
{
	TMOGUIImage *pOldImage = pImage;
	pImage = pWindow;
	if (pImage && pImage->pImage)
	{
		if (pImage->pImage->pFilters)
		{
			pImage->pImage->pFilters->dBAverage = 1.0;	// Component / Luminance
			pImage->pImage->pFilters->dRGamma = 1.0;	// Color
			pImage->pImage->pFilters->dGGamma = 0.0;	// .
			pImage->pImage->pFilters->dBGamma = 0.0;	// .
			pImage->pImage->pFilters->dMinimum = 0.0;	// Luminance
			pImage->pImage->pFilters->dMaximum = 0.0;	// .
			pImage->pImage->pFilters->dAverage = -1.0;	// enabled
			pImage->pImage->pFilters->dRAverage = -1.0; // outside
			pImage->pImage->pFilters->dGAverage = -1.0; // connected
			pImage->pImage->pFilters->dExtreme = 0.5;	// intensity
		}

		connect(pImage->pAdjust->pValues, SIGNAL(valueschanged()), this, SLOT(valueschanged()));
		valueschanged();
		refresh();
		color();
	}
}

void TMOGUIFilters::valueschanged()
{
	if (pImage)
	{
		if (pImage->pImage->pFilters->dGAverage > 0.0)
		{
			if (pImage->pImage->pFilters->dBAverage > 0.0)
			{
				pImage->pImage->pFilters->dRMinimum = pImage->pAdjust->pValues->dRMinimum;
				pImage->pImage->pFilters->dGMinimum = pImage->pAdjust->pValues->dGMinimum;
				pImage->pImage->pFilters->dBMinimum = pImage->pAdjust->pValues->dBMinimum;

				pImage->pImage->pFilters->dRMaximum = pImage->pAdjust->pValues->dRMaximum;
				pImage->pImage->pFilters->dGMaximum = pImage->pAdjust->pValues->dGMaximum;
				pImage->pImage->pFilters->dBMaximum = pImage->pAdjust->pValues->dBMaximum;
				refresh();
			}
			else
			{
				pImage->pImage->pFilters->dMinimum = pImage->pAdjust->pValues->dRMinimum;
				pImage->pImage->pFilters->dMaximum = pImage->pAdjust->pValues->dRMaximum;
				refresh();
			}
		}
	}
}

void TMOGUIFilters::refresh()
{
	QString sNum;

	if (pImage && pBlack)
	{
		sNum.setNum(pImage->pImage->pFilters->dMinimum, 'f', 2);
		pBlack->setText(sNum);
		pBlackScroll->setValue(pImage->pImage->pFilters->dMinimum / GetMaxValue());

		sNum.setNum(pImage->pImage->pFilters->dRMinimum, 'f', 2);
		pRBlack->setText(sNum);
		pRBlackScroll->setValue(100 * pImage->pImage->pFilters->dRMinimum / GetMaxValue());

		sNum.setNum(pImage->pImage->pFilters->dGMinimum, 'f', 2);
		pGBlack->setText(sNum);
		pGBlackScroll->setValue(100 * pImage->pImage->pFilters->dGMinimum / GetMaxValue());

		sNum.setNum(pImage->pImage->pFilters->dBMinimum, 'f', 2);
		pBBlack->setText(sNum);
		pBBlackScroll->setValue(100 * pImage->pImage->pFilters->dBMinimum / GetMaxValue());

		sNum.setNum(pImage->pImage->pFilters->dMaximum, 'f', 2);
		pWhite->setText(sNum);
		pWhiteScroll->setValue(100 * pImage->pImage->pFilters->dMaximum / GetMaxValue());

		sNum.setNum(pImage->pImage->pFilters->dRMaximum, 'f', 2);
		pRWhite->setText(sNum);
		pRWhiteScroll->setValue(100 * pImage->pImage->pFilters->dRMaximum / GetMaxValue());

		sNum.setNum(pImage->pImage->pFilters->dGMaximum, 'f', 2);
		pGWhite->setText(sNum);
		pGWhiteScroll->setValue(100 * pImage->pImage->pFilters->dGMaximum / GetMaxValue());

		sNum.setNum(pImage->pImage->pFilters->dBMaximum, 'f', 2);
		pBWhite->setText(sNum);
		pBWhiteScroll->setValue(100 * pImage->pImage->pFilters->dBMaximum / GetMaxValue());
	}
}

void TMOGUIFilters::component()
{
	pBlack->setEnabled(false);
	pWhite->setEnabled(false);
	pBlackScroll->setEnabled(false);
	pWhiteScroll->setEnabled(false);

	pRBlack->setEnabled(true);
	pGBlack->setEnabled(true);
	pBBlack->setEnabled(true);
	pRBlackScroll->setEnabled(true);
	pGBlackScroll->setEnabled(true);
	pBBlackScroll->setEnabled(true);

	pRWhite->setEnabled(true);
	pGWhite->setEnabled(true);
	pBWhite->setEnabled(true);
	pRWhiteScroll->setEnabled(true);
	pGWhiteScroll->setEnabled(true);
	pBWhiteScroll->setEnabled(true);

	pLuminance->setEnabled(false);
	pRed->setEnabled(true);
	pGreen->setEnabled(true);
	pBlue->setEnabled(true);
	if (pImage)
	{
		pImage->pImage->pFilters->dBAverage = 1.0;
		if (pImage->pImage->pFilters->dAverage > 0.0)
			pImage->pAdjust->pValues->RenderValues();
	}
}

void TMOGUIFilters::luminance()
{
	pBlack->setEnabled(true);
	pWhite->setEnabled(true);
	pBlackScroll->setEnabled(true);
	pWhiteScroll->setEnabled(true);

	pRBlack->setEnabled(false);

	pGBlack->setEnabled(false);
	pBBlack->setEnabled(false);
	pRBlackScroll->setEnabled(false);
	pGBlackScroll->setEnabled(false);
	pBBlackScroll->setEnabled(false);

	pRWhite->setEnabled(false);
	pGWhite->setEnabled(false);
	pBWhite->setEnabled(false);
	pRWhiteScroll->setEnabled(false);
	pGWhiteScroll->setEnabled(false);
	pBWhiteScroll->setEnabled(false);

	pLuminance->setEnabled(true);
	pRed->setEnabled(false);
	pGreen->setEnabled(false);
	pBlue->setEnabled(false);

	if (pImage)
	{
		pImage->pImage->pFilters->dBAverage = -1.0;
		if (pImage->pImage->pFilters->dAverage > 0.0)
			pImage->pAdjust->pValues->RenderValues();
	}
}

void TMOGUIFilters::color()
{
	if (pImage)
	{
		int iRed = (int)(pImage->pImage->pFilters->dRGamma * 255.0);
		int iGreen = (int)(pImage->pImage->pFilters->dGGamma * 255.0);
		int iBlue = (int)(pImage->pImage->pFilters->dBGamma * 255.0);
		pPixmap->fill(QColor(iRed, iGreen, iBlue));
		pColor->setIcon(QIcon(*pPixmap));
		pSlider->setValue((int)(pImage->pImage->pFilters->dExtreme * 100.0));
		if (pImage->pImage->pFilters->dAverage > 0.0)
			pEnable->setChecked(true);
		else
			pEnable->setChecked(false);
		if (pImage->pImage->pFilters->dRAverage > 0.0)
			pOutside->setChecked(true);
		else
			pOutside->setChecked(false);
		if (pImage->pImage->pFilters->dGAverage > 0.0)
			pConnect->setChecked(true);
		else
			pConnect->setChecked(false);

		if (pImage->pImage->pFilters->dAverage > 0.0)
			pImage->pAdjust->pValues->RenderValues();
	}
}

void TMOGUIFilters::connectvalues(bool bOn)
{
	if (pImage)
	{
		if (bOn)
			pImage->pImage->pFilters->dGAverage = 1.0;
		else
			pImage->pImage->pFilters->dGAverage = -1.0;
		valueschanged();
	}
}

void TMOGUIFilters::blackslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dMinimum = GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dMinimum, 'f', 2);
		pBlack->setText(sNum);
	}
}

void TMOGUIFilters::whiteslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dMaximum = GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dMaximum, 'f', 2);
		pWhite->setText(sNum);
	}
}

void TMOGUIFilters::rblackslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dRMinimum =
			GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dRMinimum, 'f', 2);
		pRBlack->setText(sNum);
	}
}

void TMOGUIFilters::rwhiteslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dRMaximum =
			GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dRMaximum, 'f', 2);
		pRWhite->setText(sNum);
	}
}

void TMOGUIFilters::gblackslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dGMinimum =
			GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dGMinimum, 'f', 2);
		pGBlack->setText(sNum);
	}
}

void TMOGUIFilters::gwhiteslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dGMaximum =
			GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dGMaximum, 'f', 2);
		pGWhite->setText(sNum);
	}
}

void TMOGUIFilters::bblackslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dBMinimum =
			GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dBMinimum, 'f', 2);
		pBBlack->setText(sNum);
	}
}

void TMOGUIFilters::bwhiteslider(int pos)
{
	if (pImage)
	{
		QString sNum;
		pImage->pImage->pFilters->dBMaximum =
			GetMaxValue() * pos / 100.0;
		sNum.setNum(pImage->pImage->pFilters->dBMaximum, 'f', 2);
		pBWhite->setText(sNum);
	}
}

double TMOGUIFilters::GetMaxValue()
{
	double dMax = 0;
	if (pImage)
		if (pImage->pAdjust->bLog)
			dMax = exp(pImage->pAdjust->pValues->dMaximum);
		else
			dMax = pImage->pAdjust->pValues->dMaximum;
	return dMax;
}

void TMOGUIFilters::updatesliders()
{
	if (pImage)
	{
		pImage->pImage->pFilters->dMinimum = pBlack->text().toDouble();
		pBlackScroll->setValue(100 * pImage->pImage->pFilters->dMinimum / GetMaxValue());

		pImage->pImage->pFilters->dRMinimum = pRBlack->text().toDouble();
		pRBlackScroll->setValue(100 * pImage->pImage->pFilters->dRMinimum / GetMaxValue());

		pImage->pImage->pFilters->dGMinimum = pGBlack->text().toDouble();
		pGBlackScroll->setValue(100 * pImage->pImage->pFilters->dGMinimum / GetMaxValue());

		pImage->pImage->pFilters->dBMinimum = pBBlack->text().toDouble();
		pBBlackScroll->setValue(100 * pImage->pImage->pFilters->dBMinimum / GetMaxValue());

		pImage->pImage->pFilters->dMaximum = pWhite->text().toDouble();
		pWhiteScroll->setValue(100 * pImage->pImage->pFilters->dMaximum / GetMaxValue());

		pImage->pImage->pFilters->dRMaximum = pRWhite->text().toDouble();
		pRWhiteScroll->setValue(100 * pImage->pImage->pFilters->dRMaximum / GetMaxValue());

		pImage->pImage->pFilters->dGMaximum = pGWhite->text().toDouble();
		pGWhiteScroll->setValue(100 * pImage->pImage->pFilters->dGMaximum / GetMaxValue());

		pImage->pImage->pFilters->dBMaximum = pBWhite->text().toDouble();
		pBWhiteScroll->setValue(100 * pImage->pImage->pFilters->dBMaximum / GetMaxValue());

		if (pImage->pImage->pFilters->dAverage > 0.0)
			pImage->pAdjust->pValues->RenderValues();
	}
}

void TMOGUIFilters::outsidevalues(bool bOn)
{
	if (pImage)
	{
		if (bOn)
			pImage->pImage->pFilters->dRAverage = 1.0;
		else
			pImage->pImage->pFilters->dRAverage = -1.0;
		if (pImage->pImage->pFilters->dAverage > 0.0)
			pImage->pAdjust->pValues->RenderValues();
	}
}

void TMOGUIFilters::intensity()
{
	int iValue = pSlider->value();
	if (pImage)
	{
		pImage->pImage->pFilters->dExtreme = (double)iValue / 100.0;

		if (pImage->pImage->pFilters->dAverage > 0.0)
			pImage->pAdjust->pValues->RenderValues();
	}
}

void TMOGUIFilters::colorchange()
{
	if (pImage)
	{
		QColor c;
		int iRed = (int)(pImage->pImage->pFilters->dRGamma * 255.0);
		int iGreen = (int)(pImage->pImage->pFilters->dGGamma * 255.0);
		int iBlue = (int)(pImage->pImage->pFilters->dBGamma * 255.0);
		c.setRgb(iRed, iGreen, iBlue);
		c = QColorDialog::getColor(c);
		pImage->pImage->pFilters->dRGamma = (double)c.red() / 255.0;
		pImage->pImage->pFilters->dGGamma = (double)c.green() / 255.0;
		pImage->pImage->pFilters->dBGamma = (double)c.blue() / 255.0;
		color();
	}
}

void TMOGUIFilters::enablefilters(bool bEnable)
{
	if (pImage)
	{
		if (bEnable)
			pImage->pImage->pFilters->dAverage = 1.0;
		else
			pImage->pImage->pFilters->dAverage = -1.0;
		pImage->pAdjust->pValues->RenderValues();
	}
}

void TMOGUIFilters::resetvalues()
{
	if (pImage)
	{
		pImage->pImage->pFilters->dBAverage = 0.0;	// Component
		pImage->pImage->pFilters->dRGamma = 1.0;	// Color
		pImage->pImage->pFilters->dGGamma = 0.0;	// .
		pImage->pImage->pFilters->dBGamma = 0.0;	// .
		pImage->pImage->pFilters->dMinimum = 0.0;	// Luminance
		pImage->pImage->pFilters->dMaximum = 0.0;	// .
		pImage->pImage->pFilters->dAverage = -1.0;	// enabled
		pImage->pImage->pFilters->dRAverage = -1.0; // outside
		pImage->pImage->pFilters->dGAverage = -1.0; // connected
		pImage->pImage->pFilters->dExtreme = 0.5;	// intensity

		valueschanged();
		refresh();
		color();
	}
}
