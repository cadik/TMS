// TMOGUIParameters.cpp: implementation of the TMOGUIParameters class.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"
#include "TMOGUIParameters.h"
#include "TMOGUIParametersItem.h"
#include <qlayout.h>
#include <qsizepolicy.h>
#include <qlabel.h>
#include <qstring.h>

//Added by qt3to4:
#include <QResizeEvent>

#include "lqstring.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIParameters::TMOGUIParameters(QWidget *parent, const char *name) : QScrollArea(parent)
{
	parentWidget = parent;
	pParams = 0;
	iParams = 0;
	iCurParam = 0;
	backWidth = 0;
	//this->setWidgetResizable(false);
	this->setWidgetResizable(true);
	this->setAlignment(Qt::AlignTop);
	big_box = new QWidget(this->viewport()); //Q3VBox
	bigBoxLayout = new QVBoxLayout(big_box);
	bigBoxLayout->setAlignment(Qt::AlignTop);

	//big_box->setLayout(bigBoxLayout);
	//QRect r = this->contentsRect();

	//this->setWidget(big_box);
}

TMOGUIParameters::~TMOGUIParameters()
{
}

int TMOGUIParameters::SetTechnique(TMO *pTmo)
{
	TMOParameter **pParameters;
	int i;

	if (pTmo)
	{
		if (pParams && iParams)
		{
			for (i = 0; i < iParams; i++)
				delete pParams[i];
			delete[] pParams;
		}
		iParams = pTmo->GetParameterCount();
		if (bigBoxLayout)
		{
			delete bigBoxLayout;
			bigBoxLayout = new QVBoxLayout();
			bigBoxLayout->setAlignment(Qt::AlignTop);
		}
		if (iParams)
			pParams = new TMOGUIParametersItem *[iParams];
		if (iParams)
			pParameters = new TMOParameter *[iParams];
		if (iParams)
			pTmo->EnumParameters(pParameters);
		iCurParam = 0;
		for (i = 0; i < iParams; i++)
		{
			pParams[i] = new TMOGUIParametersItem(big_box);
			pParams[i]->Create(pParameters[i], this);
			//pParams[i]->setStyleSheet("background-color: blue;");
			bigBoxLayout->addWidget(pParams[i], Qt::AlignTop);
			connect(pParams[i], SIGNAL(change()), this, SLOT(valueschanged()));
			pParams[i]->show();
		}
		if (iParams)
			delete[] pParameters;
		bigBoxLayout->addStretch();
		update();
	}

	big_box->setLayout(bigBoxLayout);
	this->setWidget(big_box);
	return 0;
}

QString TMOGUIParameters::GetString(const wchar_t *s)
{
	//	QString str;
	lqstring str;
	str.setUnicodeCodes(s, wcslen(s));
	return str;
}

void TMOGUIParameters::resetvalues()
{
	for (int i = 0; i < iParams; i++)
		pParams[i]->resetvalues();
}

void TMOGUIParameters::viewportResizeEvent(QResizeEvent *e)
{
	QSize s = e->size();
	backWidth = s.width() - 63;
	for (int i = 0; i < iParams; i++)
		pParams[i]->ResizeElements(backWidth);
}

void TMOGUIParameters::valueschanged()
{
	emit changed();
}
