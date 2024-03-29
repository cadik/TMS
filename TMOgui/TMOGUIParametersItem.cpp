// TMOGUIParametersItem.cpp: implementation of the TMOGUIParametersItem class.
//
//////////////////////////////////////////////////////////////////////
#include <qwidget.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qscrollbar.h>
#include <QScrollArea>
#include <qlayout.h>
#include <qtooltip.h>
#include <qcheckbox.h>
//Added by qt3to4:
#include <QGridLayout>
#include "../tmolib/TMO.h"
#include "TMOGUIParameters.h"
#include "TMOGUIParametersItem.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIParametersItem::TMOGUIParametersItem(QWidget *parent, const char *name) : QWidget(parent)
{
	pParameter = 0;
	pWidgets = 0;
	pLayout = 0;
}

TMOGUIParametersItem::~TMOGUIParametersItem()
{
	Destroy(0);
}

int TMOGUIParametersItem::Create(TMOParameter *pParam, TMOGUIParameters *pParentWidget)
{
	pParameter = pParam;
	if (pParam->Is(TMO_INT) || pParam->Is(TMO_DOUBLE))
	{
		pLayout = new QGridLayout(this);
		pLayout->addItem(new QSpacerItem(0, 10), 0, 0); //pLayout->setRowSpacing(0, 10);
		pLayout->addItem(new QSpacerItem(0, 5), 3, 0);	//pLayout->setRowSpacing(3, 5);
		pLayout->addItem(new QSpacerItem(5, 0), 0, 0);	//pLayout->setColSpacing(0, 5);
		pLayout->addItem(new QSpacerItem(5, 0), 0, 2);	//pLayout->setColSpacing(2, 5);
		pLayout->addItem(new QSpacerItem(5, 0), 0, 4);	//pLayout->setColSpacing(4, 5);
		pLayout->setColumnStretch(1, 20);
		pLayout->setColumnStretch(2, 1);
		pLayout->setColumnStretch(3, 8);
		pLayout->setRowStretch(1, 1);
		pLayout->setRowStretch(2, 1);
		QString s, s2;
		iWidgets = 3;
		pWidgets = new QWidget *[iWidgets];
		pWidgets[2] = new QScrollBar(Qt::Horizontal, this); //, "ParamScroll");
		QScrollBar *sb = static_cast<QScrollBar *>(pWidgets[2]);
		if (pParam->Is(TMO_INT))
		{
			int min, max;
			TMOInt *pI = static_cast<TMOInt *>(pParam);
			pWidgets[1] = new QLineEdit(s.setNum(pParam->GetInt()), this); //, "Parameter1");
			pI->GetRange(min, max);
			if ((max - min) > maxRange)
			{
				max = min + maxRange;
			}
			s.setNum(min);
			s2.setNum(max);
			sb->setRange(min, max);
		}
		if (pParam->Is(TMO_DOUBLE))
		{
			double min, max;
			TMODouble *pD = static_cast<TMODouble *>(pParam);
			pWidgets[1] = new QLineEdit(s.setNum(pParam->GetDouble()), this); //, "Parameter1");
			pD->GetRange(min, max);
			if ((max - min) > maxRange)
			{
				max = min + maxRange;
			}
			s.setNum(min, 'f', 2);
			s2.setNum(max, 'f', 2);
			sb->setRange(0, 100);
		}
		const wchar_t *temp = pParam->GetDescription();
		pWidgets[0] = new QLabel(TMOGUIParameters::GetString(pParam->GetName()) + " [" + s + ", " + s2 + "]", this); //, "Parameter");
		pWidgets[0]->setToolTip(pParentWidget->GetString(temp));
		pWidgets[0]->setFixedHeight(20);
		int index = pParentWidget->iCurParam * 2;
		//pWidgets[1]->setFixedWidth(48);
		pWidgets[1]->setFixedHeight(25);
		pWidgets[1]->setMinimumWidth(48);
		//pWidgets[1]->setMaximumWidth(50);
		//pWidgets[2]->setFixedWidth((pParentWidget->backWidth != 0) ? pParentWidget->backWidth : 130);
		pWidgets[2]->setFixedHeight(25);
		pLayout->addWidget(pWidgets[0], 1, 1, 1, 4, Qt::AlignTop);
		//pLayout->addMultiCellWidget(pWidgets[0], 1, 1, 1, 4);
		pLayout->addWidget(pWidgets[2], 2, 1);
		pLayout->addWidget(pWidgets[1], 2, 3);
		this->setFixedHeight(70);

		this->setLayout(pLayout);
		pWidgets[0]->show();
		pWidgets[1]->show();
		pWidgets[2]->show();
		pParentWidget->iCurParam++;
		resetvalues();
		connect((QLineEdit *)pWidgets[1], &QLineEdit::textChanged, this, QOverload<const QString &>::of(&TMOGUIParametersItem::valuechanged));
		connect((QScrollBar *)pWidgets[2], &QScrollBar::valueChanged, this, &TMOGUIParametersItem::scrollbarchanged);
		connect((QScrollBar *)pWidgets[2], &QScrollBar::sliderReleased, this, &TMOGUIParametersItem::paramChanged);
	}
	else if (pParam->Is(TMO_BOOL))
	{
		pLayout = new QGridLayout(this);				// , 2, 4);
		pLayout->addItem(new QSpacerItem(0, 10), 0, 0); //pLayout->setRowSpacing(0, 10);
		pLayout->addItem(new QSpacerItem(0, 5), 2, 0);	//pLayout->setRowSpacing(2, 5);
		QString s;
		QCheckBox *pCheck;
		iWidgets = 2;
		pWidgets = new QWidget *[iWidgets];
		const wchar_t *temp = pParam->GetDescription();
		pWidgets[0] = new QLabel(TMOGUIParameters::GetString(pParam->GetName()), this); //, "Parameter");
		pWidgets[0]->setToolTip(pParentWidget->GetString(temp));
		pWidgets[1] = pCheck = new QCheckBox(this); //, "Parameter1");
		pCheck->setChecked(pParam->GetBool());
		pLayout->addWidget(pWidgets[0], 1, 1);
		pLayout->addWidget(pWidgets[1], 1, 2);
		this->setLayout(pLayout);
		pWidgets[0]->show();
		pWidgets[1]->show();
		pParentWidget->iCurParam++;
		connect((QCheckBox *)pWidgets[1], QOverload<int>::of(&QCheckBox::stateChanged), this, QOverload<int>::of(&TMOGUIParametersItem::valuechanged));
	}

	return 0;
}

void TMOGUIParametersItem::valuechanged(const QString &text)
{
	if (pParameter->Is(TMO_INT))
	{
		*pParameter = text.toInt();
	}
	if (pParameter->Is(TMO_DOUBLE))
	{
		*pParameter = text.toDouble();
	}
}

void TMOGUIParametersItem::scrollbarchanged(int pos)
{
	QString s;
	QLineEdit *pEdit = static_cast<QLineEdit *>(pWidgets[1]);
	if (pParameter->Is(TMO_INT))
	{
		s.setNum(pos);
		*pParameter = pos;
	}
	if (pParameter->Is(TMO_DOUBLE))
	{
		double min, max;
		TMODouble *pD = static_cast<TMODouble *>(pParameter);
		pD->GetRange(min, max);
		if ((max - min) > maxRange)
		{
			max = min + maxRange;
		}
		double value = (max - min) * pos / 100.0;
		s.setNum(min + value, 'f', 2);
		*pParameter = value;
	}
	pEdit->setText(s);
}

void TMOGUIParametersItem::valuechanged(int state)
{
	if (pParameter->Is(TMO_BOOL))
	{
		*pParameter = state == 2;
		emit change();
	}
}

void TMOGUIParametersItem::paramChanged()
{
	emit change();
}

void TMOGUIParametersItem::resetvalues()
{
	QString s;

	pParameter->Reset();
	if (pParameter->Is(TMO_INT))
	{
		ResetScrollBar();
		QLineEdit *pEdit = static_cast<QLineEdit *>(pWidgets[1]);
		s.setNum(pParameter->GetInt());
		pEdit->setText(s);
	}
	else if (pParameter->Is(TMO_DOUBLE))
	{
		ResetScrollBar();
		QLineEdit *pEdit = static_cast<QLineEdit *>(pWidgets[1]);
		s.setNum(pParameter->GetDouble());
		pEdit->setText(s);
	}
	else if (pParameter->Is(TMO_BOOL))
	{
		QCheckBox *pCheck = static_cast<QCheckBox *>(pWidgets[1]);
		pCheck->setChecked(pParameter->GetBool());
	}

	emit change();
}

int TMOGUIParametersItem::Destroy(TMOGUIParameters *pParentWidget)
{
	for (int i = 0; i < iWidgets; i++)
	{
		//pLayout->invalidate();
		if (pWidgets && pWidgets[i])
			delete pWidgets[i];
	}
	if (pWidgets)
		delete[] pWidgets;
	if (pLayout)
		delete pLayout;
	return 0;
}

void TMOGUIParametersItem::ResizeElements(int width)
{
	if (pParameter->Is(TMO_INT) || pParameter->Is(TMO_DOUBLE))
		pWidgets[2]->setFixedWidth(width);
}

void TMOGUIParametersItem::ResetScrollBar()
{
	int pos;
	QScrollBar *pScb = static_cast<QScrollBar *>(pWidgets[2]);
	if (pParameter->Is(TMO_INT))
	{
		pos = pParameter->GetInt();
	}
	if (pParameter->Is(TMO_DOUBLE))
	{
		double min, max;
		TMODouble *pD = static_cast<TMODouble *>(pParameter);
		pD->GetRange(min, max);
		pos = (int)((pParameter->GetDouble() - min) / (max - min) * 100);
	}
	pScb->setValue(pos);
}
