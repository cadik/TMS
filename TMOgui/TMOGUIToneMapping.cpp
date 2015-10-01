// TMOGUIToneMapping.cpp: implementation of the TMOGUIToneMapping class.
//
//////////////////////////////////////////////////////////////////////
#ifndef LINUX
 #include "../TMOW32/TMOW32.h"
#else
 #include "../TMOLinux/TMOLinux.h"
#endif
#include "../tmolib/TMO.h"
#include "TMOGUIToneMapping.h"
#include "TMOGUIParameters.h"
#include <qlayout.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qcombobox.h>
#include <qscrollview.h>
#include <qpushbutton.h>
#include <qmultilineedit.h>
#include "lqstring.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIToneMapping::TMOGUIToneMapping( QWidget* parent, const char* name, WFlags f ):
	QWidget(parent, name, f)
{
	QGroupBox *pGroupBox;
	QLabel* pLabel;

	pTMO = 0;
	pParameters = 0;
	QGridLayout* pLayout = new QGridLayout(this,8,7);

	pLayout->addColSpacing(0,10);
	pLayout->addColSpacing(3,10);
	pLayout->addColSpacing(6,10);
	pLayout->addRowSpacing(1,10);
	pLayout->addRowSpacing(3,10);
	pLayout->addRowSpacing(5,10);
	pLayout->addRowSpacing(7,10);
	pLayout->setRowStretch(2,1);
	pLayout->setRowStretch(4,5);
	pLayout->setColStretch(1,1);
	pLayout->setColStretch(5,1);
	
	pGroupBox = new QGroupBox(1,Horizontal,"Tone Mapping",this, "ToneMappingGroupBox");
	pGroupBox->move(10,0);
	pLayout->addMultiCellWidget(pGroupBox,0,0,1,5);
	
	pLabel = new QLabel("Library",pGroupBox);
	pLibrary = new QComboBox( false, pGroupBox, "Library" );
	pLabel = new QLabel("Technique",pGroupBox);
	pTechnique = new QComboBox( false, pGroupBox, "Technique" );
	FillLibrary();
	
	pDescription = new QMultiLineEdit(this, "Description");
	pDescription->setReadOnly(true);
	pLayout->addMultiCellWidget(pDescription, 2,2,1,5);

	pParameters = new TMOGUIParameters(this, "Parameters");
	pLayout->addMultiCellWidget(pParameters, 4,4,1,5);
	ChangeTechnique(0);

	pOk = new QPushButton("OK", this, "RightOkButton");
	pOk->setFixedSize(64, 24);
	pLayout->addWidget(pOk, 6, 2);

	QPushButton *pCancel = new QPushButton("Reset", this, "ResetButton");
	pCancel->setFixedSize(64, 24);
	pLayout->addWidget(pCancel, 6, 4);

	connect(pLibrary, SIGNAL(activated(int)), this, SLOT(FillTechnique(int)));
	connect(pTechnique, SIGNAL(activated(int)), this, SLOT(ChangeTechnique(int)));
	connect(pCancel, SIGNAL(clicked()), pParameters, SLOT(resetvalues()));
}

TMOGUIToneMapping::~TMOGUIToneMapping()
{

}

int TMOGUIToneMapping::FillLibrary()
{
	int iLibCount;
	wchar_t* filenames[256];
	lqstring s;

	iLibCount = EnumLibraries(filenames, 256);
	pLibrary->clear();
	if (!iLibCount) 
	{
		pTMO = 0;
		s.setUnicodeCodes(L"No library available", 20);
		pLibrary->insertItem(s);
	}
	else
	{
		sLibraries = new wchar_t*[iLibCount];
		for (int i = 0; i < iLibCount; i++)
		{
		        int j;
			for (j = 0; filenames[i][j] != L'.'; j++);
			s.setUnicodeCodes(filenames[i], j);
			for (; filenames[i][j] != 0; j++);
			sLibraries[i] = new wchar_t[j+1];
			wcscpy(sLibraries[i], filenames[i]);
			pLibrary->insertItem(s);
		}
		iCurLibrary = 0;
		FillTechnique(0);
	}
	return 0;
}

void TMOGUIToneMapping::FillTechnique(int index)
{
	lqstring s;
	int i;
	int iTechCount = 0, iOpCount = 0;
	
	if (pTMO)
	{
		CloseLibrary(sLibraries[iCurLibrary], pTMO);
		delete[] pTMO;
	}
	pTechnique->clear();
	iOpCount = OperatorCount(sLibraries[index]);
	if (iOpCount) pTMO = new TMO*[iOpCount];
	else pTMO = 0;
	for (i = 0; i < iOpCount; i++) pTMO[i] = 0;
	if (pTMO)
	{
		iTechCount = OpenLibrary(sLibraries[index], pTMO);
		for (i = 0; i < iTechCount; i++)
		{
			try
			{
				s.setUnicodeCodes(pTMO[i]->GetName(), wcslen(pTMO[i]->GetName()));
			}
			catch(...)
			{
				// bad plugin format // plugin version mismatch
			}
			pTechnique->insertItem(s);
		}
		if (!iTechCount) 
		{
			s.setUnicodeCodes(L"No technique found.", 19);
			pTechnique->insertItem(s);
			ChangeTechnique(0);
		}
		else 
		{
			iCurLibrary = index;
			iCurTechnique = 0;
			ChangeTechnique(0);
		}
	}
	else
	{
		iTechCount = 0;
		s.setUnicodeCodes(L"No technique found.", 18);
		pTechnique->insertItem(s);
		ChangeTechnique(0);
	}
}

void TMOGUIToneMapping::ChangeTechnique(int index)
{
	if (pParameters == 0) return;
	if (index < 0 || pTMO == 0) 
	{
		iCurTechnique = 0;
		pDescription->setText(QString(0));
		pParameters->SetTechnique(0);
	}
	else 
	{
		iCurTechnique = index;
		if (pTMO[index])
		{
			const wchar_t* temp = pTMO[index]->GetDescription();
			pDescription->setText(pParameters->GetString(temp));
			pParameters->SetTechnique(pTMO[index]);
		}
		else
		{
			pDescription->setText("Failed to retrieve a description");
			pParameters->SetTechnique(0);
		}
	}
}
