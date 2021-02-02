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
#include "tmoguiwindow.h"
#include "TMOGUIToneMappingChooser.h"
#include "TMOGUIParameters.h"
#include "TMOGUIBitmap.h"
#include <qlayout.h>
#include <qlabel.h>
#include <QGroupBox>
#include <qcombobox.h>
#include <QScrollArea>
#include <qpushbutton.h>
#include <QTextEdit>
//Added by qt3to4:
#include <QGridLayout>
#include "lqstring.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIToneMapping::TMOGUIToneMapping( QWidget* parent, const char* name/*, Qt::WindowFlags f */):
    QWidget(parent)
{
    int i,j,ops = 0;
	QLabel* pLabel;

    //setAutoFillBackground(true);

	pTMO = 0;
	pParameters = 0;

    vBoxLayout = new QVBoxLayout(this);
    pParametersParent = new QWidget(this);
    gridLayout = new QGridLayout(pParametersParent); //,8,7);

    gridLayout->addItem(new QSpacerItem(10,0), 0, 0); //pLayout->addColSpacing(0,10);
    gridLayout->addItem(new QSpacerItem(10,0), 0, 3); //pLayout->addColSpacing(3,10);
    gridLayout->addItem(new QSpacerItem(10,0), 0, 6); //pLayout->addColSpacing(6,10);
    gridLayout->addItem(new QSpacerItem(0,10), 1, 0); //pLayout->addRowSpacing(1,10);
    gridLayout->addItem(new QSpacerItem(0,10), 3, 0); //pLayout->addRowSpacing(3,10);
    gridLayout->addItem(new QSpacerItem(0,10), 5, 0); //pLayout->addRowSpacing(5,10);
    gridLayout->addItem(new QSpacerItem(0,10), 7, 0); //pLayout->addRowSpacing(7,10);
    gridLayout->setRowStretch(2,1);
    gridLayout->setRowStretch(4,5);
    gridLayout->setColumnStretch(1,1);
    gridLayout->setColumnStretch(5,1);
	
    pGroupBox = new QGroupBox("Tone Mapping",pParametersParent);//, "ToneMappingGroupBox");
    pGroupBox->setObjectName("ToneMappingGroupBox");

    // strips 1, orientation Qt::Horizontal,
	pGroupBox->move(10,0);
    QHBoxLayout* pGroupBoxLayout = new QHBoxLayout();


    //pLayout->addMultiCellWidget(pGroupBox,0,0,1,5);
	
	pLabel = new QLabel("Library",pGroupBox);
    pGroupBoxLayout->addWidget(pLabel);
    pLibrary = new QComboBox( pGroupBox); //false,
    pLibrary->setMinimumWidth(100);
    pLibrary->setObjectName("Library");
    pGroupBoxLayout->addWidget(pLibrary);
	pLabel = new QLabel("Technique",pGroupBox);
    pGroupBoxLayout->addWidget(pLabel);
    pTechnique = new QComboBox(  pGroupBox); //false,
    pTechnique->setMinimumWidth(100);
    pTechnique->setObjectName( "Technique" );
    pGroupBoxLayout->addWidget(pTechnique);

	FillLibrary();

    pChooser = new TMOGUIToneMappingChooser(this, nullptr);

    for(i = 0; i < pLibrary->count(); i++)
    {
        ops = OperatorCount(sLibraries[i]);
        TMO** tempTMO = new TMO*[ops];
        OpenLibrary(sLibraries[i], tempTMO);
        for(j = 0; j < ops; j++)
        {
            pChooser->AddTMOPreview(tempTMO[j], i, j);
        }
    }
    mChoosers.insert("default-preview", pChooser);

    FillTechnique(0);

    pGroupBox->setLayout(pGroupBoxLayout);
    gridLayout->addWidget(pGroupBox, 0, 1, 1, 5);

    pBackButton = new QPushButton("Select TMO", this);
    pBackButton->setHidden(true);

    pDescription = new QTextEdit(pParametersParent);//, "Description");
    pDescription->setObjectName("Description");
	pDescription->setReadOnly(true);
    gridLayout->addWidget(pDescription, 2, 1, 1, 5);
    //pLayout->addMultiCellWidget(pDescription, 2,2,1,5);

    pParameters = new TMOGUIParameters(pParametersParent, "Parameters");
    gridLayout->addWidget(pParameters, 4, 1, 1, 5);
    //pLayout->addMultiCellWidget(pParameters, 4,4,1,5);
	ChangeTechnique(0);
    connect(pParameters, SIGNAL(changed()), this, SLOT(paramChanged()));

    pOk = new QPushButton("OK", pParametersParent);
    pOk->setObjectName("RightOkButton");
    pOk->setFixedSize(128, 24);
    gridLayout->addWidget(pOk, 6, 2);

    pPreview = new QPushButton("Preview", pParametersParent);
    pPreview->setObjectName("RightOkButton");
    pPreview->setFixedSize(64, 24);
    gridLayout->addWidget(pPreview, 6, 3);

    QPushButton *pCancel = new QPushButton("Reset", pParametersParent);
    pCancel->setObjectName( "ResetButton");
	pCancel->setFixedSize(64, 24);
    gridLayout->addWidget(pCancel, 6, 4);

    vBoxLayout->addWidget(pBackButton, Qt::AlignTop);
    vBoxLayout->addWidget(pChooser, Qt::AlignTop);
    vBoxLayout->addWidget(pParametersParent, Qt::AlignTop);

    this->setLayout(vBoxLayout);
    pChooser->displayAll();
    pChooser->setHidden(true);

    connect(pLibrary, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIToneMapping::FillTechnique);
    connect(pTechnique, QOverload<int>::of(&QComboBox::activated), this, &TMOGUIToneMapping::ChangeTechnique);
	connect(pCancel, SIGNAL(clicked()), pParameters, SLOT(resetvalues()));
    connect(pBackButton, SIGNAL(clicked()), this, SLOT(toggleTechniqueChooser()));
}


TMOGUIToneMapping::~TMOGUIToneMapping()
{
    if(pChooser) delete pChooser;
    mChoosers.clear();
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
        pLibrary->addItem(s);// insertItem
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
            pLibrary->addItem(s);
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
	
    /* BUG closing libraries results in SEGFAULT when computing previews
    if (pTMO)
	{
		CloseLibrary(sLibraries[iCurLibrary], pTMO);
		delete[] pTMO;
    }*/
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
            pTechnique->addItem(s); //insertItem
		}
		if (!iTechCount) 
		{
			s.setUnicodeCodes(L"No technique found.", 19);
            pTechnique->addItem(s); //insertItem
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
        pTechnique->addItem(s); //insertItem
		ChangeTechnique(0);
	}
}

void TMOGUIToneMapping::ChangeTechnique(int index)
{
    if (pParameters == nullptr) return;
    if (index < 0 || pTMO == nullptr)
	{
		iCurTechnique = 0;
        pDescription->setText(QString("0"));
        pParameters->SetTechnique(nullptr);
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
            pParameters->SetTechnique(nullptr);
		}
	}
}

void TMOGUIToneMapping::toggleTechniqueChooser()
{
    bool chooserVisible = pChooser->isVisible();

    pChooser->setHidden(chooserVisible);
    pParametersParent->setHidden(!chooserVisible);

    QString name;
    pBackButton->setDown(!chooserVisible);
    if(!chooserVisible) pBackButton->setText("Select TMO");
    else pBackButton->setText(name.fromWCharArray(pTMO[iCurTechnique]->GetName()));

}

void TMOGUIToneMapping::getAllTechniques()
{
    /*lqstring s;
    int i;
    int iTechCount = 0, iOpCount = 0;

    if (pTMO)
    {
        CloseLibrary(sLibraries[iCurLibrary], pTMO);
        delete[] pTMO;
    }
    pTechnique->clear();
    for(i = 0; i < pLibrary->count(); i++){
        iOpCount = OperatorCount(sLibraries[i]);
    }

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
            pTechnique->addItem(s); //insertItem
        }
        if (!iTechCount)
        {
            s.setUnicodeCodes(L"No technique found.", 19);
            pTechnique->addItem(s); //insertItem
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
        pTechnique->addItem(s); //insertItem
        ChangeTechnique(0);
    }*/
}

void TMOGUIToneMapping::changeWorkspace(bool advanced)
{
    pGroupBox->setHidden(!advanced);
    pBackButton->setHidden(advanced);
    pChooser->setHidden(advanced);
    pParametersParent->setHidden(!advanced);
    pBackButton->setDown(!advanced);
}

void TMOGUIToneMapping::displayPreviews(TMOGUIImage *pImg)
{
    int i, j, ops;
    TMOGUIToneMappingChooser *oldChooser = pChooser;
    //TMOGUIImage *smallImage = nullptr;

    // If not new image, show saved chooser, else create new
    if(!pImg){
        pChooser = mChoosers.value("default-preview");
    } else if(mChoosers.contains(pImg->GetImage()->GetFilename())){
        pChooser = mChoosers.value(pImg->GetImage()->GetFilename());
    } else {

        pChooser = new TMOGUIToneMappingChooser(this, pImg);

        for(i = 0; i < pLibrary->count(); i++)
        {
            ops = OperatorCount(sLibraries[i]);
            TMO** tempTMO = new TMO*[ops];
            OpenLibrary(sLibraries[i], tempTMO);
            for(j = 0; j < ops; j++)
            {
                pChooser->AddTMOPreview(tempTMO[j], i, j);
            }
        }
        pChooser->displayAll();
        mChoosers.insert(pImg->GetImage()->GetFilename(), pChooser);
    }

    // Replace chooser if different
    if(oldChooser && oldChooser != pChooser)
    {
    pChooser->setHidden(oldChooser->isHidden());
    vBoxLayout->replaceWidget(oldChooser, pChooser);
    oldChooser->setHidden(true);
    }




}

void TMOGUIToneMapping::chooseTechnique(int indexLib, int indexTMO)
{
    FillTechnique(indexLib);
    ChangeTechnique(indexTMO);
    pLibrary->setCurrentIndex(indexLib);
    toggleTechniqueChooser();
}

void::TMOGUIToneMapping::paramChanged(){
    emit change();
}

void TMOGUIToneMapping::windowChanged(TMOGUIImage* pImage){
    if(!pImage){
        displayPreviews(nullptr);
        return;
    }

    if(pImage->bPreview){
        TMOGUIWindow* parentWin = (TMOGUIWindow*)(parent()->parent()->parent());
        TMOGUIImage* origImage = parentWin->GetActiveImage();
        displayPreviews(origImage);
    } else {
        //pChooser->windowChanged(pImage);
        displayPreviews(pImage);
    }

}
