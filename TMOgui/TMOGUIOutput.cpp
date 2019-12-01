// TMOGUIOutput.cpp: implementation of the TMOGUIOutput class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIOutput.h"
#include "TMOGUIInfo.h"
#include "lqstring.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
QMap<TMOImage*, TMOGUIOutput*> TMOGUIOutput::mapLocal;
TMOGUIInfo* TMOGUIOutput::pInfo = 0;

TMOGUIOutput::TMOGUIOutput(QWidget* parent, const char * name):
	Q3MultiLineEdit(parent, name)
{
	bVisible = false;
	sText = "";
	uiLines = 0;
}

TMOGUIOutput::~TMOGUIOutput()
{
	QMap<TMOImage*, TMOGUIOutput*>::Iterator i;

	for (i = mapLocal.begin(); i != mapLocal.end(); i++) 
		if (i.data() == this) 
		{
			mapLocal.remove(i);
			break;
		}
}

int TMOGUIOutput::Assign(TMOImage *pImage)
{
	QMap<TMOImage*, TMOGUIOutput*>::Iterator i;

	for (i = mapLocal.begin(); i != mapLocal.end(); i++)
		if (i.data() == this)
		{
			mapLocal.remove(i);
			break;
		}
	mapLocal.insert(pImage, this);
	return 0;
}

int TMOGUIOutput::WriteLine(TMOImage* pImage, const wchar_t* text)
{
	QMap<TMOImage*, TMOGUIOutput*>::Iterator i;
	TMOGUIOutput* pLocal;
	lqstring s;

	i = mapLocal.find(pImage);
	if (i == mapLocal.end()) return 0;
	pLocal = i.data();

	if (!pLocal) return 0;

	if (!pLocal->isVisible()) 
	{
	}

	pLocal->bVisible = true;

	s.setUnicodeCodes(text, TextLength(text));
	pLocal->sText.append(s);
	s.setUnicodeCodes(L"\n", 1);
	pLocal->sText.append(s);
	pLocal->setText(pLocal->sText);
	pLocal->uiLines++;
	pLocal->setCursorPosition(pLocal->uiLines, 0);
	return 0;
}

int TMOGUIOutput::TextLength(const wchar_t *text)
{
	int i=0;
	for (i = 0; i < 256; i++)
		if (!text[i]) break;
	return i;
}
