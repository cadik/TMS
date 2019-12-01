// TMOGUIStatistics.cpp: implementation of the TMOGUIStatistics class.
//
//////////////////////////////////////////////////////////////////////
#include <qlabel.h>
#include "TMOGUIStatistics.h"
#include "TMOGUIImage.h"
#include "TMOGUIAdjust.h"
#include "TMOGUIInfoTool.h"
#include "../tmolib/TMOImage.h"
#include "TMOGUIInfoPanel.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIStatistics::TMOGUIStatistics(QWidget* parent, const char * name):
	Q3ScrollView(parent, name)
{
	pPanel = new TMOGUIInfoPanel(viewport(), "Panel");
	addChild(pPanel);

	resetToolStats();
	setMinimumHeight(pPanel->height()+5);
}

TMOGUIStatistics::~TMOGUIStatistics()
{

}

void TMOGUIStatistics::windowChanged(TMOGUIImage* pWindow)
{
	TMOGUIImage* pOldImage = pImage;
	pImage = pWindow;
	if (pImage)
	{
		connect(pImage->pAdjust->pValues, SIGNAL(valueschanged()), this, SLOT(valueschanged()));
		connect(pImage->pImage, SIGNAL(rendered()), this, SLOT(valueschanged()));
		if(pImage->pImage->iTool)
		{
			connect(pImage->pImage->iTool, SIGNAL(toolApllied(int, int, int, bool)), this, SLOT(printToolStats(int, int, int, bool)));
			connect(pImage->pImage->iTool, SIGNAL(toolCancelled()), this, SLOT(resetToolStats()));
		}
		valueschanged();		
	}
}

void TMOGUIStatistics::valueschanged()
{
	QString sNum, sBuffer;

	if (pImage)
	{
		sNum.setNum(0.299 * pImage->pAdjust->pValues->dRMinimum + 0.587 * pImage->pAdjust->pValues->dGMinimum + 0.114 * pImage->pAdjust->pValues->dBMinimum, 'f', 2);
		sBuffer = sNum;
		sNum.setNum(pImage->pAdjust->pValues->dRMinimum, 'f', 2);
		sBuffer += " <" + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dGMinimum, 'f', 2);
		sBuffer += "; " + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dBMinimum, 'f', 2);
		sBuffer += "; " + sNum + ">";
		pPanel->pBlack->setText(sBuffer);

		sNum.setNum(0.299 * pImage->pAdjust->pValues->dRMaximum + 0.587 * pImage->pAdjust->pValues->dGMaximum + 0.114 * pImage->pAdjust->pValues->dBMaximum, 'f', 2);
		sBuffer = sNum;
		sNum.setNum(pImage->pAdjust->pValues->dRMaximum, 'f', 2);
		sBuffer += " <" + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dGMaximum, 'f', 2);
		sBuffer += "; " + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dBMaximum, 'f', 2);
		sBuffer += "; " + sNum + ">";
		pPanel->pWhite->setText(sBuffer);

		sNum.setNum(0.299 * pImage->pAdjust->pValues->dRGamma + 0.587 * pImage->pAdjust->pValues->dGGamma + 0.114 * pImage->pAdjust->pValues->dBGamma, 'f', 2);
		sBuffer = sNum;
		sNum.setNum(pImage->pAdjust->pValues->dRGamma, 'f', 2);
		sBuffer += " <" + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dGGamma, 'f', 2);
		sBuffer += "; " + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dBGamma, 'f', 2);
		sBuffer += "; " + sNum + ">";
		pPanel->pGamma->setText(sBuffer);

		sNum.setNum(0.299 * pImage->pAdjust->pValues->dRAverage + 0.587 * pImage->pAdjust->pValues->dGAverage + 0.114 * pImage->pAdjust->pValues->dBAverage, 'f', 2);
		sBuffer = sNum;
		sNum.setNum(pImage->pAdjust->pValues->dRAverage, 'f', 2);
		sBuffer += " <" + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dGAverage, 'f', 2);
		sBuffer += "; " + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dBAverage, 'f', 2);
		sBuffer += "; " + sNum + ">";
		pPanel->pAverage->setText(sBuffer);

		sNum.setNum(0.299 * pImage->pAdjust->pValues->dRGamma + 0.587 * pImage->pAdjust->pValues->dGGamma + 0.114 * pImage->pAdjust->pValues->dBGamma, 'f', 2);
		sBuffer = sNum;
		sNum.setNum(pImage->pAdjust->pValues->dRGamma, 'f', 2);
		sBuffer += " <" + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dGGamma, 'f', 2);
		sBuffer += "; " + sNum;
		sNum.setNum(pImage->pAdjust->pValues->dBGamma, 'f', 2);
		sBuffer += "; " + sNum + ">";
		pPanel->pGamma->setText(sBuffer);

		sNum.setNum((double)pImage->pImage->iBelow * 100.0 / pImage->pImage->iTotal, 'f', 2);
		sBuffer = sNum + "%";
		pPanel->pBelow->setText(sBuffer);

		sNum.setNum((double)pImage->pImage->iAbove * 100.0 / pImage->pImage->iTotal, 'f', 2);
		sBuffer = sNum + "%";
		pPanel->pAbove->setText(sBuffer);

		sNum.setNum(100.0 - (pImage->pImage->iAbove + pImage->pImage->iBelow) * 100.0 / pImage->pImage->iTotal, 'f', 2);
		sBuffer = sNum + "%";
		pPanel->pVisible->setText(sBuffer);

		sNum.setNum(pImage->pAdjust->pValues->dExtreme, 'f', 2);
		sBuffer = sNum;
		pPanel->pExtreme->setText(sBuffer);				
	}
}

void TMOGUIStatistics::printToolStats(int x, int y, int size, bool shape)
{
	QString sNum, sBuffer;
	double ratio = pImage->pImage->GetRatioD();
	x /= ratio;
	y /= ratio;
	if(shape)
		pImage->pImage->pSrc->GetStatisticsCircleRange(x, y, size);	
	else
		pImage->pImage->pSrc->GetStatisticsSquareRange(x, y, size);	

	sNum.setNum(pImage->pImage->pSrc->statistics.avgLum, 'f', 4);
	pPanel->pAvgLum->setText(sNum);
	sNum.setNum(pImage->pImage->pSrc->statistics.minLum, 'f', 4);
	pPanel->pMinLum->setText(sNum);
	sNum.setNum(pImage->pImage->pSrc->statistics.maxLum, 'f', 4);
	pPanel->pMaxLum->setText(sNum);
	sBuffer = "<";
	sNum.setNum(pImage->pImage->pSrc->statistics.avgColor[0], 'f', 2);
	sBuffer += sNum + "; ";
	sNum.setNum(pImage->pImage->pSrc->statistics.avgColor[1], 'f', 2);
	sBuffer += sNum + "; ";
	sNum.setNum(pImage->pImage->pSrc->statistics.avgColor[2], 'f', 2);
	sBuffer += sNum + ">";
	pPanel->pAvgCol->setText(sBuffer);
	sBuffer = "<";
	sNum.setNum(pImage->pImage->pSrc->statistics.color[0], 'f', 2);
	sBuffer += sNum + "; ";
	sNum.setNum(pImage->pImage->pSrc->statistics.color[1], 'f', 2);
	sBuffer += sNum + "; ";
	sNum.setNum(pImage->pImage->pSrc->statistics.color[2], 'f', 2);
	sBuffer += sNum + ">";
	pPanel->pColor->setText(sBuffer);
	sNum.setNum(pImage->pImage->pSrc->statistics.curLum, 'f', 4);
	pPanel->pLum->setText(sNum);
}

void TMOGUIStatistics::resetToolStats()
{
	pPanel->pAvgCol->setText(" ");
	pPanel->pAvgLum->setText(" ");
	pPanel->pLum->setText(" ");
	pPanel->pMinLum->setText(" ");
	pPanel->pMaxLum->setText(" ");
	pPanel->pColor->setText(" ");
}