
#ifndef TMOGUIINFOPANEL_H
#define TMOGUIINFOPANEL_H

#include <qwidget.h>

class QLabel;
class QGroupBox;
class TMOGUILineResizer;

/*! \class TMOGUIInfoPanel	
 *	\brief Displays statistics of the current image window.
 *
 * Consists of three group boxes with lineresizers (Resizes groupbox's size).
 */
class TMOGUIInfoPanel : public QWidget
{
    Q_OBJECT

public:
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
    TMOGUIInfoPanel( QWidget* parent = 0, const char* name = 0 );
	//! Destructor
    ~TMOGUIInfoPanel();

	QGroupBox* groupBoxHisto;
    QLabel* textLabelBlack;
	QLabel* textLabelWhite;
    QLabel* textLabelGammma;
	QLabel* textLabelAvg;
	QLabel* pBlack;
	QLabel* pWhite;
	QLabel* pGamma;
    QLabel* pAverage;
	TMOGUILineResizer* line1;
    QGroupBox* groupBoxStat;
    QLabel* textLabelBelow;
    QLabel* textLabelAbove;
    QLabel* textLabelVisible;
    QLabel* textLabelExtreme;
	QLabel* pBelow;	
	QLabel* pAbove;
    QLabel* pVisible;    
    QLabel* pExtreme;
	TMOGUILineResizer* line2;
    QGroupBox* groupBoxLocalTool;
	QLabel* textLabelNeighb;
	QLabel* textLabelAvgLum;
	QLabel* textLabelMinLum;
    QLabel* textLabelMaxLum;
	QLabel* textLabelAvgCol;
	QLabel* pAvgLum;
	QLabel* pMaxLum;
    QLabel* pMinLum;
	QLabel* pAvgCol;
    QLabel* ptextLabelCursor;
	QLabel* ptextLabelLum;
    QLabel* ptextLabelColor;
    QLabel* pLum;
    QLabel* pColor;
	TMOGUILineResizer* line3;

	virtual void languageChange();

public slots:
	void changeSizeInfoFirst(int);
	void changeSizeInfoSecond(int);
	void changeSizeInfoThird(int);

};

#endif // TMOGUIINFOPANEL_H 