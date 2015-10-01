// TMOGUIFilters.h: interface for the TMOGUIFilters class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIAdjustValues.h"
#include <qwidget.h>
class QComboBox;
class TMOGUIImage;
class QPushButton;
class QLineEdit;
class QLabel;
class QGridLayout;
class QCheckBox;
class QPixmap;
class QSlider;
class QScrollBar;

class TMOGUIFilters : public QWidget 
{
	Q_OBJECT
public:
	TMOGUIFilters(QWidget* parent = 0, const char* name = 0, WFlags f = 0 );
	virtual ~TMOGUIFilters();

protected:
	int Create();
	double GetMaxValue();

	QLineEdit *pWhite, *pRWhite, *pGWhite, *pBWhite;
	QLineEdit *pBlack, *pRBlack, *pGBlack, *pBBlack;

	QScrollBar *pWhiteScroll, *pRWhiteScroll, *pGWhiteScroll, *pBWhiteScroll;
	QScrollBar *pBlackScroll, *pRBlackScroll, *pGBlackScroll, *pBBlackScroll;

	QLabel *pLuminance, *pRed, *pGreen, *pBlue;
	TMOGUIImage *pImage;
	QGridLayout* pLayout;
	QPixmap* pPixmap;
	QPushButton *pColor;
	QPushButton *pEnable;
	QCheckBox *pOutside;
	QCheckBox *pConnect;
	QSlider *pSlider;
	
public slots:
	void windowChanged(TMOGUIImage* pWindow);
	void valueschanged();
	void component();
	void luminance();
	void refresh();
	void color();
	void connectvalues(bool bOn);
	void blackslider(int pos);
	void whiteslider(int pos);
	void rblackslider(int pos);
	void rwhiteslider(int pos);
	void gblackslider(int pos);
	void gwhiteslider(int pos);
	void bblackslider(int pos);
	void bwhiteslider(int pos);
	void outsidevalues(bool bOn);
	void updatesliders();
	void intensity();
	void colorchange();
	void enablefilters(bool bEnable);
	void resetvalues();
};




