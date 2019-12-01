// TMOGUIHisto.h: interface for the TMOGUIHisto class.
//
//////////////////////////////////////////////////////////////////////
#include <qwidget.h>
//Added by qt3to4:
#include <QPaintEvent>
#include <QResizeEvent>
#include <QPixmap>

#define HISTOGRAM_WIDTH 256
#define COMPONENTS_COUNT 3

class TMOGUIProgressBar;
class TMOGUIAdjustValues;
class TMOImage;
class QPixmap;

class TMOGUIHisto : public QWidget
{
	Q_OBJECT
public:
	virtual int SetWindow(TMOImage* pImage);
 	TMOGUIHisto (QWidget* parent, const char * name);
	virtual int Create(TMOImage* pImage, TMOGUIAdjustValues* pValues);
	virtual ~TMOGUIHisto();
public slots:
	void compute();
	void setscale(int);
	void setr();
	void setg();
	void setb();
	void setl();
	void setlinear();
	void setlog();
	void reset(int);
protected:
	double mapfrom(double);
	double mapto(double);
	bool bLog;
	unsigned int iMode;
	int iMaxCount;
	double dScale;
	unsigned int pLuminance[HISTOGRAM_WIDTH];
	unsigned int pComponents[COMPONENTS_COUNT][HISTOGRAM_WIDTH];
	void paintEvent ( QPaintEvent * pe);
	void resizeEvent ( QResizeEvent * pe);
	TMOImage* pSrc;
	TMOGUIAdjustValues* pValues; 
	QPixmap* pBackBuffer;
};
