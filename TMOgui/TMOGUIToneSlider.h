// TMOGUIToneSlider.h: interface for the TMOGUIToneSlider class.
//
//////////////////////////////////////////////////////////////////////
#include <qwidget.h>
//Added by qt3to4:
#include <QPixmap>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QPaintEvent>

class TMOImage;
class TMOGUIAdjustValues;
class QSize;
class QColor;
class QPixmap;


class TMOGUIToneSlider : public QWidget  
{
	Q_OBJECT
private:
	void correct_iGamma();
public:
 	TMOGUIToneSlider(QWidget* parent, const char * name);
	virtual int Create(TMOGUIAdjustValues* pValues);
	virtual ~TMOGUIToneSlider();
	void paintEvent ( QPaintEvent * pe);
	void resizeEvent ( QResizeEvent * re );
	void mousePressEvent ( QMouseEvent * me );
	void mouseReleaseEvent ( QMouseEvent * me );
	void mouseMoveEvent ( QMouseEvent * e );
	void mouseDoubleClickEvent ( QMouseEvent * e );
	void resetsliders();
	void setToWidth();
protected:
	double mapfrom(double);
	double mapto(double);
	int UpdateValues(bool bGamma = true);
	int DrawMarker(QPainter* p,int x, QColor col);
	bool InRect(const QPoint& p, int x1, int x2, int y1, int y2);
	int iMode;
	int iGrab;
	TMOGUIAdjustValues* pValues;
	QPixmap* pBackBuffer;
	int iBlack;
	int iWhite;
	int iGamma;
	int iCaptured;
	bool bLog;
	QSize s;
public slots:
	void setlinear();
	void setlog();
	void setr();
	void setg();
	void setb();
	void setl();
};
