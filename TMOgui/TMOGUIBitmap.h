// TMOGUIBitmap.h: interface for the TMOGUIBitmap class.
//
//////////////////////////////////////////////////////////////////////
#include <qwidget.h>
//Added by qt3to4:
#include <QPixmap>
#include <QEvent>
#include <QMouseEvent>
#include <QContextMenuEvent>
#include <QPaintEvent>

class QPixmap;
class QSize;
class TMOImage;
class TMOGUIProgressBar;
class TMOGUIInfoTool;
class QString;
class TMOGUIAdjustValues;

class TMOGUIBitmap : public QWidget  
{
	Q_OBJECT
public:
	TMOImage* GetImage(TMOImage* pDst);
	virtual QSize* GetSize();
	virtual int SetSize(const QSize* s = 0);
	virtual int SetImage(TMOImage* pImage);
	virtual int AddString(const QString& s);
	virtual int GetRatio();
	virtual double GetRatioD() const { return dRatio; }
	virtual int Zoom(double ratio);
	virtual int Create(TMOImage* src, TMOGUIAdjustValues* pValues, TMOGUIAdjustValues* pFilters, QSize window, TMOGUIProgressBar* pProgress);
	virtual int Render(bool bRepaint = true);
	virtual int DrawIt(QPainter* pPainter, double *pMargins, bool bZoom);
	virtual void ActivateTool ( TMOGUIInfoTool * ptr );
	virtual void DeactivateTool();
	TMOGUIBitmap(QWidget * parent=0, const char * name=0);
	virtual ~TMOGUIBitmap();
	unsigned int iAbove;
	unsigned int iBelow;	
	unsigned int iTotal;
	TMOGUIAdjustValues* pFilters;
	TMOGUIInfoTool* iTool;
	TMOImage* pSrc;
	QPixmap* pSrcPixmap;

protected:
	virtual void mouseMoveEvent ( QMouseEvent * e );
	virtual void mousePressEvent ( QMouseEvent * e );
	virtual void contextMenuEvent( QContextMenuEvent * );
	virtual void leaveEvent ( QEvent * );
	QString sInitial;
	int iInitial;	
	TMOGUIProgressBar* pProgress;
	double dRatio;
	QSize s;	
	QPixmap* pPixmap;
	void paintEvent ( QPaintEvent * );
	TMOGUIAdjustValues* pValues;	
	bool bRendering;

public slots:
	virtual int SetInitial(int i);
	void valueschanged();

signals:
	void rendered();	
};

