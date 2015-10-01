// TMOGUICenterView.h: interface for the TMOGUICenterView class.
//
//////////////////////////////////////////////////////////////////////
#include <qscrollview.h>

class TMOGUICenterView : public QScrollView  
{
	Q_OBJECT
public:
	virtual void Zoom(bool bZoom = true);
	TMOGUICenterView(QWidget * parent=0, const char * name=0, WFlags f=0 );
	virtual ~TMOGUICenterView();
protected:
	virtual void viewportPaintEvent ( QPaintEvent * pe );
	bool bZooming;

};


