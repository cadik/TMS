// TMOGUICenterView.h: interface for the TMOGUICenterView class.
//
//////////////////////////////////////////////////////////////////////
#include <QScrollArea>
//Added by qt3to4:
#include <QPaintEvent>

class TMOGUICenterView : public QScrollArea
{
	Q_OBJECT
public:
	virtual void Zoom(bool bZoom = true);
    TMOGUICenterView(QWidget * parent=0, const char * name=0/*, Qt::WindowFlags f=0 */);
	virtual ~TMOGUICenterView();
protected:
	virtual void viewportPaintEvent ( QPaintEvent * pe );
	bool bZooming;

};


