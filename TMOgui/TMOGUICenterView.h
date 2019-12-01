// TMOGUICenterView.h: interface for the TMOGUICenterView class.
//
//////////////////////////////////////////////////////////////////////
#include <q3scrollview.h>
//Added by qt3to4:
#include <QPaintEvent>

class TMOGUICenterView : public Q3ScrollView  
{
	Q_OBJECT
public:
	virtual void Zoom(bool bZoom = true);
	TMOGUICenterView(QWidget * parent=0, const char * name=0, Qt::WFlags f=0 );
	virtual ~TMOGUICenterView();
protected:
	virtual void viewportPaintEvent ( QPaintEvent * pe );
	bool bZooming;

};


