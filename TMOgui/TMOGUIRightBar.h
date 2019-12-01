// TMOGUIRightBar.h: interface for the TMOGUIRightBar class.
//
//////////////////////////////////////////////////////////////////////

#include <qtabwidget.h>
//Added by qt3to4:
#include <QResizeEvent>

class Q3VBox;
class QPushButton;
class TMOGUIToneMapping;
class TMOGUIFilters;
class TMO;

class TMOGUIRightBar : public QTabWidget  
{
	Q_OBJECT
public:
	TMOGUIRightBar( QWidget* parent = 0, const char* name = 0, Qt::WFlags f = 0 );
	virtual ~TMOGUIRightBar();
	virtual TMO* GetTMO();
	virtual TMOGUIToneMapping* GetMapping() {return pToneMapping;}

protected:

	void resizeEvent ( QResizeEvent * );
	Q3VBox* pVBox;
	QPushButton* pRightButton;
	TMOGUIFilters* pFilters;
	QTabWidget* pTabWidget;
	TMOGUIToneMapping* pToneMapping;
	bool bVisible;

friend class TMOGUIWindow;	

public slots:
	void hideright();

signals:
	void closeBar();
};
