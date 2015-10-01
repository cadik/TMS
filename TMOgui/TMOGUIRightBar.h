// TMOGUIRightBar.h: interface for the TMOGUIRightBar class.
//
//////////////////////////////////////////////////////////////////////

#include <qtabwidget.h>

class QVBox;
class QPushButton;
class TMOGUIToneMapping;
class TMOGUIFilters;
class TMO;

class TMOGUIRightBar : public QTabWidget  
{
	Q_OBJECT
public:
	TMOGUIRightBar( QWidget* parent = 0, const char* name = 0, WFlags f = 0 );
	virtual ~TMOGUIRightBar();
	virtual TMO* GetTMO();
	virtual TMOGUIToneMapping* GetMapping() {return pToneMapping;}

protected:

	void resizeEvent ( QResizeEvent * );
	QVBox* pVBox;
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
