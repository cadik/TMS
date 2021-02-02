// TMOGUIRightBar.h: interface for the TMOGUIRightBar class.
//
//////////////////////////////////////////////////////////////////////
#ifndef TMOGUI_RBAR_H
#define TMOGUI_RBAR_H
#include <qtabwidget.h>
//Added by qt3to4:
#include <QResizeEvent>

class QWidget;
class QPushButton;
class TMOGUIToneMapping;
class TMOGUIFilters;
class TMOGUIStatistics;
class TMO;

class TMOGUIRightBar : public QTabWidget  
{
	Q_OBJECT
public:
    TMOGUIRightBar( QWidget* parent = 0, const char* name = 0/*, Qt::WindowFlags f = 0*/);
	virtual ~TMOGUIRightBar();
	virtual TMO* GetTMO();
	virtual TMOGUIToneMapping* GetMapping() {return pToneMapping;}

protected:

	void resizeEvent ( QResizeEvent * );
    void changeWorkspace(bool advanced);
    QWidget* pVBox;
	QPushButton* pRightButton;
    TMOGUIFilters* pFiltersTab;
	QTabWidget* pTabWidget;
	TMOGUIToneMapping* pToneMapping;
    TMOGUIStatistics* pStats;
    QTabBar* pTabBar;
	bool bVisible;
    bool bAdvanced;

friend class TMOGUIWindow;	

public slots:
	void hideright();
    void paramChanged();

signals:
	void closeBar();
    void change();
};

#endif // TMOGUI_RBAR_H
