// TMOGUIInfo.h: interface for the TMOGUIInfo class.
//
//////////////////////////////////////////////////////////////////////
#ifndef TMOGUI_INFO_H
#define TMOGUI_INFO_H
#include <qtabwidget.h>
//Added by qt3to4:
#include <QGridLayout>
#include <QLabel>
#include <QResizeEvent>

class TMOGUIOutput;
class TMOGUIStatistics;
class QPushButton;
class QGridLayout;
class QLabel;

/*! \class TMOGUIInfo	
 *	\brief Information tab panel bellow the workspace.
 *
 * Consits of two tabs. The first one is used for displaying image statistics. The second one displays output 
 * informations of actions to be performed.
 */
class TMOGUIInfo : public QTabWidget
{
	Q_OBJECT
public:
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
	TMOGUIInfo(QWidget *parent, const char *name);
	//! Destructor.
	virtual ~TMOGUIInfo();
	//! Pointer to statistics panel.
	TMOGUIStatistics *pStats;

	//! Sets output tab.
	/*!
	* \param pOut Output widget.
	*/
	virtual int SetOutput(QWidget *pOut);

	//! Indicates whether is tab widget visible.
	bool bVisible;

protected:
	//! Moves close button and sends event next.
	/*!
	* \param re Resize event parameter.
	*/
	virtual void resizeEvent(QResizeEvent *e);
	//! Close button in the upper-right corner.
	QPushButton *pClose;
	//! Displays initial message in the output tab.
	QLabel *pLabel;
	//! Current output widget.
	QWidget *pOutput;
public slots:
	//! Forwards a message form close button as "closeBar" signal.
	void hideinfo();

signals:
	//! Indicates close button pressed.
	/*!
	* Used for signal to view menu to know whether the tab panel is active.
	*/
	void closeBar();
};

#endif // TMOGUI_INFO_H
