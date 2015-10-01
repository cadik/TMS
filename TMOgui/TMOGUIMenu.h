// TMOGUIMenu.h: interface for the TMOGUIMenu class.
//
//////////////////////////////////////////////////////////////////////
#include <qmenubar.h>
#include <qlist.h>

class QWorkspace;
class TMOGUIImage;

/*! \class TMOGUIMenu	
 *	\brief menu bar of the application
 *
 * Inherits from QMenuBar.
 */
class TMOGUIMenu : public QMenuBar  
{
	Q_OBJECT
public:
	//! Gets recent file from the list according to the id.
	/*!  	
	* \param ID identificator of a recent file.
	* \return Name of the recent file.
	*/
	virtual QString GetRecent(int ID);
	//! Disables a menu-item in the menu bar.
	/*!  	
	* \param menu Ordinal number of menu.
	* \param item Ordinal number of menu item.
	* \return Success indicator.
	*/
	virtual int Disable(int menu, int item);
	//! Enables a menu-item in the menu bar.
	/*!  	
	* \param menu Ordinal number of menu.
	* \param item Ordinal number of menu item.
	* \return Success indicator.
	*/
	virtual int Enable(int menu, int item);
	//! Sets menu according to workspace.
	/*!  
	* Disables or enables menu bar.
	*
	* \param w State of workspace.
	* \return Success indicator.
	*/
	virtual int SetWindows(QWorkspace* w);
	//! Creates menu-items according to list of recent files.
	/*! 
	* \return Success indicator.
	*/
	virtual int SetRecent();	
	//! Adds a recent file to list.
	/*!
	* \param s new file.
	* \return Success indicator.
	*/
	virtual int AddRecent(QString& s);
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
	TMOGUIMenu(QWidget * parent=0, const char * name=0);
	//! Destructor.
	virtual ~TMOGUIMenu();
	//! Public because of the histogram menu item.
	QPopupMenu* pView;

public slots:
	//! Slot for signal when window is changed.
	/*!
	* Checks the active window in the window menu.
	*
	* \param pImg Current active window.
	*/
	void windowChanged(TMOGUIImage* pImg);
	//! Shows dialog containing informations about the application.
	void about();

protected:
	//! Creates menu items and makes appropriate connections.
	virtual int Create();
	//! Loads recent files to the list from file "recent.dat".
	virtual int LoadRecent();
	//! Saves list of recent files to file "recent.dat".
	virtual int SaveRecent();

	QPopupMenu* pFile;
	QPopupMenu* pRecent;	
	QPopupMenu* pEdit;
	QPopupMenu* pCommand;
	QPopupMenu* pWindows;
	QPopupMenu* pHelp;
	QPopupMenu* pComponent;
	QPopupMenu* pHelpIt;
	QWidget* pParent;
	//! List of names of recent files.
	QList<QString> listRecent;
	//! Current active window.
	TMOGUIImage* pImage;
};

