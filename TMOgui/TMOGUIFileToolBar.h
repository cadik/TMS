// TMOGUIFileToolBar.h: interface for the TMOGUIFileToolBar class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __TMOGUIFileToolBar_H
#define __TMOGUIFileToolBar_H

#include <QToolBar>
class QMdiArea;
class QToolButton;
class Switch;

/*! \class TMOGUIFileToolBar	
 *	\brief file group in toolbar
 *
 * Group of tools for file manipulations. Inherits from QToolBar.
 */
class TMOGUIFileToolBar : public QToolBar
{
	Q_OBJECT	
public:

	//! Sets element according to workspace.
	/*!  
	* Disables or enables group in the toolbar.
	*
	* \param w State of workspace.
	* \return Success indicator.
	*/

    virtual int SetWindows(QMdiArea* w);
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/

	TMOGUIFileToolBar(QWidget * parent=0, const char * name=0);

	//! Destructor
	virtual ~TMOGUIFileToolBar();

    Switch* switchWorkspace;


protected:
	//! Initializes all elements.
	/*!
	* \return Success indicator.
	*/
	virtual int Create();

	//! new file creation tool button.
    QAction * fileNewBtn;
	//! file open tool button.
    QAction * fileOpenBtn;
	//! file save tool button.
    QAction * fileSaveBtn;
	//! file save as tool button.
    QAction * fileSaveAllBtn;
	//! file print tool button.
    QAction * filePrintBtn;
    //! toggle switch for workspace
    QAction * fileWorkspace;
	//! Pointer to parent widget.
	QWidget * pParent;

private:
    QAction *addButton(const QIcon& s, const QString &textLabel,
                          const QString& grouptext,
                          QObject * receiver, const char* slot,
                          QWidget * parent, QAction *before);
};

#endif //__TMOGUIFileToolBar_H
