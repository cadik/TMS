// TMOGUIFileToolBar.h: interface for the TMOGUIFileToolBar class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __TMOGUIFileToolBar_H
#define __TMOGUIFileToolBar_H

#include <q3toolbar.h>
class QWorkspace;
class QToolButton;

/*! \class TMOGUIFileToolBar	
 *	\brief file group in toolbar
 *
 * Group of tools for file manipulations. Inherits from QToolBar.
 */
class TMOGUIFileToolBar : public Q3ToolBar
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

	virtual int SetWindows(QWorkspace* w);
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/

	TMOGUIFileToolBar(QWidget * parent=0, const char * name=0);

	//! Destructor
	virtual ~TMOGUIFileToolBar();


protected:
	//! Initializes all elements.
	/*!
	* \return Success indicator.
	*/
	virtual int Create();

	//! new file creation tool button.
	QToolButton * fileNewBtn;
	//! file open tool button.
	QToolButton * fileOpenBtn;
	//! file save tool button.
	QToolButton * fileSaveBtn;	
	//! file save as tool button.
	QToolButton * fileSaveAllBtn;	
	//! file print tool button.
	QToolButton * filePrintBtn;
	//! Pointer to parent widget.
	QWidget * pParent;
};

#endif //__TMOGUIFileToolBar_H