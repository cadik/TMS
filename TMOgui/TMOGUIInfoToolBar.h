// TMOGUIZoomTool.h: interface for the TMOGUIZoomTool class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __TMOGUIINFOTOOLBAR_H
#define __TMOGUIINFOTOOLBAR_H

#include <q3toolbar.h>

class QWorkspace;
class QToolButton;

/*! \class TMOGUIInfoToolBar	
 *	\brief local information tool group in toolbar
 *
 * Consists of two buttons. The first one activates the tool, the second shows tool's setting panel. 
 * Inherits from QToolBar.
 */
class TMOGUIInfoToolBar : public Q3ToolBar
{
	Q_OBJECT
public:
	//! Return whether the tool is activated.
	/*!
	* \return Bool value if tool is activated.
	*/
	virtual bool IsActivated() const;
	//! Sets elements according to workspace.
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
	TMOGUIInfoToolBar(QWidget * parent=0, const char * name=0);
	//! Destructor
	virtual ~TMOGUIInfoToolBar();	

protected:
	//! Creates buttons and makes appropriate connections.
	virtual int Create();

	//! Tool button that activates the tool.
	QToolButton * toolBtn;
	//! Tool button that shows setting panel of the tool.
	QToolButton * toolSettingBtn;
	//! Pointer to parent widget.
	QWidget * pParent;	

};

#endif //__TMOGUIINFOTOOLBAR_H