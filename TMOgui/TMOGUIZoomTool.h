// TMOGUIZoomTool.h: interface for the TMOGUIZoomTool class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __TMOGUIZoomTool_H
#define __TMOGUIZoomTool_H

#include <QToolBar>
#include <QSignalMapper>
class TMOGUIImage;
class QMdiArea;
class QToolButton;
class QLineEdit;
class QComboBox;

//! zoom step during zooming
#define ZOOMSTEP 10	

/*! \class TMOGUIZoomTool	
 *	\brief zoom group in toolbar
 *
 * Consists of two buttons, editbox and roll down menu to select desirable zoom value. Inherits from QToolBar.
 */

class TMOGUIZoomTool : public QToolBar
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

	//! Sets the value into the editBox to indicate proper zoom state.
	virtual void SetState();

	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
	TMOGUIZoomTool(QWidget * parent=0, const char * name=0);

	//! Destructor.
	virtual ~TMOGUIZoomTool();

public slots:

	//! Slot for signal when window is changed.
	/*!  
	* \param pWindow Currnet active window.
	*/
	void windowChanged(TMOGUIImage* pWindow);


protected:
	//! Initializes all elements.
	/*!
	* \return Success indicator.
	*/
	virtual int Create();

	//! Zoom out button.
    QToolButton * zoomOutBtn;
    QAction * zoomOutBtnAct;
	//! Editable control to enter value from keyboard.
	QLineEdit * zoomEdit;
    QAction * zoomEditAct;
	//! Zoom in button.
    QToolButton * zoomInBtn;
    QAction * zoomInBtnAct;
	//! Button that invokes roll down menu.
    QToolButton * zoomChoice;
    QAction * zoomChoiceAct;
	//! Pointer to active image window.
	TMOGUIImage * pImage;
	//! Pointer to parent widget.
	QWidget * pParent;

protected slots:
	//! Slot to editbox change.
	void zoom();

	//! Slot to roll down menu item selection.
	/*!  
	* Roll down menuitem's id is used to set zoom value.
	*
	* \param value New zoom value.
	*/
	void zoomWithValue(int value);
private:
    QSignalMapper* signalMapper;

};

#endif //__TMOGUIZoomTool_H
