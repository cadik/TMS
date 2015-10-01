#ifndef __INFOTOOL_H
#define __INFOTOOL_H

#include <qobject.h>

class TMOGUIBitmap;
class TMOGUITool;

/*! \class TMOGUIInfoTool	
 *	\brief zoom group in toolbar
 *
 * Consists of two buttons, editbox and roll down menu to select desirable zoom value. Inherits from QToolBar.
 */
class TMOGUIInfoTool : public QObject
{
	Q_OBJECT

public:	
	//! Draws the tool's shape on the image.
	/*!
	* It is drawn accordint ot mouse position.
	*
	* \param paint QPainter parameter.
	*/
	virtual void DrawTool(QPainter & paint);
	//! Changes when mouse move over the image is performed.
	/*!
	* Sets the current mouse position. Emits signal to statistics and updates the current image.
	*
	* \param e Mouse move event.
	*/
	virtual void MouseAction(QMouseEvent * e);
	//! Creates context dialog to set the tool's parameters.
	virtual void CreateContextMenu();
	//! Disables tool.
	/*!
	* Emits signal "toolCancelled". Disables tool.
	*/
	virtual void LeaveImage();
	//! Returns whether tool is active.
	/*!
	* \return Bool value whether tool is active.
	*/
	virtual bool IsEnabled() const { return enableTool; }
	//! Activates the tool.
	/*!
	* \param on Bool value.
	*/
	virtual void SetEnabled(bool on);
	//! Sets tool's parameters to default values.
	virtual void Reset();
	//! Sets current image to apply the tool.
	/*!
	* \param parentBp New selected window.
	*/
	virtual void setToolPtr(TMOGUIBitmap * parentBp) { parentBitmap = parentBp; }
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
	TMOGUIInfoTool(QObject * parent = 0, const char * name = 0);
	//! Destructor.
	virtual ~TMOGUIInfoTool();
	//! X coordinate of mouse.
	unsigned int x;
	//! Y coordinate of mouse.
	unsigned int y; 
	//! Size of the shape to be drawn.
	unsigned int toolSize;
	//! Tool setting dialog pointer.
	TMOGUITool * toolContext;

protected:
	//! Manages if the tool's shape have to be drawn.
	/*!
	* It used when mouse leaves the current image. Otherwise tool is drawn on the inactive windows.
	*/
	bool drawToolIcon;
	//! Indicates if tool is active.
	bool enableTool;
	//! Determines tool's shape { true - circle, false - squere }
	bool circleShape;
	//! It enables drawing when context dialog is visible to better see parameters results.
	bool contextDialogEnable;	
	//! Pointer to active image window.
	TMOGUIBitmap * parentBitmap;

public slots:
	//! Size of tool is changed.
	/*!
	* \param s strig value of edit box (TMOGUITool).
	*/
	void changeTextSize(const QString & s);
	//! Sets size of tool according to scrollbar (TMOGUITool).
	/*!
	* \param value Position in scrollbar.
	*/
	void changeSize(int value);
	//! Circle shape is selected in radio buttons (TMOGUITool).
	/*!
	* \param on Bool value.
	*/
	void changeCircle(bool on);
	//! Square shape is selected in radio buttons (TMOGUITool).
	/*!
	* \param on Bool value.
	*/
	void changeSquare(bool on);

signals:
	//! Sends data to statistics information panel.
	/*!
	* \param x Position in the image.
	* \param y Position in the image.
	* \param size Size of the tool.
	* \param shape Shape of the tool { circle, square }.
	*/
	void toolApllied(int x, int y, int size, bool shape);
	//! Tools leaves the active window.
	void toolCancelled();
};


#endif	//__INFOTOOL_H