
#ifndef TMOGUILINERESIZER_H
#define TMOGUILINERESIZER_H

#include <q3frame.h>
//Added by qt3to4:
#include <QMouseEvent>

class QWidget;

/*! \class TMOGUILineResizer	
 *	\brief resizes data rows in information statistcs
 *
 * It's a new element used for resizing rows in info panel. Resizes whole groupBox. Inherits from QFrame.
 */
class TMOGUILineResizer : public Q3Frame
{
    Q_OBJECT

public:
	//! Constructor.
	/*!  
	* Typical construcotor
	*
	* \param parent parent widget, optional, but without correct parent doesn't work properly.
	* \param name Optinonal name
	* \param fl Optional flags
	*/
    TMOGUILineResizer( QWidget* parent = 0, const char* name = 0, Qt::WFlags fl = 0 );
	//! Destructor	
    ~TMOGUILineResizer();
	//! Sets minimal width of parent element
	/*! 
	* Overwritten operator < to compare events according to this key: 
	* \par (position1 < position2 || (position == position2 && type1 < type2)) 
	*
	* \param e Right operand - Event		 
	* \return Returns result of comparing.
	*/
	virtual void SetMinWidth(int min);

protected:
	//! Generates signal resizeInfoElem.
	/*!
	* Generates signal resizeInfoElem only when is left mouse button pressed and mouse's displacement 
	* is in bounderies.
	*
	* \param e Mouse move event parameter.
	*/
    virtual void mouseMoveEvent( QMouseEvent * e );
	//! Sets variable pressed to true.
	/*!	
	* Indicates that left mouse button is pressed.
	*
	* \param e Mouse move event parameter.
	*/
	virtual void mousePressEvent ( QMouseEvent * e );
	//! Sets variable pressed to false.
	/*!	
	* Indicates that left mouse button is released.
	*
	* \param e Mouse move event parameter.
	*/
	virtual void mouseReleaseEvent ( QMouseEvent * e );

	//! Pointer to parent widget.
	QWidget* pParent;
	//! Indicates that left mouse button is pressed.
	bool pressed;
	//! minimum width of the parent widget.
	int minWidth;

signals:
	//! Signals that parent widget has to be resized.
	/*!
	* \param x Relative displacement of parent widget in x-axis.
	*/
	void resizeInfoElem(int x);
};

#endif // TMOGUILINERESIZER_H
