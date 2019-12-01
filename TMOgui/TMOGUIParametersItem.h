// TMOGUIParametersItem.h: interface for the TMOGUIParametersItem class.
//
//////////////////////////////////////////////////////////////////////
#include <qwidget.h>
//Added by qt3to4:
#include <Q3GridLayout>

class QWidget;
class TMOParameter;
class TMOGUIParameters;
class Q3GridLayout;

/*! \class TMOGUIParametersItem	
 *	\brief zoom group in toolbar
 *
 * Consists of two buttons, editbox and roll down menu to select desirable zoom value. Inherits from QToolBar.
 */
class TMOGUIParametersItem : public QWidget
{
	Q_OBJECT
public:
	//! Creates child widgets according to the type of the parameter.
	/*!
	* Layout is created. There are three type of parameters(TMO_BOOL, TMO_INT, TMO_DOUBLE)
	* \param pParam Parameter to be visualized.
	* \param pParentWidget Parent's widget where to place widgets.
	*/
	int Create(TMOParameter* pParam, TMOGUIParameters* pParentWidget);
	//! Resizes scroll bar's width.
	/*!  
	* \param width Width of the scroll bar.
	*/
	void ResizeElements(int width);
	//! Sets scroll bar according to parameter's value.
	void ResetScrollBar();
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
	TMOGUIParametersItem( QWidget* parent = 0, const char* name = 0 );
	//! Destructor.
	virtual ~TMOGUIParametersItem();
	//! Poublic pointer to the parameter that is displayed.
	TMOParameter* pParameter;

protected:
	//! Flexible amount of widgets according to type of the parameter
	QWidget** pWidgets;
	//! number of widgets
	int iWidgets;
	//! layout for widgets
	Q3GridLayout* pLayout;

public slots:
	//! Deletes child widgets.
	/*!  
	* \param pParentWidget Parent's widget where widgets are placed.
	* \return Success
	*/
	virtual int Destroy(TMOGUIParameters* pParentWidget);
	//! Slot to editBox's change.
	/*!  
	* \param text Value from editBox.
	*/
	void valuechanged(const QString& text);
	//! Slot to checkbox's change. (TMO_BOOL parameter)
	/*!  
	* \param state State of the checkbox.
	*/
	void valuechanged(int state);
	//! Initial setting of parameters.	
	void resetvalues();	
	//! Slot to scroll bar's change.
	/*!  
	* \param value Position in the scroll bar.
	*/
	void scrollbarchanged(int value);
};

