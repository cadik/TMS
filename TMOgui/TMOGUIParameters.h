// TMOGUIParameters.h: interface for the TMOGUIParameters class.
//
//////////////////////////////////////////////////////////////////////

#include <qscrollview.h>
#include <string>

class TMOGUIParametersItem;
class QGridLayout;
class TMO;
class QVBox;
class QString;

/*! \class TMOGUIParameters	
 *	\brief place for tone mapping method's parameters.
 *
 * There are placed instances of TMOParameterItem. Inherits from QScrollView.
 */
class TMOGUIParameters : public QScrollView
{
	Q_OBJECT
public:
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
	TMOGUIParameters(QWidget * parent=0, const char * name=0);
	//! Destructor.	
	virtual ~TMOGUIParameters();
	//! Sets technique and creates his parameters.
	/*!  
	* \param PTmo.
	* \return Success.
	*/
	virtual int SetTechnique(TMO* pTmo);
	static QString GetString(const wchar_t* s);

	//! Pointer to parent widget.
	QWidget* parentWidget;
	//! panel for parameters.
	QVBox* big_box;
	//! Used for parameterItems. Indicates last created parameter;
	int iCurParam;
	int backWidth;

protected:
	//! Number of parameters to be created.
	int iParams;
	//! Widgets for parameters
	TMOGUIParametersItem** pParams;

	//! Resizes parameters' widgets.
	/*!  
	* \param e Resize event.
	*/
	virtual void viewportResizeEvent ( QResizeEvent * e );

public slots:
	//! Resets all parameters to default values.
	void resetvalues();
};
