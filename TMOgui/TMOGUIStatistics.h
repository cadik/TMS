// TMOGUIStatistics.h: interface for the TMOGUIStatistics class.
//
//////////////////////////////////////////////////////////////////////
#include <q3scrollview.h>

class TMOGUIImage;
class TMOGUIInfoPanel;

/*! \class TMOGUIStatistics	
 *	\brief Information panel for image statistics.
 *
 * Consists of panel that is placed in viewport of this panel. It enables optional size of panel.
 * Inherits from QScrollView.
 */
class TMOGUIStatistics : public Q3ScrollView  
{
	Q_OBJECT
public:
	//! Constructor.
	/*!  
	* \param parent Optional parameter.
	* \param name Optional parameter.
	*/
	TMOGUIStatistics(QWidget* parent=0, const char * name=0);

	//! Destructor.
	virtual ~TMOGUIStatistics();

protected:
	//! Pointer to active image window.
	TMOGUIImage* pImage;
	//! Panel that is content of this scrollView.
	TMOGUIInfoPanel* pPanel;
	
public slots:
	//! Slot for signal when window is changed.
	/*! 
	* Connects active image to information panel to view his statistics.
	* 
	* \param pWindow Currnet active window.
	*/
	void windowChanged(TMOGUIImage* pWindow);

	//! Slot for signal when window statistics are changed.
	/*! 
	* \param pWindow Currnet active window.
	*/
	void valueschanged();

	//! Accepts Local information tool statistics changes.
	/*! 	
	* \param x x - coordinate in the image.
	* \param y y - coordinate in the image.
	* \param size Size of the shape.
	* \param shape Shape of the tool (square, circle).
	*/
	void printToolStats(int x, int y, int size, bool shape);

	//! Slot for tool cancel signal.
	/*! 	
	*  Clears all local's tool statistics.
	*/
	void resetToolStats();
};
