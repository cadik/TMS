// TMOGUIStatus.h: interface for the TMOGUIStatus class.
//
//////////////////////////////////////////////////////////////////////
#include <qstatusbar.h>
//Added by qt3to4:
#include <QLabel>
//#include "TMOGUIProgressBar.h"

class QLabel;

class TMOGUIStatus : public QStatusBar  
{
	Q_OBJECT
public:
//	virtual int RemovePosition();
//	virtual int AddPosition();
//	virtual int RemoveProgress();
//	virtual int AddProgress();
	TMOGUIStatus( QWidget * parent = 0, const char * name = 0 ) ;
	virtual ~TMOGUIStatus();
	
//	TMOGUIProgressBar* pProgress;
protected:
	virtual int Create();

	QLabel* pPosition;
};

