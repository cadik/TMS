// TMOGUIProgressBar.h: interface for the TMOGUIProgressBar class.
//
//////////////////////////////////////////////////////////////////////
#include <qhbox.h>
#include <qmap.h>

class QProgressBar;
class QPushButton;
class QLabel;
class QStatusBar;
class TMOImage;

class TMOGUIProgressBar : public QHBox  
{
	Q_OBJECT
public:
	int SetProgress(int part, int all);
	//virtual int Assign(TMOImage* pImage);
	virtual int SetLabel(QString text);
	virtual bool isVisible();
	TMOGUIProgressBar(QStatusBar * parent=0, const char * name=0);
	virtual ~TMOGUIProgressBar();
	//static int ProgressBar(TMOImage* pImage, int part, int all);
protected:
	bool bCancel;
	bool bVisible;
	int iLast;
	QProgressBar *pProgress;
	QPushButton *pButton;
	QLabel *pLabel;
	QStatusBar *pParent;
	//static QMap<TMOImage*, TMOGUIProgressBar*> mapLocal;
public slots:
	void cancel();
signals:
	void cancelsignal();
};
