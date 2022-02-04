// TMOGUIStatus.h: interface for the TMOGUIStatus class.
//
//////////////////////////////////////////////////////////////////////
#ifndef TMOGUI_STATUS_H
#define TMOGUI_STATUS_H
#include <qstatusbar.h>
//Added by qt3to4:
#include <QLabel>
#include "TMOGUIProgressBar.h"

class QLabel;

class TMOGUIStatus : public QStatusBar
{
	Q_OBJECT
public:
	//	virtual int RemovePosition();
	//	virtual int AddPosition();
	//	virtual int RemoveProgress();
	//	virtual int AddProgress();
	TMOGUIStatus(QWidget *parent = nullptr, const char *name = nullptr);
	virtual ~TMOGUIStatus();

	TMOGUIProgressBar *pProgress;

protected:
	virtual int Create();

	QLabel *pPosition;
};

#endif // TMOGUI_STATUS_H
