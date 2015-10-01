#include "TMOGUITransformation.h"
#include "TMOGUIImage.h"
#include "../tmolib/TMO.h"
#include <qevent.h>
#include <qwaitcondition.h>
#include <qmutex.h>
#include <qapplication.h>

QMap<TMOImage*, TMOGUITransformation*> TMOGUITransformation::mapLocal;

TMOGUITransformation::TMOGUITransformation(TMOGUIImage *pImg)
	:pImage(pImg)
{
	TMOImage *pSrc = pImg->GetImage();
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;
	mutex.lock();
	retval = 0;
	for (i = mapLocal.begin(); i != mapLocal.end(); i++)
		if (i.data() == this)
		{
			mapLocal.remove(i);
			break;
		}
	pTMO = 0;
	iOperation = 0;
	refresh = new QWaitCondition;
	pSrc->SetProgress(ProgressBar);
	mapLocal.insert(pSrc, this);
	mutex.unlock();
	start(QThread::LowPriority);
}

int TMOGUITransformation::SetTMO(TMO* pToneMap)
{
	mutex.lock();
	if (pTMO) 
	{
		mutex.unlock();
		return -1;
	}
	pTMO = pToneMap;
	iOperation = 1;
	mutex.unlock();	
	condition.wakeOne();
	return 0;
}
 
TMOGUITransformation::~TMOGUITransformation(void)
{
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;

	mutex.lock();
	delete refresh;
	refresh = 0;
	for (i = mapLocal.begin(); i != mapLocal.end(); i++) 
		if (i.data() == this) 
		{
			mapLocal.remove(i);
			break;
		}
	iOperation = -1;
	mutex.unlock();
	condition.wakeOne();	// Forcing thread to terminate
	wait();
}

void TMOGUITransformation::run()
{
	TMOImage *pSrc = 0;
	bool bActive = true;
	
	while(bActive)
	{
		condition.wait();
		mutex.lock();
		switch (iOperation)
		{
		case 1:
			retval = 0;
			pTMO->GetSource(&pSrc);
			pDst = new TMOImage;
			pDst->New(*pSrc);
			pTMO->SetDestination(*pDst);
			mutex.unlock();
			try
			{
				pTMO->Transform();
			}
			catch (int e)
			{
				retval = e;
			}
			QApplication::postEvent( pImage, new QCustomEvent(QEvent::User, this) );
			mutex.lock();
			iOperation = 0;
			pTMO = 0;
			break;
		case -1:
			bActive = false;
			break;
		}
		mutex.unlock();
	}
}

int TMOGUITransformation::ProgressBar(TMOImage* pImage, int part, int all)
{
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;
	TMOGUITransformation* pLocal;
	int* iValue = new(int);
		
	i = mapLocal.find(pImage);
	if (i == mapLocal.end()) return 0;
	pLocal = i.data();
	if (all) *iValue = (part * 100) / all;
	else *iValue = 100;
	
	if (!pLocal) 
	{
		return 0;
	}

	QApplication::postEvent( pLocal->pImage, new QCustomEvent((QEvent::Type)(QEvent::User + 1), (void*)iValue) );
	pLocal->RefreshGUI();	
	return pLocal->retval;
}

int TMOGUITransformation::WriteLine(TMOImage* pImage, const wchar_t* text)
{
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;
	TMOGUITransformation* pLocal;
	
	i = mapLocal.find(pImage);
	if (i == mapLocal.end()) return 0;
	pLocal = i.data();

	if (!pLocal) 
	{
		return 0;
	}
	QApplication::postEvent( pLocal->pImage, new QCustomEvent((QEvent::Type)(QEvent::User + 2), (void*)text ) );
	return 0;
}

int TMOGUITransformation::Assign(TMOImage *pImage)
{
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;

	for (i = mapLocal.begin(); i != mapLocal.end(); i++)
		if (i.data() == this)
		{
			mapLocal.remove(i);
			break;
		}
	mapLocal.insert(pImage, this);
	return 0;
}

TMOImage *TMOGUITransformation::GetDest()
{
	return pDst;
}
int TMOGUITransformation::RefreshGUI()
{
	mutex.lock();
	if (iOperation && refresh)
	{
		mutex.unlock();
		//if (refresh) refresh->wait(10);
		return 0;
	}
	else 
	{
		mutex.unlock();
		qApp->processEvents();
		return 1;
	}
}

int TMOGUITransformation::Cancel()
{
	mutex.lock();
	retval = 1;
	mutex.unlock();
	return 0;
}
