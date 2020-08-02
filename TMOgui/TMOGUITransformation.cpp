#include "TMOGUITransformation.h"
#include "TMOGUIImage.h"
#include "TMOGUICustomEvents.h"
#include "../tmolib/TMO.h"
#include <qevent.h>
#include <qwaitcondition.h>
#include <qmutex.h>
#include <qapplication.h>
//Added by qt3to4:
#include <QEvent>

QMap<TMOImage*, TMOGUITransformation*> TMOGUITransformation::mapLocal;

TMOGUITransformation::TMOGUITransformation(TMOGUIImage *pImg)
	:pImage(pImg)
{
	TMOImage *pSrc = pImg->GetImage();
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;
	mutex.lock();
	retval = 0;
	for (i = mapLocal.begin(); i != mapLocal.end(); i++)
        if (i.value() == this) // data()
		{
            mapLocal.remove(i.key());
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
    runningMutex.lock();
    stopWaiting = true;
    condition.wakeOne();	// Forcing thread to terminate
    runningMutex.unlock();
	return 0;
}
 
TMOGUITransformation::~TMOGUITransformation(void)
{
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;

	mutex.lock();
	delete refresh;
	refresh = 0;
	for (i = mapLocal.begin(); i != mapLocal.end(); i++) 
        if (i.value() == this)
		{
            mapLocal.remove(i.key());
			break;
		}
	iOperation = -1;
	mutex.unlock();
    runningMutex.lock();
    stopWaiting = true;
    condition.wakeOne();	// Forcing thread to terminate
    runningMutex.unlock();
    QThread::wait();                 // TODO check
}

void TMOGUITransformation::run()
{
	TMOImage *pSrc = 0;
	bool bActive = true;
	
	while(bActive)
	{
        runningMutex.lock();
        stopWaiting = false;
        while(!stopWaiting)
        {
            condition.wait(&runningMutex);
        }
        runningMutex.unlock();

		mutex.lock();
		switch (iOperation)
		{
        case 1:{
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
            TMOGUICustomEvent *ev = new TMOGUICustomEvent((QEvent::User), this );
            QApplication::postEvent( pImage, reinterpret_cast<QEvent*>(ev) );
			mutex.lock();
			iOperation = 0;
			pTMO = 0;
			break;
        }
        case -1:{
			bActive = false;
			break;
        }
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
    pLocal = i.value();
	if (all) *iValue = (part * 100) / all;
	else *iValue = 100;
	
	if (!pLocal) 
	{
		return 0;
	}
    TMOGUICustomEvent *ev = new TMOGUICustomEvent((QEvent::Type)(QEvent::User + 1), (void*)iValue );
    QApplication::postEvent( pLocal->pImage, reinterpret_cast<QEvent*>(ev) );
	pLocal->RefreshGUI();	
	return pLocal->retval;
}

int TMOGUITransformation::WriteLine(TMOImage* pImage, const wchar_t* text)
{
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;
	TMOGUITransformation* pLocal;
    TMOGUICustomEvent* ev;
	
	i = mapLocal.find(pImage);
	if (i == mapLocal.end()) return 0;
    pLocal = i.value(); //data()

	if (!pLocal) 
	{
		return 0;
	}

    ev = new TMOGUICustomEvent((QEvent::Type)(QEvent::User + 2), (void*)text );
    QApplication::postEvent( pLocal->pImage, reinterpret_cast<QEvent*>(ev));
	return 0;
}

int TMOGUITransformation::Assign(TMOImage *pImage)
{
	QMap<TMOImage*, TMOGUITransformation*>::Iterator i;

	for (i = mapLocal.begin(); i != mapLocal.end(); i++)
        if (i.value() == this)
		{
            mapLocal.erase(i);
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
