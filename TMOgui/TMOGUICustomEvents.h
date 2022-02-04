#ifndef TMOGUICUSTOMEVENTS_H
#define TMOGUICUSTOMEVENTS_H

#include <QEvent>

class TMOGUICustomEvent : QEvent
{
public:
    TMOGUICustomEvent(QEvent::Type type, void *event_data) : QEvent(type)
    {
        eventdata = event_data;
    }

    void *data()
    {
        return eventdata;
    }

protected:
    void *eventdata;
};

#endif // TMOGUICUSTOMEVENTS_H
