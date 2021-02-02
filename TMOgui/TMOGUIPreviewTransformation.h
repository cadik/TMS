#ifndef TMOGUIPREVIEWTRANSFORMATION_H
#define TMOGUIPREVIEWTRANSFORMATION_H

#include <QObject>
#include <QThread>

class TMOGUIPreviewImage;

class TMOGUIPreviewTransformation : public QThread
{
    Q_OBJECT
public:
    TMOGUIPreviewTransformation(TMOGUIPreviewImage *);
   /* virtual ~TMOGUIPreviewTransformation(void);
protected:
        void run();*/
};

#endif // TMOGUIPREVIEWTRANSFORMATION_H
