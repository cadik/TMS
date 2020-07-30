#ifndef TMOGUITONEMAPPINGCHOOSER_H
#define TMOGUITONEMAPPINGCHOOSER_H

#include <QScrollArea>
#include <QGridLayout>
#include <QMap>

class TMO;
class TMOGUIImage;
class TMOGUIBitmap;
class QTemporaryDir;
class TMOGUIToneMapping;

class TMOGUIToneMappingChooser : public QScrollArea
{
    Q_OBJECT
public:
    TMOGUIToneMappingChooser(QWidget *parent = nullptr);
    //! Add preview of TMO
    void AddTMOPreview(TMO* pTmo, int indexLib, int indexTMO);
    virtual ~TMOGUIToneMappingChooser(void);
    //! Display all previews
    void display();
    //! Pointer to parent widget.
    TMOGUIToneMapping* parentWidget;
    //! panel for parameters.
    QWidget* big_box; // TODO Q3VBox
    QVBoxLayout* bigBoxLayout;
    int backWidth;

protected:
    //! Resizes parameters' widgets.
    /*!
    * \param e Resize event.
    */
    virtual void viewportResizeEvent ( QResizeEvent * e );
    int savePreview(TMOGUIImage* pImg);
    int openPreview(TMOGUIImage* pImage);
    TMOGUIImage* pPreviewImage;
    QString* path;
    QMap<TMOGUIImage*, TMO*>::Iterator listPreviewIt;
    QMap<TMOGUIImage*, TMO*> listPreview;
    QSize imageWidgetSize;
    QSize imageFitSize;

public slots:
    void finishTransform();

};

#endif // TMOGUITONEMAPPINGCHOOSER_H
