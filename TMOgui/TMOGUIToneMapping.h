// TMOGUIToneMapping.h: interface for the TMOGUIToneMapping class.
//
//////////////////////////////////////////////////////////////////////

#include <qwidget.h>
#include "TMOGUIImage.h"

class QScrollArea;
class QComboBox;
class QGroupBox;
class QGridLayout;
class TMOGUIToneMappingChooser;
class TMOGUIParameters;
class TMOGUIBitmap;
class TMOImage;
class TMO;
class QTextEdit;
class QPushButton;

class TMOGUIToneMapping : public QWidget
{
    Q_OBJECT
public:
    //	int Create(QWidget* pWindow);
    TMOGUIToneMapping(QWidget *parent = 0, const char *name = 0 /*, Qt::WindowFlags f = 0*/);
    virtual ~TMOGUIToneMapping();
    void changeWorkspace(bool advanced, TMOGUIImage *pImg);
    void displayPreviews(TMOGUIImage *pImg);
    TMO **pTMO;
    int iCurTechnique;
    QPushButton *pOk;
    QPushButton *pPreview;
    QScrollArea *pScrollView;
    TMOGUIToneMappingChooser *pChooser;

protected:
    int FillLibrary();

    QWidget *pParametersParent;
    QComboBox *pLibrary;
    QComboBox *pTechnique;
    QGroupBox *pGroupBox;
    QScrollArea *pView;
    QTextEdit *pDescription;
    TMOGUIParameters *pParameters;
    QPushButton *pBackButton;
    QGridLayout *gridLayout;
    QVBoxLayout *vBoxLayout;
    QMap<QString, TMOGUIToneMappingChooser *> mChoosers;
    bool bIsSwitching = false;
    bool bAdvanced;

    wchar_t **sLibraries;
    int iCurLibrary;
public slots:
    void FillTechnique(int);
    void ChangeTechnique(int);
    void getAllTechniques();
    void paramChanged();
    void toggleTechniqueChooser();
    void chooseTechnique(int indexLib, int indexTMO);
    void windowChanged(TMOGUIImage *pImage);

signals:
    void change();
};
