// TMOGUIToneMapping.h: interface for the TMOGUIToneMapping class.
//
//////////////////////////////////////////////////////////////////////

#include <qwidget.h>

class QScrollArea;
class QComboBox;
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
    TMOGUIToneMapping( QWidget* parent = 0, const char* name = 0, Qt::WindowFlags f = 0 );
	virtual ~TMOGUIToneMapping();
	TMO** pTMO;
	int iCurTechnique;
	QPushButton *pOk;
    QPushButton *pPreview;
    QScrollArea *pScrollView;

protected:
	int FillLibrary();

	QComboBox *pLibrary;
    QComboBox *pTechnique;
    QScrollArea *pView;
    QTextEdit* pDescription;
	TMOGUIParameters *pParameters;
	wchar_t** sLibraries;
	int iCurLibrary;
public slots:
	void FillTechnique(int);
	void ChangeTechnique(int);
    void paramChanged();
signals:
    void change();
};


