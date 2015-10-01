// TMOGUIToneMapping.h: interface for the TMOGUIToneMapping class.
//
//////////////////////////////////////////////////////////////////////

#include <qwidget.h>

class QScrollView;
class QComboBox;
class TMOGUIParameters;
class TMO;
class QMultiLineEdit;
class QPushButton;

class TMOGUIToneMapping : public QWidget  
{
	Q_OBJECT
public:
//	int Create(QWidget* pWindow);
	TMOGUIToneMapping( QWidget* parent = 0, const char* name = 0, WFlags f = 0 );
	virtual ~TMOGUIToneMapping();
	TMO** pTMO;
	int iCurTechnique;
	QPushButton *pOk;

protected:
	int FillLibrary();

	QComboBox *pLibrary;
    QComboBox *pTechnique;
	QScrollView *pView;
	QMultiLineEdit* pDescription;
	TMOGUIParameters *pParameters;
	wchar_t** sLibraries;
	int iCurLibrary;
public slots:
	void FillTechnique(int);
	void ChangeTechnique(int);
};


