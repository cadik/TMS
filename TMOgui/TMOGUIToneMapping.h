// TMOGUIToneMapping.h: interface for the TMOGUIToneMapping class.
//
//////////////////////////////////////////////////////////////////////

#include <qwidget.h>

class Q3ScrollView;
class QComboBox;
class TMOGUIParameters;
class TMO;
class Q3MultiLineEdit;
class QPushButton;

class TMOGUIToneMapping : public QWidget  
{
	Q_OBJECT
public:
//	int Create(QWidget* pWindow);
	TMOGUIToneMapping( QWidget* parent = 0, const char* name = 0, Qt::WFlags f = 0 );
	virtual ~TMOGUIToneMapping();
	TMO** pTMO;
	int iCurTechnique;
	QPushButton *pOk;

protected:
	int FillLibrary();

	QComboBox *pLibrary;
    QComboBox *pTechnique;
	Q3ScrollView *pView;
	Q3MultiLineEdit* pDescription;
	TMOGUIParameters *pParameters;
	wchar_t** sLibraries;
	int iCurLibrary;
public slots:
	void FillTechnique(int);
	void ChangeTechnique(int);
};


