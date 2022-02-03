// TMOGUIOutput.h: interface for the TMOGUIOutput class.
//
//////////////////////////////////////////////////////////////////////
#include <QTextEdit>
#include <qmap.h>

class TMOImage;
class QString;
class TMOGUIInfo;

class TMOGUIOutput : public QTextEdit
{
	Q_OBJECT
public:
	static int TextLength(const wchar_t *text);

	virtual int Assign(TMOImage *pImage);
	static int WriteLine(TMOImage *pImage, const wchar_t *text);

	TMOGUIOutput(QWidget *parent, const char *name);
	virtual ~TMOGUIOutput();

	static TMOGUIInfo *pInfo;

protected:
	static QMap<TMOImage *, TMOGUIOutput *> mapLocal;
	QString sText;
	unsigned int uiLines;
	bool bVisible;
};
