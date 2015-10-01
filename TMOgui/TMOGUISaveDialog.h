#ifndef TMOGUISAVEDIALOG
#define TMOGUISAVEDIALOG

#include <qfiledialog.h>
#include <qmap.h>
#include <qstring.h>

struct fileType
{
 int type;
 QString ext;
 fileType(int _type, const QString& _ext):type(_type), ext(_ext){}
 fileType(){}
};

typedef QMap<QString, fileType> filterMap;

class TMOGUISaveDialog : public QFileDialog
{
 Q_OBJECT
public:
 private:
  QLineEdit * fname;
  filterMap * filters;
 protected:
 protected slots:
  void change_ext(const QString& filter);
 public:
  TMOGUISaveDialog( const QString & dirName, filterMap * fm, QWidget * parent = 0, const char * name = 0, bool modal = FALSE );
  int getSelectedFileType();

};



#endif
