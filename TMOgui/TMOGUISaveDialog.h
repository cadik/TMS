#ifndef TMOGUISAVEDIALOG
#define TMOGUISAVEDIALOG

#include <QFileDialog>
#include <QLineEdit>
#include <qmap.h>
#include <qstring.h>

struct fileType
{
  int type;
  QString ext;
  fileType(int _type, const QString &_ext) : type(_type), ext(_ext) {}
  fileType() {}
};

typedef QMap<QString, fileType> filterMap;

class TMOGUISaveDialog : public QFileDialog
{
  Q_OBJECT
public:
  filterMap *filters;
  QString *filtersString;
  QString *file;

private:
  QString fname;

protected:
protected slots:
  void change_ext(const QString &filter);
  void setSelectedFile(const QString &file);

public:
  TMOGUISaveDialog(const QString &dirName, filterMap *fm, QWidget *parent = 0, const char *name = 0, bool modal = false);
  int getSelectedFileType();
};

#endif
