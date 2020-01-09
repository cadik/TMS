#include "TMOGUISaveDialog.h"

#include <iostream>
#include <qobject.h>
#include <qlineedit.h>



TMOGUISaveDialog::TMOGUISaveDialog( const QString & dirName, filterMap * fm, QWidget * parent, const char * name, bool modal):Q3FileDialog(dirName, "", parent, name, modal)
{
 setMode(Q3FileDialog::AnyFile);
 setSelection(dirName);
 QObjectList  ch = findChildren<QObject*>("QLineEdit");//, "name/filter editor");

 for(QObject* w : ch )
 {
 // std::cout << w->className() << " " << w->name() << std::endl;
  if(w->metaObject()->className() == "QLineEdit") fname=((QLineEdit *)w);
 }
 filterMap::Iterator it;
 for ( it = fm->begin(); it != fm->end(); ++it ) 
 {
  addFilter(it.key());
 }
 filters = fm;
 connect(this, SIGNAL(filterSelected (const QString &)), this, SLOT(change_ext(const QString&)));
}

void TMOGUISaveDialog::change_ext(const QString& filter)
{
 QString fileName = fname->text();
 int iFound=0;
 if ((iFound = fileName.lastIndexOf(".", -1)) > 0)
 {
  fileName = fileName.remove(iFound,fileName.length()-iFound);
 }
 fname->setText(fileName+"."+(*filters)[filter].ext);
}

int TMOGUISaveDialog::getSelectedFileType()
{
 return (*filters)[selectedFilter()].type;
}
