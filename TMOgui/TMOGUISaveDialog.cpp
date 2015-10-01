#include "TMOGUISaveDialog.h"

//#include <iostream>
#include <qobjectlist.h>
#include <qlineedit.h>



TMOGUISaveDialog::TMOGUISaveDialog( const QString & dirName, filterMap * fm, QWidget * parent, const char * name, bool modal):QFileDialog(dirName, "", parent, name, modal)
{
 setMode(QFileDialog::AnyFile);
 setSelection(dirName);
 QObjectList * ch = queryList("QLineEdit", "name/filter editor"); 
 QObject * w;
 for(w = ch->first(); w; w = ch->next() )
 {
 // std::cout << w->className() << " " << w->name() << std::endl;
  if(w->isA("QLineEdit"))fname=((QLineEdit *)w);
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
 if ((iFound = fileName.findRev(".", -1)) > 0)
 {
  fileName = fileName.remove(iFound,fileName.length()-iFound);
 }
 fname->setText(fileName+"."+(*filters)[filter].ext);
}

int TMOGUISaveDialog::getSelectedFileType()
{
 return (*filters)[selectedFilter()].type;
}
