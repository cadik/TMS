#include "TMOGUISaveDialog.h"

#include <iostream>
#include <qobject.h>
#include <qlineedit.h>



TMOGUISaveDialog::TMOGUISaveDialog( const QString & dirName, filterMap * fm, QWidget * parent, const char * name, bool modal):QFileDialog(  parent, "", dirName, "")
{
 setFileMode(QFileDialog::FileMode::AnyFile);
 setAcceptMode(QFileDialog::AcceptMode::AcceptSave);
 file = new QString(dirName);
 selectFile(*file);
 /*QObjectList ch = findChildren<QObject*>("QLineEdit");//, "name/filter editor");

 for(QObject* w : ch )
 {
 //std::cout << w->metaObject()->className() << " " << w->objectName() << std::endl;
  if(strcmp(w->metaObject()->className(), "QLineEdit")) fname=((QLineEdit *)w);
 }

 if(fname == nullptr) fname = getSaveFileName();*/

 filterMap::Iterator it;
 QStringList filtersList;
 filtersString = new QString();

 for ( it = fm->begin(); it != fm->end(); ++it ) 
 {
  filtersList.append(it.key());
  if(filtersString->size() > 0){
      *filtersString = filtersString->append(";;");
  }
  *filtersString = filtersString->append(it.key());
 }
 setNameFilters(filtersList);
 filters = fm;


 connect(this, QOverload<const QString &>::of(&QFileDialog::filterSelected), this, QOverload<const QString &>::of(&TMOGUISaveDialog::change_ext));
 connect(this, QOverload<const QString &>::of(&QFileDialog::filterSelected), this, QOverload<const QString &>::of(&QFileDialog::setDefaultSuffix));
 connect(this, QOverload<const QString &>::of(&QFileDialog::fileSelected), this, QOverload<const QString &>::of(&TMOGUISaveDialog::setSelectedFile));

 selectNameFilter(QString("Tiff image HDR(*.tif *.tiff)"));

}

void TMOGUISaveDialog::change_ext(const QString& filter)
{
 QString fileName = *file;
 int iFound=0;
 if ((iFound = fileName.lastIndexOf(".", -1)) > 0)
 {
  fileName = fileName.remove(iFound,fileName.length()-iFound);
 }
 //setNameFilter(fileName+"."+(*filters)[filter].ext);
 selectFile(fileName+"."+(*filters)[filter].ext);
}

int TMOGUISaveDialog::getSelectedFileType()
{
    return (*filters)[selectedNameFilter()].type;
}

void TMOGUISaveDialog::setSelectedFile(const QString& sfile){
    *file = sfile;
}


