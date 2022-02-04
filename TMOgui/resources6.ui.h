/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you want to add, delete, or rename functions or slots, use
** Qt Designer to update this file, preserving your code.
**
** You should not define a constructor or destructor in this file.
** Instead, write your code in functions called init() and destroy().
** These will automatically be called by the form's constructor and
** destructor.
*****************************************************************************/

void TMOGUITool::scrollBar_valueChanged(int value)
{
    QString s;
    editSize->setText(s.setNum(value));
}

void TMOGUITool::editSize_textChanged(const QString &s)
{
    bool ok;
    int size = s.toInt(&ok);
    if (!ok || 1 > size || size > 150)
        return;
    scrollBar->setValue(size);
}
