#ifndef TMOGUITOOL_H
#define TMOGUITOOL_H

#include <qvariant.h>


#include <Qt3Support/Q3ButtonGroup>
#include <Qt3Support/Q3GroupBox>
#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QRadioButton>
#include <QtGui/QScrollBar>

QT_BEGIN_NAMESPACE

class Ui_TMOGUITool
{
public:
    Q3GroupBox *size;
    QLabel *labelMin;
    QLabel *radMax;
    QLineEdit *editSize;
    QScrollBar *scrollBar;
    Q3ButtonGroup *shape;
    QRadioButton *radioCircle;
    QRadioButton *radioSquare;

    void setupUi(QDialog *TMOGUITool)
    {
        if (TMOGUITool->objectName().isEmpty())
            TMOGUITool->setObjectName(QString::fromUtf8("TMOGUITool"));
        TMOGUITool->resize(249, 159);
        size = new Q3GroupBox(TMOGUITool);
        size->setObjectName(QString::fromUtf8("size"));
        size->setGeometry(QRect(10, 10, 230, 80));
        labelMin = new QLabel(size);
        labelMin->setObjectName(QString::fromUtf8("labelMin"));
        labelMin->setGeometry(QRect(20, 20, 20, 20));
        labelMin->setWordWrap(false);
        radMax = new QLabel(size);
        radMax->setObjectName(QString::fromUtf8("radMax"));
        radMax->setGeometry(QRect(140, 20, 30, 20));
        radMax->setWordWrap(false);
        editSize = new QLineEdit(size);
        editSize->setObjectName(QString::fromUtf8("editSize"));
        editSize->setGeometry(QRect(170, 40, 50, 20));
        scrollBar = new QScrollBar(size);
        scrollBar->setObjectName(QString::fromUtf8("scrollBar"));
        scrollBar->setGeometry(QRect(20, 40, 140, 20));
        scrollBar->setMinimum(1);
        scrollBar->setMaximum(150);
        scrollBar->setSingleStep(2);
        scrollBar->setValue(30);
        scrollBar->setOrientation(Qt::Horizontal);
        shape = new Q3ButtonGroup(TMOGUITool);
        shape->setObjectName(QString::fromUtf8("shape"));
        shape->setGeometry(QRect(10, 100, 230, 50));
        radioCircle = new QRadioButton(shape);
        radioCircle->setObjectName(QString::fromUtf8("radioCircle"));
        radioCircle->setGeometry(QRect(40, 20, 50, 20));
        radioCircle->setChecked(true);
        radioSquare = new QRadioButton(shape);
        radioSquare->setObjectName(QString::fromUtf8("radioSquare"));
        radioSquare->setGeometry(QRect(130, 20, 60, 20));

        retranslateUi(TMOGUITool);
        QObject::connect(scrollBar, SIGNAL(valueChanged(int)), TMOGUITool, SLOT(scrollBar_valueChanged(int)));
        QObject::connect(editSize, SIGNAL(textChanged(QString)), TMOGUITool, SLOT(editSize_textChanged(QString)));

        QMetaObject::connectSlotsByName(TMOGUITool);
    } // setupUi

    void retranslateUi(QDialog *TMOGUITool)
    {
        TMOGUITool->setWindowTitle(QApplication::translate("TMOGUITool", "Info Tool", 0, QApplication::UnicodeUTF8));
        size->setTitle(QApplication::translate("TMOGUITool", "Size", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        size->setProperty("toolTip", QVariant(QApplication::translate("TMOGUITool", "Circle - radius, box - square", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
        labelMin->setText(QApplication::translate("TMOGUITool", "1px", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        labelMin->setProperty("toolTip", QVariant(QApplication::translate("TMOGUITool", "Minimum valid value", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
        radMax->setText(QApplication::translate("TMOGUITool", "150px", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        radMax->setProperty("toolTip", QVariant(QApplication::translate("TMOGUITool", "Maximum valid value", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        editSize->setProperty("toolTip", QVariant(QApplication::translate("TMOGUITool", "You can write size from keyboard", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        scrollBar->setProperty("toolTip", QVariant(QApplication::translate("TMOGUITool", "Set Size of tool", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
        shape->setTitle(QApplication::translate("TMOGUITool", "Shape", 0, QApplication::UnicodeUTF8));
        radioCircle->setText(QApplication::translate("TMOGUITool", "circle", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        radioCircle->setProperty("toolTip", QVariant(QApplication::translate("TMOGUITool", "circle shape of tool", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
        radioSquare->setText(QApplication::translate("TMOGUITool", "square", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        radioSquare->setProperty("toolTip", QVariant(QApplication::translate("TMOGUITool", "square shape of tool", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
    } // retranslateUi

};

namespace Ui {
    class TMOGUITool: public Ui_TMOGUITool {};
} // namespace Ui

QT_END_NAMESPACE

class TMOGUITool : public QDialog, public Ui::TMOGUITool
{
    Q_OBJECT

public:
    TMOGUITool(QWidget* parent = 0, const char* name = 0, bool modal = false, Qt::WindowFlags fl = 0);
    ~TMOGUITool();

public slots:
    virtual void scrollBar_valueChanged( int value );
    virtual void editSize_textChanged( const QString & s );

protected slots:
    virtual void languageChange();

};

#endif // TMOGUITOOL_H
