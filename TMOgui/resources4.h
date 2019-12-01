#ifndef TMOGUIPAGESETUP_H
#define TMOGUIPAGESETUP_H

#include <qvariant.h>


#include <Qt3Support/Q3GroupBox>
#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_TMOGUIPageSetup
{
public:
    QPushButton *PushButton1;
    QPushButton *PushButton2;
    QCheckBox *CheckBox1;
    Q3GroupBox *GroupBox1;
    QLabel *TextLabel3;
    QLabel *TextLabel4;
    QLabel *TextLabel2;
    QLabel *TextLabel5;
    QLineEdit *LineEdit1;
    QLineEdit *LineEdit3;
    QLineEdit *LineEdit4;
    QLineEdit *LineEdit2;

    void setupUi(QDialog *TMOGUIPageSetup)
    {
        if (TMOGUIPageSetup->objectName().isEmpty())
            TMOGUIPageSetup->setObjectName(QString::fromUtf8("TMOGUIPageSetup"));
        TMOGUIPageSetup->resize(385, 200);
        PushButton1 = new QPushButton(TMOGUIPageSetup);
        PushButton1->setObjectName(QString::fromUtf8("PushButton1"));
        PushButton1->setGeometry(QRect(270, 20, 101, 31));
        PushButton2 = new QPushButton(TMOGUIPageSetup);
        PushButton2->setObjectName(QString::fromUtf8("PushButton2"));
        PushButton2->setGeometry(QRect(270, 60, 101, 31));
        CheckBox1 = new QCheckBox(TMOGUIPageSetup);
        CheckBox1->setObjectName(QString::fromUtf8("CheckBox1"));
        CheckBox1->setGeometry(QRect(270, 120, 101, 21));
        GroupBox1 = new Q3GroupBox(TMOGUIPageSetup);
        GroupBox1->setObjectName(QString::fromUtf8("GroupBox1"));
        GroupBox1->setGeometry(QRect(10, 10, 250, 180));
        TextLabel3 = new QLabel(GroupBox1);
        TextLabel3->setObjectName(QString::fromUtf8("TextLabel3"));
        TextLabel3->setGeometry(QRect(10, 60, 30, 21));
        TextLabel3->setWordWrap(false);
        TextLabel4 = new QLabel(GroupBox1);
        TextLabel4->setObjectName(QString::fromUtf8("TextLabel4"));
        TextLabel4->setGeometry(QRect(200, 60, 40, 21));
        TextLabel4->setWordWrap(false);
        TextLabel2 = new QLabel(GroupBox1);
        TextLabel2->setObjectName(QString::fromUtf8("TextLabel2"));
        TextLabel2->setGeometry(QRect(110, 20, 30, 21));
        TextLabel2->setWordWrap(false);
        TextLabel5 = new QLabel(GroupBox1);
        TextLabel5->setObjectName(QString::fromUtf8("TextLabel5"));
        TextLabel5->setGeometry(QRect(100, 140, 50, 21));
        TextLabel5->setWordWrap(false);
        LineEdit1 = new QLineEdit(GroupBox1);
        LineEdit1->setObjectName(QString::fromUtf8("LineEdit1"));
        LineEdit1->setGeometry(QRect(70, 40, 101, 22));
        LineEdit1->setAlignment(Qt::AlignRight);
        LineEdit3 = new QLineEdit(GroupBox1);
        LineEdit3->setObjectName(QString::fromUtf8("LineEdit3"));
        LineEdit3->setGeometry(QRect(130, 80, 101, 22));
        LineEdit3->setAlignment(Qt::AlignRight);
        LineEdit4 = new QLineEdit(GroupBox1);
        LineEdit4->setObjectName(QString::fromUtf8("LineEdit4"));
        LineEdit4->setGeometry(QRect(70, 120, 100, 22));
        LineEdit4->setAlignment(Qt::AlignRight);
        LineEdit2 = new QLineEdit(GroupBox1);
        LineEdit2->setObjectName(QString::fromUtf8("LineEdit2"));
        LineEdit2->setGeometry(QRect(10, 80, 101, 22));
        LineEdit2->setAlignment(Qt::AlignRight);

        retranslateUi(TMOGUIPageSetup);
        QObject::connect(PushButton1, SIGNAL(clicked()), TMOGUIPageSetup, SLOT(accept()));
        QObject::connect(PushButton2, SIGNAL(clicked()), TMOGUIPageSetup, SLOT(reject()));

        QMetaObject::connectSlotsByName(TMOGUIPageSetup);
    } // setupUi

    void retranslateUi(QDialog *TMOGUIPageSetup)
    {
        TMOGUIPageSetup->setWindowTitle(QApplication::translate("TMOGUIPageSetup", "Page Setup", 0, QApplication::UnicodeUTF8));
        PushButton1->setText(QApplication::translate("TMOGUIPageSetup", "OK", 0, QApplication::UnicodeUTF8));
        PushButton2->setText(QApplication::translate("TMOGUIPageSetup", "Cancel", 0, QApplication::UnicodeUTF8));
        CheckBox1->setText(QApplication::translate("TMOGUIPageSetup", "Zoom To Fit", 0, QApplication::UnicodeUTF8));
        GroupBox1->setTitle(QApplication::translate("TMOGUIPageSetup", "Margins [%]", 0, QApplication::UnicodeUTF8));
        TextLabel3->setText(QApplication::translate("TMOGUIPageSetup", "Left", 0, QApplication::UnicodeUTF8));
        TextLabel4->setText(QApplication::translate("TMOGUIPageSetup", "Right", 0, QApplication::UnicodeUTF8));
        TextLabel2->setText(QApplication::translate("TMOGUIPageSetup", "Top", 0, QApplication::UnicodeUTF8));
        TextLabel5->setText(QApplication::translate("TMOGUIPageSetup", "Bottom", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class TMOGUIPageSetup: public Ui_TMOGUIPageSetup {};
} // namespace Ui

QT_END_NAMESPACE

class TMOGUIPageSetup : public QDialog, public Ui::TMOGUIPageSetup
{
    Q_OBJECT

public:
    TMOGUIPageSetup(QWidget* parent = 0, const char* name = 0, bool modal = false, Qt::WindowFlags fl = 0);
    ~TMOGUIPageSetup();

protected slots:
    virtual void languageChange();

};

#endif // TMOGUIPAGESETUP_H
