#ifndef TMOGUINEWFILE_H
#define TMOGUINEWFILE_H

#include <qvariant.h>


#include <Qt3Support/Q3ButtonGroup>
#include <Qt3Support/Q3GroupBox>
#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>

QT_BEGIN_NAMESPACE

class Ui_TMOGUINewFile
{
public:
    Q3GroupBox *GroupBox1;
    QLabel *TextLabel2;
    QLabel *TextLabel3;
    QLabel *TextLabel4;
    QLineEdit *LineEdit2;
    QLineEdit *LineEdit3;
    QLineEdit *LineEdit1;
    QLabel *TextLabel5;
    QLabel *TextLabel9;
    QLabel *TextLabel10;
    QLabel *TextLabel8;
    QLabel *TextLabel6;
    QLabel *TextLabel13;
    QLabel *TextLabel12;
    QLabel *TextLabel11;
    QLabel *TextLabel14;
    QPushButton *PushButton2;
    QPushButton *PushButton1;
    Q3ButtonGroup *ButtonGroup1;
    QRadioButton *RadioButton1;
    QRadioButton *RadioButton2;
    QLineEdit *LineEdit5;
    QLineEdit *LineEdit6;
    QLineEdit *LineEdit7;
    QLineEdit *LineEdit8;
    QLineEdit *LineEdit9;
    QComboBox *ComboBox1;
    QLineEdit *LineEdit10;

    void setupUi(QDialog *TMOGUINewFile)
    {
        if (TMOGUINewFile->objectName().isEmpty())
            TMOGUINewFile->setObjectName(QString::fromUtf8("TMOGUINewFile"));
        TMOGUINewFile->resize(371, 345);
        GroupBox1 = new Q3GroupBox(TMOGUINewFile);
        GroupBox1->setObjectName(QString::fromUtf8("GroupBox1"));
        GroupBox1->setGeometry(QRect(10, 10, 251, 121));
        TextLabel2 = new QLabel(GroupBox1);
        TextLabel2->setObjectName(QString::fromUtf8("TextLabel2"));
        TextLabel2->setGeometry(QRect(10, 20, 101, 21));
        TextLabel2->setWordWrap(false);
        TextLabel3 = new QLabel(GroupBox1);
        TextLabel3->setObjectName(QString::fromUtf8("TextLabel3"));
        TextLabel3->setGeometry(QRect(10, 60, 100, 21));
        TextLabel3->setWordWrap(false);
        TextLabel4 = new QLabel(GroupBox1);
        TextLabel4->setObjectName(QString::fromUtf8("TextLabel4"));
        TextLabel4->setGeometry(QRect(10, 90, 100, 21));
        TextLabel4->setWordWrap(false);
        LineEdit2 = new QLineEdit(GroupBox1);
        LineEdit2->setObjectName(QString::fromUtf8("LineEdit2"));
        LineEdit2->setGeometry(QRect(110, 60, 131, 22));
        LineEdit3 = new QLineEdit(GroupBox1);
        LineEdit3->setObjectName(QString::fromUtf8("LineEdit3"));
        LineEdit3->setGeometry(QRect(110, 90, 131, 22));
        LineEdit1 = new QLineEdit(GroupBox1);
        LineEdit1->setObjectName(QString::fromUtf8("LineEdit1"));
        LineEdit1->setGeometry(QRect(110, 20, 131, 22));
        TextLabel5 = new QLabel(TMOGUINewFile);
        TextLabel5->setObjectName(QString::fromUtf8("TextLabel5"));
        TextLabel5->setGeometry(QRect(20, 220, 60, 21));
        TextLabel5->setWordWrap(false);
        TextLabel9 = new QLabel(TMOGUINewFile);
        TextLabel9->setObjectName(QString::fromUtf8("TextLabel9"));
        TextLabel9->setGeometry(QRect(20, 280, 60, 21));
        TextLabel9->setWordWrap(false);
        TextLabel10 = new QLabel(TMOGUINewFile);
        TextLabel10->setObjectName(QString::fromUtf8("TextLabel10"));
        TextLabel10->setGeometry(QRect(20, 310, 60, 21));
        TextLabel10->setWordWrap(false);
        TextLabel8 = new QLabel(TMOGUINewFile);
        TextLabel8->setObjectName(QString::fromUtf8("TextLabel8"));
        TextLabel8->setGeometry(QRect(20, 250, 61, 21));
        TextLabel8->setWordWrap(false);
        TextLabel6 = new QLabel(TMOGUINewFile);
        TextLabel6->setObjectName(QString::fromUtf8("TextLabel6"));
        TextLabel6->setGeometry(QRect(210, 220, 70, 21));
        TextLabel6->setWordWrap(false);
        TextLabel13 = new QLabel(TMOGUINewFile);
        TextLabel13->setObjectName(QString::fromUtf8("TextLabel13"));
        TextLabel13->setGeometry(QRect(210, 310, 50, 21));
        TextLabel13->setWordWrap(false);
        TextLabel12 = new QLabel(TMOGUINewFile);
        TextLabel12->setObjectName(QString::fromUtf8("TextLabel12"));
        TextLabel12->setGeometry(QRect(210, 280, 50, 21));
        TextLabel12->setWordWrap(false);
        TextLabel11 = new QLabel(TMOGUINewFile);
        TextLabel11->setObjectName(QString::fromUtf8("TextLabel11"));
        TextLabel11->setGeometry(QRect(210, 250, 50, 21));
        TextLabel11->setWordWrap(false);
        TextLabel14 = new QLabel(TMOGUINewFile);
        TextLabel14->setObjectName(QString::fromUtf8("TextLabel14"));
        TextLabel14->setGeometry(QRect(180, 140, 140, 21));
        TextLabel14->setWordWrap(false);
        PushButton2 = new QPushButton(TMOGUINewFile);
        PushButton2->setObjectName(QString::fromUtf8("PushButton2"));
        PushButton2->setGeometry(QRect(270, 60, 91, 31));
        PushButton1 = new QPushButton(TMOGUINewFile);
        PushButton1->setObjectName(QString::fromUtf8("PushButton1"));
        PushButton1->setGeometry(QRect(270, 20, 91, 31));
        ButtonGroup1 = new Q3ButtonGroup(TMOGUINewFile);
        ButtonGroup1->setObjectName(QString::fromUtf8("ButtonGroup1"));
        ButtonGroup1->setGeometry(QRect(10, 140, 150, 70));
        RadioButton1 = new QRadioButton(ButtonGroup1);
        RadioButton1->setObjectName(QString::fromUtf8("RadioButton1"));
        RadioButton1->setGeometry(QRect(10, 20, 80, 20));
        RadioButton1->setChecked(false);
        RadioButton2 = new QRadioButton(ButtonGroup1);
        RadioButton2->setObjectName(QString::fromUtf8("RadioButton2"));
        RadioButton2->setGeometry(QRect(10, 40, 80, 21));
        RadioButton2->setChecked(true);
        LineEdit5 = new QLineEdit(TMOGUINewFile);
        LineEdit5->setObjectName(QString::fromUtf8("LineEdit5"));
        LineEdit5->setGeometry(QRect(80, 250, 90, 22));
        LineEdit6 = new QLineEdit(TMOGUINewFile);
        LineEdit6->setObjectName(QString::fromUtf8("LineEdit6"));
        LineEdit6->setGeometry(QRect(80, 280, 90, 22));
        LineEdit7 = new QLineEdit(TMOGUINewFile);
        LineEdit7->setObjectName(QString::fromUtf8("LineEdit7"));
        LineEdit7->setGeometry(QRect(80, 310, 90, 22));
        LineEdit8 = new QLineEdit(TMOGUINewFile);
        LineEdit8->setObjectName(QString::fromUtf8("LineEdit8"));
        LineEdit8->setGeometry(QRect(270, 250, 91, 22));
        LineEdit9 = new QLineEdit(TMOGUINewFile);
        LineEdit9->setObjectName(QString::fromUtf8("LineEdit9"));
        LineEdit9->setGeometry(QRect(270, 280, 91, 22));
        ComboBox1 = new QComboBox(TMOGUINewFile);
        ComboBox1->setObjectName(QString::fromUtf8("ComboBox1"));
        ComboBox1->setGeometry(QRect(180, 170, 180, 21));
        LineEdit10 = new QLineEdit(TMOGUINewFile);
        LineEdit10->setObjectName(QString::fromUtf8("LineEdit10"));
        LineEdit10->setGeometry(QRect(270, 310, 91, 22));

        retranslateUi(TMOGUINewFile);
        QObject::connect(PushButton1, SIGNAL(clicked()), TMOGUINewFile, SLOT(accept()));
        QObject::connect(PushButton2, SIGNAL(clicked()), TMOGUINewFile, SLOT(reject()));
        QObject::connect(RadioButton1, SIGNAL(toggled(bool)), ComboBox1, SLOT(setDisabled(bool)));
        QObject::connect(RadioButton1, SIGNAL(toggled(bool)), LineEdit8, SLOT(setDisabled(bool)));
        QObject::connect(RadioButton1, SIGNAL(toggled(bool)), LineEdit9, SLOT(setDisabled(bool)));
        QObject::connect(RadioButton1, SIGNAL(toggled(bool)), LineEdit10, SLOT(setDisabled(bool)));

        QMetaObject::connectSlotsByName(TMOGUINewFile);
    } // setupUi

    void retranslateUi(QDialog *TMOGUINewFile)
    {
        TMOGUINewFile->setWindowTitle(QApplication::translate("TMOGUINewFile", "New File", 0, QApplication::UnicodeUTF8));
        GroupBox1->setTitle(QApplication::translate("TMOGUINewFile", "Properties", 0, QApplication::UnicodeUTF8));
        TextLabel2->setText(QApplication::translate("TMOGUINewFile", "File Name :", 0, QApplication::UnicodeUTF8));
        TextLabel3->setText(QApplication::translate("TMOGUINewFile", "Width :", 0, QApplication::UnicodeUTF8));
        TextLabel4->setText(QApplication::translate("TMOGUINewFile", "Height :", 0, QApplication::UnicodeUTF8));
        TextLabel5->setText(QApplication::translate("TMOGUINewFile", "Minimum :", 0, QApplication::UnicodeUTF8));
        TextLabel9->setText(QApplication::translate("TMOGUINewFile", "Green", 0, QApplication::UnicodeUTF8));
        TextLabel10->setText(QApplication::translate("TMOGUINewFile", "Blue", 0, QApplication::UnicodeUTF8));
        TextLabel8->setText(QApplication::translate("TMOGUINewFile", "Red", 0, QApplication::UnicodeUTF8));
        TextLabel6->setText(QApplication::translate("TMOGUINewFile", "Maximum :", 0, QApplication::UnicodeUTF8));
        TextLabel13->setText(QApplication::translate("TMOGUINewFile", "Blue", 0, QApplication::UnicodeUTF8));
        TextLabel12->setText(QApplication::translate("TMOGUINewFile", "Green", 0, QApplication::UnicodeUTF8));
        TextLabel11->setText(QApplication::translate("TMOGUINewFile", "Red", 0, QApplication::UnicodeUTF8));
        TextLabel14->setText(QApplication::translate("TMOGUINewFile", "Maximum Placement :", 0, QApplication::UnicodeUTF8));
        PushButton2->setText(QApplication::translate("TMOGUINewFile", "Cancel", 0, QApplication::UnicodeUTF8));
        PushButton1->setText(QApplication::translate("TMOGUINewFile", "OK", 0, QApplication::UnicodeUTF8));
        ButtonGroup1->setTitle(QApplication::translate("TMOGUINewFile", "Luminance", 0, QApplication::UnicodeUTF8));
        RadioButton1->setText(QApplication::translate("TMOGUINewFile", "Constant", 0, QApplication::UnicodeUTF8));
        RadioButton2->setText(QApplication::translate("TMOGUINewFile", "Gradient", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class TMOGUINewFile: public Ui_TMOGUINewFile {};
} // namespace Ui

QT_END_NAMESPACE

class TMOGUINewFile : public QDialog, public Ui::TMOGUINewFile
{
    Q_OBJECT

public:
    TMOGUINewFile(QWidget* parent = 0, const char* name = 0, bool modal = false, Qt::WindowFlags fl = 0);
    ~TMOGUINewFile();

protected slots:
    virtual void languageChange();

};

#endif // TMOGUINEWFILE_H
