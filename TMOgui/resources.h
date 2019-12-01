#ifndef TMOGUIIMAGESIZE_H
#define TMOGUIIMAGESIZE_H

#include <qvariant.h>


#include <Qt3Support/Q3GroupBox>
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

class Ui_TMOGUIImageSize
{
public:
    QCheckBox *CheckBox1;
    Q3GroupBox *GroupBox1;
    QLabel *TextLabel2;
    QLabel *TextLabel1;
    QLineEdit *LineEdit1;
    QLineEdit *LineEdit2;
    QLabel *ratioLabel;
    QLabel *textLabel1_2;
    QLabel *textLabel1;
    QPushButton *PushButton2;
    QPushButton *PushButton1;

    void setupUi(QDialog *TMOGUIImageSize)
    {
        if (TMOGUIImageSize->objectName().isEmpty())
            TMOGUIImageSize->setObjectName(QString::fromUtf8("TMOGUIImageSize"));
        TMOGUIImageSize->resize(285, 151);
        CheckBox1 = new QCheckBox(TMOGUIImageSize);
        CheckBox1->setObjectName(QString::fromUtf8("CheckBox1"));
        CheckBox1->setGeometry(QRect(20, 120, 150, 21));
        GroupBox1 = new Q3GroupBox(TMOGUIImageSize);
        GroupBox1->setObjectName(QString::fromUtf8("GroupBox1"));
        GroupBox1->setGeometry(QRect(10, 10, 170, 100));
        TextLabel2 = new QLabel(GroupBox1);
        TextLabel2->setObjectName(QString::fromUtf8("TextLabel2"));
        TextLabel2->setGeometry(QRect(10, 60, 40, 21));
        TextLabel2->setWordWrap(false);
        TextLabel1 = new QLabel(GroupBox1);
        TextLabel1->setObjectName(QString::fromUtf8("TextLabel1"));
        TextLabel1->setGeometry(QRect(10, 30, 40, 21));
        TextLabel1->setWordWrap(false);
        LineEdit1 = new QLineEdit(GroupBox1);
        LineEdit1->setObjectName(QString::fromUtf8("LineEdit1"));
        LineEdit1->setGeometry(QRect(60, 30, 60, 22));
        LineEdit2 = new QLineEdit(GroupBox1);
        LineEdit2->setObjectName(QString::fromUtf8("LineEdit2"));
        LineEdit2->setGeometry(QRect(60, 60, 60, 22));
        ratioLabel = new QLabel(GroupBox1);
        ratioLabel->setObjectName(QString::fromUtf8("ratioLabel"));
        ratioLabel->setGeometry(QRect(143, 32, 20, 50));
        ratioLabel->setPixmap(qt_get_icon(image0_ID));
        ratioLabel->setWordWrap(false);
        textLabel1_2 = new QLabel(GroupBox1);
        textLabel1_2->setObjectName(QString::fromUtf8("textLabel1_2"));
        textLabel1_2->setGeometry(QRect(126, 60, 20, 20));
        textLabel1_2->setWordWrap(false);
        textLabel1 = new QLabel(GroupBox1);
        textLabel1->setObjectName(QString::fromUtf8("textLabel1"));
        textLabel1->setGeometry(QRect(126, 30, 20, 20));
        textLabel1->setWordWrap(false);
        PushButton2 = new QPushButton(TMOGUIImageSize);
        PushButton2->setObjectName(QString::fromUtf8("PushButton2"));
        PushButton2->setGeometry(QRect(193, 60, 81, 31));
        PushButton1 = new QPushButton(TMOGUIImageSize);
        PushButton1->setObjectName(QString::fromUtf8("PushButton1"));
        PushButton1->setGeometry(QRect(193, 20, 81, 31));

        retranslateUi(TMOGUIImageSize);
        QObject::connect(PushButton1, SIGNAL(clicked()), TMOGUIImageSize, SLOT(accept()));
        QObject::connect(PushButton2, SIGNAL(clicked()), TMOGUIImageSize, SLOT(reject()));
        QObject::connect(CheckBox1, SIGNAL(toggled(bool)), TMOGUIImageSize, SLOT(CheckBox1_toggled(bool)));

        QMetaObject::connectSlotsByName(TMOGUIImageSize);
    } // setupUi

    void retranslateUi(QDialog *TMOGUIImageSize)
    {
        TMOGUIImageSize->setWindowTitle(QApplication::translate("TMOGUIImageSize", "Image Size", 0, QApplication::UnicodeUTF8));
        CheckBox1->setText(QApplication::translate("TMOGUIImageSize", "Constrain proportions ", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        CheckBox1->setProperty("toolTip", QVariant(QApplication::translate("TMOGUIImageSize", "Constrain aspect ratio", 0, QApplication::UnicodeUTF8)));
#endif // QT_NO_TOOLTIP
        GroupBox1->setTitle(QApplication::translate("TMOGUIImageSize", "Dimensions", 0, QApplication::UnicodeUTF8));
        TextLabel2->setText(QApplication::translate("TMOGUIImageSize", "Height :", 0, QApplication::UnicodeUTF8));
        TextLabel1->setText(QApplication::translate("TMOGUIImageSize", "Width :", 0, QApplication::UnicodeUTF8));
        ratioLabel->setText(QString());
        textLabel1_2->setText(QApplication::translate("TMOGUIImageSize", "px", 0, QApplication::UnicodeUTF8));
        textLabel1->setText(QApplication::translate("TMOGUIImageSize", "px", 0, QApplication::UnicodeUTF8));
        PushButton2->setText(QApplication::translate("TMOGUIImageSize", "Cancel", 0, QApplication::UnicodeUTF8));
        PushButton1->setText(QApplication::translate("TMOGUIImageSize", "OK", 0, QApplication::UnicodeUTF8));
    } // retranslateUi


protected:
    enum IconID
    {
        image0_ID,
        unknown_ID
    };
    static QPixmap qt_get_icon(IconID id)
    {
    static const unsigned char image0_data[] = { 
    0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
    0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x14, 0x00, 0x00, 0x00, 0x32,
    0x08, 0x06, 0x00, 0x00, 0x00, 0x5c, 0x7c, 0xfb, 0x26, 0x00, 0x00, 0x00,
    0x4a, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0xed, 0xd7, 0xb1, 0x0e, 0x00,
    0x10, 0x0c, 0x45, 0xd1, 0x12, 0xff, 0xff, 0xcb, 0x6c, 0x22, 0x2f, 0x16,
    0xd5, 0xa1, 0xe2, 0xde, 0xad, 0xcb, 0x49, 0x6b, 0x63, 0xf6, 0x5d, 0x45,
    0xe6, 0x7e, 0xeb, 0x54, 0xff, 0x2e, 0xfb, 0xc2, 0x41, 0x3d, 0xf9, 0xa4,
    0xf5, 0x79, 0x1e, 0x3a, 0x19, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10,
    0x10, 0x10, 0x10, 0x30, 0xb6, 0x26, 0xb3, 0xf7, 0x27, 0x35, 0xcb, 0x7f,
    0x72, 0xfe, 0x06, 0xc2, 0x32, 0x03, 0x5a, 0xb3, 0xf9, 0xa1, 0xc3, 0x00,
    0x00, 0x00, 0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
};

    switch (id) {
        case image0_ID:  { QImage img; img.loadFromData(image0_data, sizeof(image0_data), "PNG"); return QPixmap::fromImage(img); }
        default: return QPixmap();
    } // switch
    } // icon

};

namespace Ui {
    class TMOGUIImageSize: public Ui_TMOGUIImageSize {};
} // namespace Ui

QT_END_NAMESPACE

class TMOGUIImageSize : public QDialog, public Ui::TMOGUIImageSize
{
    Q_OBJECT

public:
    TMOGUIImageSize(QWidget* parent = 0, const char* name = 0, bool modal = false, Qt::WindowFlags fl = 0);
    ~TMOGUIImageSize();

public slots:
    virtual void CheckBox1_toggled( bool on );

protected slots:
    virtual void languageChange();

};

#endif // TMOGUIIMAGESIZE_H
