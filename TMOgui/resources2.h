#ifndef TMOGUIOPERATION_H
#define TMOGUIOPERATION_H

#include <qvariant.h>


#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_TMOGUIOperation
{
public:
    QPushButton *PushButton3;
    QPushButton *PushButton4;
    QLabel *TextLabel1;
    QComboBox *ComboBox1;
    QComboBox *ComboBox2;
    QComboBox *ComboBox3;
    QLabel *PixmapLabel1;
    QLabel *PixmapLabel2;

    void setupUi(QDialog *TMOGUIOperation)
    {
        if (TMOGUIOperation->objectName().isEmpty())
            TMOGUIOperation->setObjectName(QString::fromUtf8("TMOGUIOperation"));
        TMOGUIOperation->resize(489, 282);
        PushButton3 = new QPushButton(TMOGUIOperation);
        PushButton3->setObjectName(QString::fromUtf8("PushButton3"));
        PushButton3->setGeometry(QRect(390, 10, 91, 31));
        PushButton4 = new QPushButton(TMOGUIOperation);
        PushButton4->setObjectName(QString::fromUtf8("PushButton4"));
        PushButton4->setGeometry(QRect(390, 50, 91, 31));
        TextLabel1 = new QLabel(TMOGUIOperation);
        TextLabel1->setObjectName(QString::fromUtf8("TextLabel1"));
        TextLabel1->setGeometry(QRect(10, 230, 101, 21));
        TextLabel1->setWordWrap(false);
        ComboBox1 = new QComboBox(TMOGUIOperation);
        ComboBox1->setObjectName(QString::fromUtf8("ComboBox1"));
        ComboBox1->setGeometry(QRect(10, 200, 181, 21));
        ComboBox2 = new QComboBox(TMOGUIOperation);
        ComboBox2->setObjectName(QString::fromUtf8("ComboBox2"));
        ComboBox2->setGeometry(QRect(200, 200, 181, 21));
        ComboBox3 = new QComboBox(TMOGUIOperation);
        ComboBox3->setObjectName(QString::fromUtf8("ComboBox3"));
        ComboBox3->setGeometry(QRect(10, 250, 181, 21));
        PixmapLabel1 = new QLabel(TMOGUIOperation);
        PixmapLabel1->setObjectName(QString::fromUtf8("PixmapLabel1"));
        PixmapLabel1->setGeometry(QRect(20, 20, 161, 161));
        PixmapLabel1->setScaledContents(false);
        PixmapLabel1->setAlignment(Qt::AlignCenter);
        PixmapLabel1->setWordWrap(false);
        PixmapLabel2 = new QLabel(TMOGUIOperation);
        PixmapLabel2->setObjectName(QString::fromUtf8("PixmapLabel2"));
        PixmapLabel2->setGeometry(QRect(210, 20, 161, 161));
        PixmapLabel2->setScaledContents(false);
        PixmapLabel2->setAlignment(Qt::AlignCenter);
        PixmapLabel2->setWordWrap(false);

        retranslateUi(TMOGUIOperation);
        QObject::connect(PushButton3, SIGNAL(clicked()), TMOGUIOperation, SLOT(accept()));
        QObject::connect(PushButton4, SIGNAL(clicked()), TMOGUIOperation, SLOT(reject()));

        QMetaObject::connectSlotsByName(TMOGUIOperation);
    } // setupUi

    void retranslateUi(QDialog *TMOGUIOperation)
    {
        TMOGUIOperation->setWindowTitle(QApplication::translate("TMOGUIOperation", "Image Operation", 0, QApplication::UnicodeUTF8));
        PushButton3->setText(QApplication::translate("TMOGUIOperation", "OK", 0, QApplication::UnicodeUTF8));
        PushButton4->setText(QApplication::translate("TMOGUIOperation", "Cancel", 0, QApplication::UnicodeUTF8));
        TextLabel1->setText(QApplication::translate("TMOGUIOperation", "Operation", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class TMOGUIOperation: public Ui_TMOGUIOperation {};
} // namespace Ui

QT_END_NAMESPACE

class TMOGUIOperation : public QDialog, public Ui::TMOGUIOperation
{
    Q_OBJECT

public:
    TMOGUIOperation(QWidget* parent = 0, const char* name = 0, bool modal = false, Qt::WindowFlags fl = 0);
    ~TMOGUIOperation();

protected slots:
    virtual void languageChange();

};

#endif // TMOGUIOPERATION_H
