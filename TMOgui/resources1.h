#ifndef TMOGUIMERGECOMPONENTS_H
#define TMOGUIMERGECOMPONENTS_H

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

class Ui_TMOGUIMergeComponents
{
public:
    QLabel *TextLabel1;
    QComboBox *ComboBox2;
    QLabel *TextLabel2;
    QLabel *TextLabel3;
    QComboBox *ComboBox3;
    QPushButton *PushButton1;
    QPushButton *PushButton2;
    QComboBox *ComboBox1;
    QLabel *PixmapLabel1;
    QLabel *PixmapLabel2;
    QLabel *PixmapLabel3;

    void setupUi(QDialog *TMOGUIMergeComponents)
    {
        if (TMOGUIMergeComponents->objectName().isEmpty())
            TMOGUIMergeComponents->setObjectName(QString::fromUtf8("TMOGUIMergeComponents"));
        TMOGUIMergeComponents->resize(636, 245);
        TextLabel1 = new QLabel(TMOGUIMergeComponents);
        TextLabel1->setObjectName(QString::fromUtf8("TextLabel1"));
        TextLabel1->setGeometry(QRect(70, 10, 30, 20));
        TextLabel1->setWordWrap(false);
        ComboBox2 = new QComboBox(TMOGUIMergeComponents);
        ComboBox2->setObjectName(QString::fromUtf8("ComboBox2"));
        ComboBox2->setGeometry(QRect(190, 210, 161, 21));
        TextLabel2 = new QLabel(TMOGUIMergeComponents);
        TextLabel2->setObjectName(QString::fromUtf8("TextLabel2"));
        TextLabel2->setGeometry(QRect(250, 10, 40, 21));
        TextLabel2->setWordWrap(false);
        TextLabel3 = new QLabel(TMOGUIMergeComponents);
        TextLabel3->setObjectName(QString::fromUtf8("TextLabel3"));
        TextLabel3->setGeometry(QRect(430, 10, 30, 21));
        TextLabel3->setWordWrap(false);
        ComboBox3 = new QComboBox(TMOGUIMergeComponents);
        ComboBox3->setObjectName(QString::fromUtf8("ComboBox3"));
        ComboBox3->setGeometry(QRect(360, 210, 161, 21));
        PushButton1 = new QPushButton(TMOGUIMergeComponents);
        PushButton1->setObjectName(QString::fromUtf8("PushButton1"));
        PushButton1->setGeometry(QRect(540, 30, 91, 30));
        PushButton2 = new QPushButton(TMOGUIMergeComponents);
        PushButton2->setObjectName(QString::fromUtf8("PushButton2"));
        PushButton2->setGeometry(QRect(540, 70, 91, 31));
        ComboBox1 = new QComboBox(TMOGUIMergeComponents);
        ComboBox1->setObjectName(QString::fromUtf8("ComboBox1"));
        ComboBox1->setGeometry(QRect(20, 210, 161, 21));
        PixmapLabel1 = new QLabel(TMOGUIMergeComponents);
        PixmapLabel1->setObjectName(QString::fromUtf8("PixmapLabel1"));
        PixmapLabel1->setGeometry(QRect(20, 40, 160, 160));
        QSizePolicy sizePolicy(static_cast<QSizePolicy::Policy>(1), static_cast<QSizePolicy::Policy>(1));
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(PixmapLabel1->sizePolicy().hasHeightForWidth());
        PixmapLabel1->setSizePolicy(sizePolicy);
        PixmapLabel1->setScaledContents(false);
        PixmapLabel1->setAlignment(Qt::AlignCenter);
        PixmapLabel1->setWordWrap(false);
        PixmapLabel2 = new QLabel(TMOGUIMergeComponents);
        PixmapLabel2->setObjectName(QString::fromUtf8("PixmapLabel2"));
        PixmapLabel2->setGeometry(QRect(190, 40, 160, 160));
        PixmapLabel2->setScaledContents(false);
        PixmapLabel2->setAlignment(Qt::AlignCenter);
        PixmapLabel2->setWordWrap(false);
        PixmapLabel3 = new QLabel(TMOGUIMergeComponents);
        PixmapLabel3->setObjectName(QString::fromUtf8("PixmapLabel3"));
        PixmapLabel3->setGeometry(QRect(360, 40, 160, 160));
        PixmapLabel3->setScaledContents(false);
        PixmapLabel3->setAlignment(Qt::AlignCenter);
        PixmapLabel3->setWordWrap(false);

        retranslateUi(TMOGUIMergeComponents);
        QObject::connect(PushButton1, SIGNAL(clicked()), TMOGUIMergeComponents, SLOT(accept()));
        QObject::connect(PushButton2, SIGNAL(clicked()), TMOGUIMergeComponents, SLOT(reject()));

        QMetaObject::connectSlotsByName(TMOGUIMergeComponents);
    } // setupUi

    void retranslateUi(QDialog *TMOGUIMergeComponents)
    {
        TMOGUIMergeComponents->setWindowTitle(QApplication::translate("TMOGUIMergeComponents", "Merge Components", 0, QApplication::UnicodeUTF8));
        TextLabel1->setText(QApplication::translate("TMOGUIMergeComponents", "Red", 0, QApplication::UnicodeUTF8));
        TextLabel2->setText(QApplication::translate("TMOGUIMergeComponents", "Green", 0, QApplication::UnicodeUTF8));
        TextLabel3->setText(QApplication::translate("TMOGUIMergeComponents", "Blue", 0, QApplication::UnicodeUTF8));
        PushButton1->setText(QApplication::translate("TMOGUIMergeComponents", "OK", 0, QApplication::UnicodeUTF8));
        PushButton2->setText(QApplication::translate("TMOGUIMergeComponents", "Cancel", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class TMOGUIMergeComponents: public Ui_TMOGUIMergeComponents {};
} // namespace Ui

QT_END_NAMESPACE

class TMOGUIMergeComponents : public QDialog, public Ui::TMOGUIMergeComponents
{
    Q_OBJECT

public:
    TMOGUIMergeComponents(QWidget* parent = 0, const char* name = 0, bool modal = false, Qt::WindowFlags fl = 0);
    ~TMOGUIMergeComponents();

protected slots:
    virtual void languageChange();

};

#endif // TMOGUIMERGECOMPONENTS_H
