#ifndef TMOGUISTYLE_H
#define TMOGUISTYLE_H

#include <QProxyStyle>
#include <QObject>

class TMOGUIStyle : public QProxyStyle
{
public:
    TMOGUIStyle();

    QPalette standardPalette() const override;



    void polish(QWidget *widget) override;
    void unpolish(QWidget *widget) override;
    int pixelMetric(PixelMetric metric, const QStyleOption *option,
                    const QWidget *widget) const override;
    int styleHint(StyleHint hint, const QStyleOption *option,
                  const QWidget *widget, QStyleHintReturn *returnData) const override;
    void drawPrimitive(PrimitiveElement element, const QStyleOption *option,
                       QPainter *painter, const QWidget *widget) const override;
    void drawControl(ControlElement control, const QStyleOption *option,
                     QPainter *painter, const QWidget *widget) const override;




private:
    static void setTexture(QPalette &palette, QPalette::ColorRole role,
                               const QImage &image);
    static QPainterPath roundRectPath(const QRect &rect);
    mutable QPalette myStandardPalette;
};

#endif // TMOGUISTYLE_H
