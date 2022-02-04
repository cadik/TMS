#include "TMOGUIStyle.h"

#include <QComboBox>
#include <QPainter>
#include <QPainterPath>
#include <QPushButton>
#include <QStyleFactory>

TMOGUIStyle::TMOGUIStyle() : QProxyStyle(QStyleFactory::create("fusion"))
{
    setObjectName("TMOGUIStyle");
}

QPalette TMOGUIStyle::standardPalette() const
{
    if (!myStandardPalette.isBrushSet(QPalette::Disabled, QPalette::Mid))
    {
        QColor primary_color(96, 125, 139);
        QColor accent_green(205, 220, 57);
        QColor light_primary(207, 216, 220);
        QColor dark_primary(69, 90, 100);
        QColor darker_primary(35, 45, 51);
        QColor slightlyOpaqueBlack(0, 0, 0, 63);

        QColor text_primary(33, 33, 33);
        QColor text_secondary(117, 117, 117);
        QColor text_button(Qt::white);

        QColor divider(189, 189, 189);

        QPalette palette(primary_color);

        palette.setBrush(QPalette::Window, light_primary);
        palette.setBrush(QPalette::WindowText, text_primary);

        palette.setBrush(QPalette::Base, light_primary.lighter());
        palette.setBrush(QPalette::AlternateBase, light_primary);

        palette.setBrush(QPalette::PlaceholderText, text_secondary);
        palette.setBrush(QPalette::Text, text_primary);

        palette.setBrush(QPalette::BrightText, Qt::white);

        palette.setBrush(QPalette::Button, primary_color);
        palette.setBrush(QPalette::ButtonText, text_button);

        palette.setBrush(QPalette::Light, light_primary.lighter());
        palette.setBrush(QPalette::Midlight, light_primary);
        palette.setBrush(QPalette::Mid, dark_primary);
        palette.setBrush(QPalette::Dark, darker_primary);

        palette.setBrush(QPalette::Highlight, dark_primary);
        palette.setBrush(QPalette::HighlightedText, text_secondary.lighter());

        QBrush brush = palette.window();
        brush.setColor(brush.color().darker());
        //brush.setColor(accent_green);

        palette.setBrush(QPalette::Disabled, QPalette::ToolTipBase, light_primary);
        palette.setBrush(QPalette::Disabled, QPalette::ToolTipText, text_secondary);

        palette.setBrush(QPalette::Disabled, QPalette::WindowText, brush);
        palette.setBrush(QPalette::Disabled, QPalette::Text, brush);
        palette.setBrush(QPalette::Disabled, QPalette::ButtonText, brush);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush);
        palette.setBrush(QPalette::Disabled, QPalette::Button, brush);
        palette.setBrush(QPalette::Disabled, QPalette::Mid, brush);

        QPalette darkPalette;
        darkPalette.setColor(QPalette::Window, QColor(53, 53, 53));
        darkPalette.setColor(QPalette::WindowText, Qt::white);
        darkPalette.setColor(QPalette::Base, QColor(25, 25, 25));
        darkPalette.setColor(QPalette::AlternateBase, QColor(53, 53, 53));
        darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
        darkPalette.setColor(QPalette::ToolTipText, Qt::white);
        darkPalette.setColor(QPalette::Text, Qt::white);
        darkPalette.setColor(QPalette::Button, QColor(53, 53, 53));
        darkPalette.setColor(QPalette::ButtonText, Qt::white);
        darkPalette.setColor(QPalette::BrightText, QColor(0, 51, 51));
        darkPalette.setColor(QPalette::Link, QColor(42, 130, 218));
        darkPalette.setColor(QPalette::Shadow, QColor(10, 10, 10));

        darkPalette.setColor(QPalette::Highlight, QColor(0, 100, 100));
        darkPalette.setColor(QPalette::HighlightedText, Qt::black);

        darkPalette.setColor(QPalette::Disabled, QPalette::Shadow, Qt::black);

        //myStandardPalette = palette;
        myStandardPalette = darkPalette;
    }

    return myStandardPalette;
}

void TMOGUIStyle::setTexture(QPalette &palette, QPalette::ColorRole role,
                             const QImage &image)
{
    for (int i = 0; i < QPalette::NColorGroups; ++i)
    {
        QBrush brush(image);
        brush.setColor(palette.brush(QPalette::ColorGroup(i), role).color());
        palette.setBrush(QPalette::ColorGroup(i), role, brush);
    }
}

void TMOGUIStyle::polish(QWidget *widget)
{
    widget->setPalette(standardPalette());
    if (qobject_cast<QPushButton *>(widget) || qobject_cast<QComboBox *>(widget))
        widget->setAttribute(Qt::WA_Hover, true);
}

void TMOGUIStyle::unpolish(QWidget *widget)
{
    if (qobject_cast<QPushButton *>(widget) || qobject_cast<QComboBox *>(widget))
        widget->setAttribute(Qt::WA_Hover, false);
}

int TMOGUIStyle::pixelMetric(PixelMetric metric,
                             const QStyleOption *option,
                             const QWidget *widget) const
{
    switch (metric)
    {
    case QStyle::PM_MenuBarItemSpacing:
        return 7;
    case PM_MenuBarHMargin:
        return 13;
    case PM_MenuBarVMargin:
        return 0;
    case PM_MenuBarPanelWidth:
        return 1;
    case PM_ComboBoxFrameWidth:
        return 8;
    /*case PM_ScrollBarExtent:
            return QProxyStyle::pixelMetric(metric, option, widget) + 4;*/
    default:
        return QProxyStyle::pixelMetric(metric, option, widget);
    }
}

int TMOGUIStyle::styleHint(StyleHint hint, const QStyleOption *option,
                           const QWidget *widget,
                           QStyleHintReturn *returnData) const
{
    switch (hint)
    {
    case SH_DitherDisabledText:
        return int(false);
    case SH_EtchDisabledText:
        return int(true);
    default:
        return QProxyStyle::styleHint(hint, option, widget, returnData);
    }
}

void TMOGUIStyle::drawPrimitive(PrimitiveElement element,
                                const QStyleOption *option,
                                QPainter *painter,
                                const QWidget *widget) const
{
    switch (element)
    {
    case PE_PanelButtonCommand:
    {
        int delta = (option->state & State_MouseOver) ? 64 : 0;
        QColor slightlyOpaqueBlack(0, 0, 0, 63);
        QColor semiTransparentWhite(255, 255, 255, 127 + delta);
        QColor semiTransparentBlack(0, 0, 0, 127 - delta);

        int x, y, width, height;
        option->rect.getRect(&x, &y, &width, &height);

        QPainterPath roundRect = roundRectPath(option->rect);

        int radius = qMin(width, height) / 4;

        QBrush brush;

        bool darker;

        const QStyleOptionButton *buttonOption =
            qstyleoption_cast<const QStyleOptionButton *>(option);
        if (buttonOption && (buttonOption->features & QStyleOptionButton::Flat))
        {
            brush = option->palette.window();
            darker = (option->state & (State_Sunken | State_On));
        }
        else
        {
            if (option->state & (State_Sunken | State_On))
            {
                brush = option->palette.mid();
                darker = !(option->state & State_Sunken);
            }
            else
            {
                brush = option->palette.button();
                darker = false;
            }
        }

        painter->save();
        painter->setRenderHint(QPainter::Antialiasing, true);
        painter->fillPath(roundRect, brush);
        if (darker)
            painter->fillPath(roundRect, slightlyOpaqueBlack);

        int penWidth;
        if (radius < 10)
            penWidth = 3;
        else if (radius < 20)
            penWidth = 5;
        else
            penWidth = 7;

        QPen topPen(semiTransparentWhite, penWidth);
        QPen bottomPen(semiTransparentBlack, penWidth);

        if (option->state & (State_Sunken | State_On))
            qSwap(topPen, bottomPen);

        int x1 = x;
        int x2 = x + radius;
        int x3 = x + width - radius;
        int x4 = x + width;

        if (option->direction == Qt::RightToLeft)
        {
            qSwap(x1, x4);
            qSwap(x2, x3);
        }

        QPolygon topHalf;
        topHalf << QPoint(x1, y)
                << QPoint(x4, y)
                << QPoint(x3, y + radius)
                << QPoint(x2, y + height - radius)
                << QPoint(x1, y + height);

        painter->setClipPath(roundRect);
        painter->setClipRegion(topHalf, Qt::IntersectClip);
        painter->setPen(topPen);
        painter->drawPath(roundRect);

        QPolygon bottomHalf = topHalf;
        bottomHalf[0] = QPoint(x4, y + height);

        painter->setClipPath(roundRect);
        painter->setClipRegion(bottomHalf, Qt::IntersectClip);
        painter->setPen(bottomPen);
        painter->drawPath(roundRect);

        painter->setPen(option->palette.windowText().color());
        painter->setClipping(false);
        painter->drawPath(roundRect);

        painter->restore();
    }
    break;
    default:
        QProxyStyle::drawPrimitive(element, option, painter, widget);
    }
}

void TMOGUIStyle::drawControl(ControlElement element,
                              const QStyleOption *option,
                              QPainter *painter,
                              const QWidget *widget) const
{
    switch (element)
    {
    case CE_PushButtonLabel:
    {
        QStyleOptionButton myButtonOption;
        const QStyleOptionButton *buttonOption =
            qstyleoption_cast<const QStyleOptionButton *>(option);
        if (buttonOption)
        {
            myButtonOption = *buttonOption;
            if (myButtonOption.palette.currentColorGroup() != QPalette::Disabled)
            {
                if (myButtonOption.state & (State_Sunken | State_On))
                {
                    myButtonOption.palette.setBrush(QPalette::ButtonText,
                                                    myButtonOption.palette.brightText());
                }
            }
        }
        QProxyStyle::drawControl(element, &myButtonOption, painter, widget);
    }
    break;
    default:
        QProxyStyle::drawControl(element, option, painter, widget);
    }
}

QPainterPath TMOGUIStyle::roundRectPath(const QRect &rect)

{
    int shortSide = qMin(rect.width(), rect.height());
    int radius = shortSide / 4;
    int diam = 2 * radius;

    int x1, y1, x2, y2;
    rect.getCoords(&x1, &y1, &x2, &y2);

    QPainterPath path;
    path.moveTo(x2, y1 + radius);
    path.arcTo(QRect(x2 - diam, y1, diam, diam), 0.0, +90.0);
    path.lineTo(x1 + radius, y1);
    path.arcTo(QRect(x1, y1, diam, diam), 90.0, +90.0);
    path.lineTo(x1, y2 - radius);
    path.arcTo(QRect(x1, y2 - diam, diam, diam), 180.0, +90.0);
    path.lineTo(x1 + radius, y2);
    path.arcTo(QRect(x2 - diam, y2 - diam, diam, diam), 270.0, +90.0);
    path.closeSubpath();
    return path;
}
