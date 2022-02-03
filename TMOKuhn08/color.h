#ifndef COLOR_H
#define COLOR_H

#include <ostream>
#include <iostream>
namespace color
{
    enum ColorCode
    {
        RED = 31,
        GREEN = 32,
        BLUE = 34,
        DEFAULT = 39,
    };

    enum TypeCode
    {
        BOLD = 1,
        NORMAL = 0
    };

    class CColor
    {
        ColorCode code;
        TypeCode tcode;

    public:
        CColor(ColorCode c, TypeCode tc)
        {
            code = c;
            tcode = tc;
        }
        CColor(ColorCode c)
        {
            code = c;
            tcode = color::NORMAL;
        }
        friend std::ostream &
        operator<<(std::ostream &os, const CColor &mod)
        {
            return os << "\033[" << mod.tcode << ";" << mod.code << "m";
        }
    };

    static CColor red(color::RED);
    static CColor blue(color::BLUE);
    static CColor green(color::GREEN);

    static CColor bred(color::RED, BOLD);
    static CColor bgreen(color::GREEN, BOLD);
    static CColor bblue(color::BLUE, BOLD);

    static CColor reset(color::DEFAULT);
}

#endif