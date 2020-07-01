#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void RGBtoLab(double R, double G, double B, double &L, double &a, double &b);
void XYZtoLab(double X, double Y, double Z, double &L, double &a, double &b);
void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z);
double F(double input);

