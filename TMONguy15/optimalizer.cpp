#include <iomanip>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <nlopt.hpp>

typedef struct {
    double a, b, c;
} my_constraint_data;

typedef struct {
    double Rs,Gs,Bs,Ds,Rv,Gv,Bv,Dv,Rj,Gj,Bj,Dj,Rz,Gz,Bz,Dz;
} my_constraint_data1;

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) // milization of a function
{
    my_constraint_data1 *d = reinterpret_cast<my_constraint_data1*>(my_func_data);
    double Rs = d->Rs, Gs = d->Gs, Bs = d->Bs,  Rj = d->Rj, Gj = d->Gj, Bj = d->Bj,  Rv = d->Rv, Gv = d->Gv, Bv = d->Bv,  Rz = d->Rz, Gz = d->Gz, Bz = d->Bz, Ds = d->Ds, Dv = d->Dv, Dj = d->Dj, Dz = d->Dz;
/*if (!grad.empty()) {
grad[0] = 0.0;
grad[1] = 0.5 / sqrt(x[1]);
}*/
return pow(( (x[0]*Rs) + (x[1]*Gs) + (x[2]*Bs) )- Ds,2) + pow(((x[0]*Rj) + (x[1]*Gj) + (x[2]*Bj)) - Dj,2) + pow(((x[0]*Rv) + (x[1]*Gv) + (x[2]*Bv)) - Dv,2) + pow(((x[0]*Rz) + (x[1]*Gz) + (x[2]*Bz)) - Dz,2);
}
double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)// minimalization constraints r+g+b=1
{
my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
double a = d->a, b = d->b, c = d->c;
/*if (!grad.empty()) {
grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
grad[1] = -1.0;
}*/
return (a*x[0]) + (b*x[1]) + (x[2]*c) - 1;
}

double myvconstraint1(const std::vector<double> &x, std::vector<double> &grad, void *data)  // minimalization constraints r+g+b=1
{
    my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
    double a = d->a, b = d->b, c = d->c;
/*if (!grad.empty()) {
grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
grad[1] = -1.0;
}*/
    return (a*x[0]) + (b*x[1]) + (x[2]*c) - 1;
}

void optimalize(double Rs, double Gs, double Bs, double Rv, double Gv, double Bv, double Rj, double Gj, double Bj, double Rz, double Gz, double Bz, double &Rres, double &Gres, double &Bres, double ds, double dv, double dj, double dz){
    //std::cout << "zde" << "\n";
    //printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n", Rs,Gs,Bs,ds,Rv,Gv,Bv,dv,Rj,Gj,Bj,dj,Rz,Gz,Bz,dz);
    nlopt::opt opt(nlopt::LN_COBYLA, 3);
    std::vector<double> lb(3);
    lb[0] = 0; lb[1] = 0; lb[2] = 0;
    opt.set_lower_bounds(lb);

    my_constraint_data1 data1[1];
    data1[0].Rs = Rs;
    data1[0].Gs = Gs;
    data1[0].Bs = Bs;
    data1[0].Ds = ds;
    data1[0].Rj = Rj;
    data1[0].Gj = Gj;
    data1[0].Bj = Bj;
    data1[0].Dj = dj;
    data1[0].Rv = Rv;
    data1[0].Gv = Gv;
    data1[0].Bv = Bv;
    data1[0].Dv = dv;
    data1[0].Rz = Rz;
    data1[0].Gz = Gz;
    data1[0].Bz = Bz;
    data1[0].Dz = dz;

    opt.set_min_objective(myvfunc, &data1[0]);
    my_constraint_data data[1] = {{1,1,1}};
    opt.add_equality_constraint(myvconstraint, &data[0], 1e-8);
    //opt.add_inequality_constraint(myvconstraint1, &data[0], 1e-8);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(3);
    x[0] = 0.3; x[1] = 0.59; x[2] = 0.11;
    double minf;

    try{
    nlopt::result result = opt.optimize(x, minf);
    std::cerr << "found minimum at f(" << x[0] << "," << x[1] <<  "," << x[2] << ") = " << std::setprecision(10) << minf << std::endl;
    Rres = x[0];
    Gres = x[1];
    Bres = x[2];
    //printf("End \n");
    }
    catch(std::exception &e) {
    //std::cerr << "nlopt failed: " << e.what() << std::endl;
    }

}

/*void optimalize(double Rs, double Gs, double Bs, double Rv, double Gv, double Bv, double Rj, double Gj, double Bj, double Rz, double Gz, double Bz, double &Rres, double &Gres, double &Bres, double ds, double dv, double dj, double dz)
{

}*/
