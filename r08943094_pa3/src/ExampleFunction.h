#define _GLIBCXX_USE_CXX11_ABI 0
#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Placement.h"
#include "Rectangle.h"

class ExampleFunction : public NumericalOptimizerInterface
{
    // density
    double lambda; // density weight
    double Mb; // desired density for a bin
    double alpha; // sigmoid parameter
    double bin_res; // bin resolution, num of bins per dimension
    Rectangle chip;
    vector<Rectangle> modules;
    vector<Rectangle> bins;
    vector<double> total_OxOy;
    // wl
    double gamma; // lse parameter
    vector<Net*> nets;
    vector<double> net_total_exp_x;
    vector<double> net_total_exp_minx;
    vector<double> net_total_exp_y;
    vector<double> net_total_exp_miny;
    vector<double> wl_x_grad;
    vector<double> wl_y_grad;
    // others
    int iter;
public:
    ExampleFunction(Placement &placement, int bin_res_in);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();

    // terms for objective function
    // density
    double sigmoid(double x);
    double Ox(Rectangle &bin, Rectangle &module);
    double Oy(Rectangle &bin, Rectangle &module);
    double total_density();
    // density gradient
    double dOxdx(Rectangle &bin, Rectangle &module);
    double dOydy(Rectangle &bin, Rectangle &module);
    double total_density_x_grad(size_t idx);
    double total_density_y_grad(size_t idx);
    // wirelength, suppose all pins are at the center of the module
    double lse(vector<double> vals, double scale);
    double total_wirelength();
    // wirelength gradient
    double dLSEdx(int net_id, double x);
    double dLSEdy(int net_id, double y);
    double total_wl_x_grad(size_t idx);
    double total_wl_y_grad(size_t idx);
    void total_wl_grad();
    // others
    bool overlap(Rectangle &bin, Rectangle &module);
    vector<int> get_overlapped_bins(Rectangle &module);
    double get_avg_dist();
    void print();
};
#endif // EXAMPLEFUNCTION_H
