#include "ExampleFunction.h"
#include <cmath>
#include <cassert>
#include <string>
#include <cstdlib>

ExampleFunction::ExampleFunction(Placement &placement, int bin_res_in)
{
    // density
    bin_res = bin_res_in;
    chip = placement.rectangleChip();
    for(int i = 0; i<placement.numModules(); i++){
        modules.push_back(placement.module(i).rectangle());
    }
    // build bins
    double bin_width = chip.width() / bin_res;
    double bin_height = chip.height() / bin_res;
    double chip_x = chip.left();
    double chip_y = chip.bottom();
    for(int r = 0; r<bin_res; r++){
        for(int c = 0; c<bin_res; c++){
            bins.push_back(
                Rectangle(bin_width*c + chip_x,
                            bin_height*r + chip_y,
                            bin_width*(c+1) + chip_x,
                            bin_height*(r+1) + chip_y)
                );
        }
    }
    assert(bins.size() == (int)(bin_res*bin_res));
    lambda = pow(10, -18);
    // lambda = 0.;
    // alpha = 0.005;
    alpha = 1. / pow(log((bin_width + bin_height)/2.), 3.);
    // alpha = 1. / sqrt((bin_width + bin_height)/2.);
    double total_cell_area = 0.;
    for(int m = 0; m<modules.size(); m++){
        total_cell_area += modules[m].width()*modules[m].height();
    }
    // Mb = 0.5;
    Mb = modules.size() / pow(bin_res, 2);
    // Mb = total_cell_area / (bin_width*bin_height*pow(bin_res, 2));
    
    // wl
    nets.clear();
    for(int n = 0; n<placement.numNets(); n++){
        nets.push_back(&(placement.net(n)));
    }
    double ub = pow(10, 1); // upper bound for exp(x/gamma)
    // gamma = 1000;
    gamma = max(chip.width(), chip.height()) / log(ub);

    iter = 1;
    print();
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    cout<<"\n#Iter: "<<iter++<<"\n";
    // evaluate objective function
    evaluateF(x, f);
    cout<<"Average distance from chip center = "<<get_avg_dist()<<endl;

    // evaluate gradient
    // density
    // cout<<"computing bin OxOy\n"<<flush;
    total_OxOy = vector<double>(bins.size(), 0.);
    for(int m = 0; m<modules.size(); m++){
        Rectangle module = modules[m];
        vector<int> overlapped_bins = get_overlapped_bins(module);
        for(int i = 0; i<overlapped_bins.size(); i++){
            int b = overlapped_bins[i];
            Rectangle &bin = bins[b];
            total_OxOy[b] += Ox(bin, module)*Oy(bin, module);
        }
    }
    // wl
    // cout<<"computing net wl\n"<<flush;
    net_total_exp_x = vector<double>(nets.size(), 0.);
    net_total_exp_minx = vector<double>(nets.size(), 0.);
    net_total_exp_y = vector<double>(nets.size(), 0.);
    net_total_exp_miny = vector<double>(nets.size(), 0.);
    for(int n = 0; n<nets.size(); n++){
        Net *net = nets[n];
        for(int p = 0; p<net->numPins(); p++){
            int m_id = net->pin(p).moduleId();
            double xp = modules[m_id].centerX();
            double yp = modules[m_id].centerY();
            net_total_exp_x[n] += exp(xp/gamma);
            net_total_exp_minx[n] += exp(-xp/gamma);
            net_total_exp_y[n] += exp(yp/gamma);
            net_total_exp_miny[n] += exp(-yp/gamma);
        }
    }
    wl_x_grad = vector<double>(modules.size(), 0.);
    wl_y_grad = vector<double>(modules.size(), 0.);
    total_wl_grad();

    // grad
    // cout<<"computing gradient\n"<<flush;
    for(int i = 0; i<modules.size(); i++){
        // g[2*i] = total_wl_x_grad(i) + lambda*total_density_x_grad(i);
        // g[2*i+1] = total_wl_y_grad(i) + lambda*total_density_y_grad(i);
        g[2*i] = wl_x_grad[i] + lambda*total_density_x_grad(i);
        g[2*i+1] = wl_y_grad[i] + lambda*total_density_y_grad(i);
        // g[2*i] = lambda*total_density_x_grad(i);
        // g[2*i+1] = lambda*total_density_y_grad(i);
    }
    lambda *= 2.;
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
    // cout<<"calculate objective\n"<<flush;
    // update coordinate for all modules
    // module bottom-left coordinate = (x[2*i], x[2*i+1])
    for(int i = 0; i<modules.size(); i++){
        modules[i] = Rectangle(x[2*i] - 0.5*modules[i].width(),
                                x[2*i+1] - 0.5*modules[i].height(),
                                x[2*i] + 0.5*modules[i].width(),
                                x[2*i+1] + 0.5*modules[i].height());
    }
    // objective function
    f = total_wirelength() + lambda*total_density();
    // f = lambda*total_density();
}

unsigned ExampleFunction::dimension()
{
    return modules.size()*2; // num_blocks*2 
    // each two dimension represent the X and Y dimensions of each block
}

inline double ExampleFunction::sigmoid(double x)
{
    return 1./(1. + exp(-alpha*x));
}

inline double ExampleFunction::Ox(Rectangle &bin, Rectangle &module)
{
    double x = module.centerX() - bin.centerX();
    double w = bin.width();
    return sigmoid(0.5*w + x)*sigmoid(0.5*w - x);
}

inline double ExampleFunction::Oy(Rectangle &bin, Rectangle &module)
{
    double y = module.centerY() - bin.centerY();
    double h = bin.height();
    return sigmoid(0.5*h + y)*sigmoid(0.5*h - y);
}

double ExampleFunction::total_density()
{
    double total_density = 0.;
    for(int m = 0; m<modules.size(); m++){
        Rectangle &module = modules[m];
        vector<int> overlapped_bins = get_overlapped_bins(module);
        for(int i = 0; i<overlapped_bins.size(); i++){
            int b = overlapped_bins[i];
            Rectangle &bin = bins[b];
            total_density += pow(Ox(bin, module)*Oy(bin, module) - Mb, 2);
        }
    }
    return total_density;
}

inline double ExampleFunction::dOxdx(Rectangle &bin, Rectangle &module)
{
    double x = module.centerX() - bin.centerX();
    double w = bin.width();
    double grad = alpha * pow(Ox(bin, module), 2) * (1./sigmoid(0.5*w+x) - 1./sigmoid(0.5*w-x));
    return grad;
}

inline double ExampleFunction::dOydy(Rectangle &bin, Rectangle &module)
{
    double y = module.centerY() - bin.centerY();
    double h = bin.height();
    double grad = alpha * pow(Oy(bin, module), 2) * (1./sigmoid(0.5*h+y) - 1./sigmoid(0.5*h-y));
    return grad;
}

double ExampleFunction::total_density_x_grad(size_t idx)
{
    // total x gradient of module #idx
    // partial derivatives in terms of x-idx 
    double total_grad = 0.;
    Rectangle &module = modules.at(idx);
    vector<int> overlapped_bins = get_overlapped_bins(module);
    for(int i = 0; i<overlapped_bins.size(); i++){
        int b = overlapped_bins[i];
        Rectangle &bin = bins[b];
        total_grad += 2*(total_OxOy[b]-Mb)
                    *dOxdx(bin, module)*Oy(bin, module);
    }
    return total_grad;
}

double ExampleFunction::total_density_y_grad(size_t idx)
{
    // total y gradient of module #idx
    // partial derivatives in terms of x-idx 
    double total_grad = 0.;
    Rectangle &module = modules.at(idx);
    vector<int> overlapped_bins = get_overlapped_bins(module);
    for(int i = 0; i<overlapped_bins.size(); i++){
        int b = overlapped_bins[i];
        Rectangle &bin = bins[b];
        total_grad += 2*(total_OxOy[b]-Mb)
                    *dOydy(bin, module)*Ox(bin, module);
    }
    return total_grad;
}

inline double ExampleFunction::lse(vector<double> vals, double scale)
{
    double lse = 0.;
    for(int i = 0; i<vals.size(); i++){
        lse += exp(scale*vals[i]/gamma);
    }
    return log(lse);
}

double ExampleFunction::total_wirelength()
{
    double total_wl = 0.;
    for(int n = 0; n<nets.size(); n++){
        Net *net = nets[n];
        vector<double> xs(net->numPins(), 0.);
        vector<double> ys(net->numPins(), 0.);
        for(int p = 0; p<net->numPins(); p++){
            int m_id = net->pin(p).moduleId();
            xs[p] = modules[m_id].centerX();
            ys[p] = modules[m_id].centerY();
            total_wl += lse(xs, 1) + lse(xs, -1) + lse(ys, 1) + lse(ys, -1);
        }
    }
    return gamma * total_wl;
}

double ExampleFunction::dLSEdx(int net_id, double x)
{
    return exp(x/gamma) / (gamma*net_total_exp_x[net_id]) + -exp(-x/gamma) / (gamma*net_total_exp_minx[net_id]);
}

double ExampleFunction::dLSEdy(int net_id, double y)
{
    return exp(y/gamma) / (gamma*net_total_exp_y[net_id]) + -exp(-y/gamma) / (gamma*net_total_exp_miny[net_id]);
}

double ExampleFunction::total_wl_x_grad(size_t idx)
{
    double x = modules[idx].centerX();
    double total_grad = 0.;
    for(int n = 0; n<nets.size(); n++){
        total_grad += dLSEdx(n, x);
    }
    return gamma*total_grad;
}

double ExampleFunction::total_wl_y_grad(size_t idx)
{
    double y = modules[idx].centerY();
    double total_grad = 0.;
    for(int n = 0; n<nets.size(); n++){
        total_grad += dLSEdy(n, y);
    }
    return gamma*total_grad;
}

void ExampleFunction::total_wl_grad()
{
    for(int n = 0; n<nets.size(); n++){
        Net *net = nets[n];
        for(int p = 0; p<net->numPins(); p++){
            Pin pin = net->pin(p);
            int m_id = pin.moduleId();
            double x = modules[m_id].centerX();
            double y = modules[m_id].centerY();
            wl_x_grad[m_id] += dLSEdx(n, x);
            wl_y_grad[m_id] += dLSEdy(n, y);
        }
    }
}

inline bool ExampleFunction::overlap(Rectangle &bin, Rectangle &module)
{
    bool horizontal_overlap1 = module.left() > bin.left() && module.left() < bin.right()
                            || module.right() > bin.left() && module.right() < bin.right();
    bool horizontal_overlap2 = bin.left() > module.left() && bin.left() < module.right()
                            || bin.right() > module.left() && bin.right() < module.right();
    bool vertical_overlap1 = module.bottom() > bin.bottom() && module.bottom() < bin.top()
                            || module.top() > bin.bottom() && module.top() < bin.top();
    bool vertical_overlap2 = bin.bottom() > module.bottom() && bin.bottom() < module.top()
                            || bin.top() > module.bottom() && bin.top() < module.top();
    // return horizontal_overlap1 && vertical_overlap1 || horizontal_overlap2 && vertical_overlap2;                                 
    return true;
}

vector<int> ExampleFunction::get_overlapped_bins(Rectangle &module)
{
    vector<int> overlapped_bins;
    int l = max(floor(bin_res * (module.left() - chip.left()) / chip.width()), 0.);
    int r = min(floor(bin_res * (module.right() - chip.left()) / chip.width()), bin_res-1.);
    int b = max(floor(bin_res * (module.bottom() - chip.bottom()) / chip.height()), 0.);
    int t = min(floor(bin_res * (module.top() - chip.bottom()) / chip.height()), bin_res-1.);
    for(int i = b; i<=t; i++){ // row
        for(int j = l; j<=r; j++){ // col
            overlapped_bins.push_back(i*(int)bin_res + j);
        }
    }
    return overlapped_bins;
}

double ExampleFunction::get_avg_dist()
{
    double dist = 0.;
    for(int m = 0; m<modules.size(); m++){
        Rectangle &module = modules[m];
        dist += sqrt(pow(module.centerX() - chip.centerX(), 2) + pow(module.centerY() - chip.centerY(), 2));
    }
    return dist / (double)modules.size();
}

void ExampleFunction::print(){
    cout<<"Mb = "<<Mb<<"\n";
    cout<<"lambda = "<<lambda<<"\n";
    cout<<"alpha = "<<alpha<<"\n";
    cout<<"gamma = "<<gamma<<"\n";
}
