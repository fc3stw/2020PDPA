#define _GLIBCXX_USE_CXX11_ABI 0
#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <ctime>

class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement);
	void place();
    void plotPlacementResult( const string outfilename, bool isPrompt = false );
    // SA
    void SA_place();
    void initial_place();
    void choose_modules();
    void swap();
    double get_time() const {return (double)(clock() - start_time) / CLOCKS_PER_SEC;}

private:
    Placement& _placement;
    void plotBoxPLT( ofstream& stream, double x1, double y1, double x2, double y2 );
    int num_iters;
    double step_size;
    int bin_res;
    // SA
    Rectangle chip;
    double T;
    double Tmin;
    double decay;
    int m1, m2;
    vector<Module> best_modules;
    double start_time;
};

#endif // GLOBALPLACER_H
