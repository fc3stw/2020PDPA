#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <cstdlib>
#include <cmath>

GlobalPlacer::GlobalPlacer(Placement &placement)
	:_placement(placement)
{
    num_iters = 80;
    double max_module_dim = 0.;
    for(int m = 0; m<placement.numModules(); m++){
        max_module_dim = max(max_module_dim, placement.module(m).width());
        max_module_dim = max(max_module_dim, placement.module(m).height());
    }
    step_size = max_module_dim*2.;
    // step_size = (placement.rectangleChip().width() + placement.rectangleChip().height()) / 100.;
    // bin_res = 100;
    bin_res = sqrt(placement.numModules());
    // cout<<"Max step size = "<<step_size<<"\n";
    // cout<<"Max iterations = "<<num_iters<<"\n";
    // cout<<"bin resolution = "<<bin_res<<"x"<<bin_res
    //     <<", dimension = ("<<placement.rectangleChip().width() / (double)bin_res
    //     <<","<<placement.rectangleChip().height() / (double)bin_res<<")\n";
    
    // SA
    chip = placement.rectangleChip();
    T = 0.5;
    Tmin = 0.1;
    decay = 0.99;
    start_time = clock();
}

void GlobalPlacer::place()
{
	///////////////////////////////////////////////////////////////////
	// The following example is only for analytical methods.
	// if you use other methods, you can skip and delete it directly.
	//////////////////////////////////////////////////////////////////
    srand(0);
	ExampleFunction ef(_placement, bin_res); // require to define the object function and gradient function

    // module mi center coordinate = (x[2*i], x[2*i+1])
    vector<double> x(ef.dimension(), 0.);
    // initialize the solution vector
    double radius = 1.*max(_placement.rectangleChip().width(), _placement.rectangleChip().height())/2.;
    // double radius = 1000.;
    for(int m = 0; m<_placement.numModules(); m++){
        x[2*m] = rand()%(int)(2*radius) -radius + _placement.rectangleChip().centerX();
        x[2*m+1] = rand()%(int)(2*radius) -radius + _placement.rectangleChip().centerY();
    }

    NumericalOptimizer no(ef);
    no.setX(x); // set initial solution
    no.setNumIteration(num_iters); // user-specified parameter
    no.setStepSizeBound(step_size); // user-specified parameter
    no.solve(); // Conjugate Gradient solver

    cout << "Objective: " << no.objective() << endl;
	////////////////////////////////////////////////////////////////

    // place all modules
    for(int m = 0; m<_placement.numModules(); m++){
        _placement.module(m).setCenterPosition(no.x(2*m), no.x(2*m+1));
    }
}

void GlobalPlacer::SA_place()
{
    srand(0);
    initial_place();
    // backup(_placement.computeHpwl(), true);
    int num_moves = _placement.numModules() / 50.;
    // int num_moves = 1;
    int iter = 1;
    while(T>Tmin){        
        cout<<"#Iter "<<iter<<", \t";
        int failure = 0;
        for(int i = 0; i < num_moves; i++){
            double cost = _placement.computeHpwl();
            choose_modules();
            swap();
            if(_placement.computeHpwl() > cost){
                swap();
                failure++;
            }
        }
        double cost = _placement.computeHpwl();
        cout<<"HPWL = "<<cost<<",\t runtime = "<<get_time()<<"\t";
        if(failure > num_moves*0.9){
            T *= decay;
            cout<<"decay";
        }
        cout<<"\n"<<flush;
        iter++;
    }
    cout<<"Done\n";
}

void GlobalPlacer::initial_place()
{
    double usage_rate = 0.95;
    int num_modules = _placement.numModules();
    int num_cols = (double)num_modules / (double)_placement.numRows();
    double col_width = chip.width()*usage_rate / (double)num_cols;
    double row_height = _placement.row(0).height()*usage_rate;
    double x_offset = chip.width()*(1.-usage_rate)/2.;
    double y_offset = chip.height()*(1.-usage_rate)/2.;
    int row = 0;
    int col = 0;
    bool module_placed[num_modules] = {};
    for(int m = 0; m<num_modules; m++){
        if(module_placed[m]) continue;
        module_placed[m] = true;
        Module &module = _placement.module(m);
        double x = chip.left() + (double)(col+1)*col_width + x_offset;
        double y = chip.bottom() + (double)(row+1)*row_height + y_offset;
        module.setCenterPosition(x, y);
        col++;
        if(col == num_cols){
            col = 0;
            row++;
        }
    }
}

void GlobalPlacer::choose_modules()
{
    // choose 2 modules under distance constraint
    double x1;
    double y1;
    m2;
    double x2;
    double y2;
    do{
        m1 = rand()%_placement.numModules();
        m2 = rand()%_placement.numModules();
        x1 = _placement.module(m1).rectangle().centerX();
        y1 = _placement.module(m1).rectangle().centerY();
        x2 = _placement.module(m2).rectangle().centerX();
        y2 = _placement.module(m2).rectangle().centerY();
        // cout<<m1<<":("<<x1<<","<<y1<<"), "<<m2<<":("<<x2<<","<<y2<<")\n"<<endl;
    }while(abs(x1-x2) < chip.width()*T && abs(y1-y2) < chip.height()*T);
}

void GlobalPlacer::swap()
{
    double x1 = _placement.module(m1).rectangle().centerX();
    double y1 = _placement.module(m1).rectangle().centerY();
    double x2 = _placement.module(m2).rectangle().centerX();
    double y2 = _placement.module(m2).rectangle().centerY();
    _placement.module(m1).setCenterPosition(x2, y2);
    _placement.module(m2).setCenterPosition(x1, y1);
}

void GlobalPlacer::plotPlacementResult( const string outfilename, bool isPrompt )
{
    ofstream outfile( outfilename.c_str() , ios::out );
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT( outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop() );
    outfile << "EOF" << endl;
    outfile << "# modules" << endl << "0.00, 0.00" << endl << endl;
    for( size_t i = 0; i < _placement.numModules(); ++i ){
        Module &module = _placement.module(i);
        plotBoxPLT( outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height() );
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if( isPrompt ){
        char cmd[ 200 ];
        sprintf( cmd, "gnuplot %s", outfilename.c_str() );
        if( !system( cmd ) ) { cout << "Fail to execute: \"" << cmd << "\"." << endl; }
    }
}

void GlobalPlacer::plotBoxPLT( ofstream& stream, double x1, double y1, double x2, double y2 )
{
    stream << x1 << ", " << y1 << endl << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl << endl;
}
