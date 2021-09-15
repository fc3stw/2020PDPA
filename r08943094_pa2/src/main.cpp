#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include "floorplanner.h"
using namespace std;

int main(int argc, char** argv)
{
    fstream block_file, net_file, output_file;
    bool gui_flag = false;
    string plot_name = "plot.png";
    cout<<fixed;

    if (argc >= 5) {
        block_file.open(argv[2], ios::in);
        net_file.open(argv[3], ios::in);
        output_file.open(argv[4], ios::out);
        if (!block_file) {
            cerr << "Cannot open the .block file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!net_file) {
            cerr << "Cannot open the .nets file \"" << argv[3]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output_file) {
            cerr << "Cannot open the output file \"" << argv[4]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if(argc >= 6 && argv[5]==string("-gui")){
            gui_flag = true;
            cout<<"Plot result will be saved\n";
            if(argc>=7){
                plot_name = argv[6];
            }
        }
    }
    else {
        cerr << "Usage: ./fp <α value> <input.block name> <input.net name> <output file name> [-gui [plot_name.png]]" << endl;
        exit(1);
    }

    Floorplanner* floorplanner = new Floorplanner();
    floorplanner->set_alpha(stod(argv[1]));
    floorplanner->parse_block(block_file);
    floorplanner->parse_net(net_file);
    floorplanner->floorplan();
    if(gui_flag) floorplanner->plot(plot_name);
    floorplanner->print_summary(1);
    floorplanner->write_result(output_file);

    return 0;
}
