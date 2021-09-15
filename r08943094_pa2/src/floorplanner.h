#include "BSTree.h"
#include <fstream>
#include <map>
#include <cassert>
#include <ctime>
#include <queue>
#include <random>
#include <cmath>
#include <tuple>

using namespace std;

// Definition of dimensions: width->horizontal, height->vertical

typedef map<string, int> Name2Id;
typedef vector<Block*> BlockList;
typedef tuple<int, int, int, bool> OP;

class Floorplanner
{
    double alpha;   // evaluation: alpha*A + (1-alpha)*WL
    int outline_width;
    int outline_height;
    int num_blocks;
    int num_terminals;
    int num_nets;
    int chip_width;
    int chip_height;
    double WL;
    int num_samples;
    double norm_area;
    double norm_wl;
    BlockList block_list;
    vector<Terminal*> terminal_list;
    vector<Net*> net_list;
    Name2Id blkname2id;
    Name2Id terminalname2id;
    Contour contour;
    BST *current_bst;
    BST *best_bst;
    BST *previous_bst;
    clock_t start_time;
    default_random_engine rand_gen;
    uniform_real_distribution<double> unif;

    // query
    bool place_legal(Node *parent, Block *block, bool LR=false){
        if(parent==nullptr) return block->get_width() <= outline_width;
        else if(LR) return block_list[parent->id]->get_low_xy().first + block->get_width() <= outline_width;
        else return block_list[parent->id]->get_top_xy().first + block->get_width() <= outline_width;
    }
    vector<int> get_highest_placed_block_ids() const;
    int get_chip_width() const{return chip_width;}
    int get_chip_height() const{return chip_height;}
    bool is_legal() const{return get_chip_width() <= outline_width && get_chip_height() <= outline_height;}
    int get_area() const{return chip_width * chip_height;}
    double get_wl() const{return WL;}
    double get_cost() const{return alpha*(double)get_area() + (1.-alpha)*get_wl();}
    double get_outline_ar() const{return (double)outline_height/(double)outline_width;}
    double get_chip_ar() const{return (double)chip_height/(double)chip_width;}
    double get_normalized_cost() const{
        if(get_chip_width() > outline_width || get_chip_height() > outline_height)
            return alpha*(double)get_area()/norm_area + (1.-alpha)*get_wl()/norm_wl + 10.*abs(get_outline_ar()-get_chip_ar());
        else
            return alpha*(double)get_area()/norm_area + (1.-alpha)*get_wl()/norm_wl + abs(get_outline_ar()-get_chip_ar());
    }
    // method
    void place_block(Node *node, Block *block, bool LR=false);
    void update_all_blocks(bool sampling_node=false);
    void update_block_dfs(Node *node);
    void calculate_wl();    
    void backup_floorplan(bool force=false);
    void restore_floorplan(BST *tree);
    void restore_block_node(Node *node);
    void random_operation();
    void sample_normalized_cost();
    void update_normalized_cost();
    // algorithm
    void initial_floorplan(BlockList &unplaced_blocks);
    void SA();
    // sanity check
    void print_placed_blocks();
    void print_bst_dfs(Node *node);
    bool check_outline();
    void check_overlap();
    bool block_overlap(Block *b1, Block *b2);
public:
    Floorplanner():
    alpha(0.),
    outline_width(0),
    outline_height(0),
    num_blocks(0),
    num_terminals(0),
    num_nets(0),
    chip_width(0),
    chip_height(0),
    WL(0.),
    num_samples(0),
    norm_area(0.),
    norm_wl(0.),
    current_bst(new BST),
    best_bst(new BST),
    previous_bst(new BST),
    start_time(clock()),
    rand_gen(default_random_engine(0)),
    unif(uniform_real_distribution<double>(0., 0.999999))
    {
        contour.reset();
    }
    ~Floorplanner(){clear();}
    // necessary inputs
    void parse_block(fstream &block_file);
    void parse_net(fstream &net_file);
    void set_alpha(double val){
        if(val>1||val<0){
            cerr<<"Alpha value error: has to be within interval [0,1]\n";
            exit(1);
        }
        alpha = val;
        cout<<"Set alpha value="<<alpha<<endl;
    }
    // floorplanning
    void floorplan();
    // stream out files
    void plot(string file_name);
    void write_result(fstream &output_file);
    // show info
    void print_summary(int verbose=0);
    double get_time() const {return (double)(clock() - start_time) / CLOCKS_PER_SEC;}
    void clear();
};
