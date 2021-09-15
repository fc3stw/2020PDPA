#include "module.h"
#include <climits>

using namespace std;

// B*-tree and node orientation alone will decide the entire floorplan
class BST
{
    Node *root;
    double cost;
    int num_nodes;
public:
    BST() : root(nullptr), cost(INT_MAX) {}
    // getter
    Node* get_root() const{return root;}
    Node* get_rightmost() const{
        Node *n = root;
        while(n->right){n = n->right;}
        return n;
    }
    double get_cost() const{return cost;}
    int get_num_nodes() const{return num_nodes;}
    // setter
    void set_root(Node *node){root = node;}
    void set_cost(double val){cost = val;}
    // method
    bool remove(Node *node);
    bool append(Node *parent, Node *node, bool LR=false);
    void rotate(Node *node);
    void swap(Node *n1, Node *n2);
    void replace(BST *tree);
    void delete_tree(Node *node);
    void copy(Node *ori_node, Node *new_parent, bool LR=false);
    void reset_num_nodes(){num_nodes = 0;}
    void inc_node(){num_nodes++;}
    void print(Node *node);
};

class Contour
{
    vector<XY> _hor; // horizontal contour
    vector<XY> _ver; // vertical contour
    int max_height;
public:
    Contour(){
        reset();
    }
    int get_max_height() const{return max_height;}
    int find_blk_y(Block *b, int x);
    void update(Block *b);
    void print(){
        cout<<"Horizontal contour:\n";
        for(XY c : _hor){
            cout<<"("<<c.first<<","<<c.second<<") ";
        }
        cout<<endl;
    }
    // sanity check
    void check();
    void reset(){
        _hor.clear();
        _ver.clear();
        _hor.push_back(XY(0,0));
        _hor.push_back(XY(INT_MAX,0));
        _ver.push_back(XY(0,0));
        _ver.push_back(XY(0,INT_MAX));
    }
};