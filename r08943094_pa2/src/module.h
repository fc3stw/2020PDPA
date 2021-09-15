#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <tuple>

using namespace std;

typedef pair<int, int> XY;
typedef pair<double, double> CenterXY;

struct Node;

class Terminal
{
    int _id;
    string _name;
    int _x;
    int _y;
public:
    Terminal(int id, string name, int x, int y):
    _id(id), _name(name), _x(x), _y(y)
    {}
    // getter
    int get_id() const{return _id;}
    string get_name() const{return _name;}
    int get_x() const{return _x;}
    int get_y() const{return _y;}
    virtual CenterXY get_center_xy() const{return CenterXY(_x, _y);}
    // setter
    void set_x(int val){_x = val;}
    void set_y(int val){_y = val;}
    virtual void print(){
        cout<<"Terminal "<<_name<<"\t"
        <<"coordinate: ("<<_x<<","<<_y<<")\n";
    }
};

class Block : public Terminal
{
    Node *_node;
    int _width;
    int _height;
public:
    Block(int id, string name, int w, int h, XY low_xy=XY(0,0)):
    Terminal(id, name, low_xy.first, low_xy.second),
    _width(max(w,h)),
    _height(min(w,h))
    {}
    // getter
    Node* get_node() const{return _node;}
    int get_width() const{return _width;}
    int get_height() const{return _height;}
    bool get_orien() const{return _width < _height;}
    XY get_low_xy() const{return XY(get_x(), get_y());}
    XY get_top_xy() const{return XY(get_x()+_width, get_y()+_height);}
    CenterXY get_center_xy() const{return CenterXY((double)get_x()+(double)_width/2.0,(double)get_y()+(double)_height/2.0);}
    // setter
    void set_node(Node *node){_node = node;}
    // method
    void rotate(){int tmp = _width; _width = _height; _height = tmp;}
    void print(){
        cout<<"Block "<<get_name()<<" #"<<get_id()<<"\t"
        <<"dimension: ("<<_width<<","<<_height<<")\t"
        <<"lower-left: ("<<get_x()<<","<<get_y()<<")       \t"
        <<"upper-right: ("<<get_top_xy().first<<","<<get_top_xy().second<<")\n";
    }
};

class Net
{
    int _id;
    vector<Terminal*> _terminals;
public:
    Net(int id) : _id(id) {}
    // method
    void add_terminal(Terminal *term){_terminals.push_back(term);}
    // getter
    int get_id() const{return _id;}
    const vector<Terminal*>& get_terminals() const{return _terminals;}
};

struct Node
{
    int id; // the same as block id
    bool orien; // orientation, 0 for larger width, 1 for larger height
    Node *left;
    Node *right;
    Node *parent;
    Node(Block *block):
    id(block->get_id()),
    orien(block->get_orien()),
    left(nullptr),
    right(nullptr),
    parent(nullptr)
    {
        block->set_node(this);
    }
    Node(Node *node):
    id(node->id),
    orien(node->orien),
    left(nullptr),
    right(nullptr),
    parent(nullptr)
    {}
};
