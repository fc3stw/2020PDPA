#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>
#include <cassert>
#include "cell.h"
#include "net.h"
using namespace std;

#define VERBOSE 0
#define MAX_EXTRA_ITERS 5

class BucketList{
    vector<Node*> blist;
    size_t size;
    size_t max_size;
    size_t offset;
    bool find_gain_by_idx(int idx){return blist[idx]!=nullptr;}
public:
    BucketList(){
        blist.clear();
        offset = 0;
        max_size = 0;
        offset = 0;
    }
    BucketList(int pmax){
        blist.clear();
        offset = pmax;
        max_size = 2*pmax+1;
        blist = vector<Node*>(max_size, nullptr);
    }
    bool find_gain(int g){return blist[g+offset]!=nullptr;}
    size_t get_size() const {return size;}
    Node* get_node(int g) const {
        int idx = g + offset;
        return blist[idx];
    }
    void append(Node* n, int g){
        int idx = g + offset;
        if(find_gain_by_idx(idx)){
            blist[idx]->setNext(n);
            n->setPrev(blist[idx]);
            blist[idx] = n;
        }
        else{
            blist[idx] = n;
        }
        size++;
    }
    void remove(Node* n, int g){
        int idx = g + offset;
        assert(find_gain_by_idx(idx) && "ERROR: cannot remove cell not in the bucketlist");
        // only one node in the list
        if(n->getPrev()==nullptr && n->getNext()==nullptr){
            blist[idx] = nullptr;
        }
        // remove the last node in the list
        else if(n->getNext()==nullptr){
            blist[idx] = n->getPrev();
            n->setPrev(nullptr);
            blist[idx]->setNext(nullptr);
        }
        // remove the first node
        else if(n->getPrev()==nullptr){
            n->getNext()->setPrev(nullptr);
            n->setNext(nullptr);
        }
        else{
            n->getNext()->setPrev(n->getPrev());
            n->getPrev()->setNext(n->getNext());
            n->setNext(nullptr);
            n->setPrev(nullptr);
        }
        size--;
    }
    void clear(){
        for(int i = 0; i<max_size; i++){
            blist[i] = nullptr;
        }
        size = 0;
    }
};

class Partitioner
{
public:
    // constructor and destructor
    Partitioner(fstream& inFile) :
        _cutSize(0), _netNum(0), _cellNum(0), _maxPinNum(0), _bFactor(0),
        _accGain(0), _maxAccGain(0), _iterNum(0) {
        parseInput(inFile);
        _partSize[0] = 0;
        _partSize[1] = 0;
        _verbose = VERBOSE;
        _max_extra_iters = MAX_EXTRA_ITERS;
    }
    ~Partitioner() {
        clear();
    }

    // basic access methods
    int getCutSize() const          { return _cutSize; }
    int getNetNum() const           { return _netNum; }
    int getCellNum() const          { return _cellNum; }
    double getBFactor() const       { return _bFactor; }
    int getPartSize(int part) const { return _partSize[part]; }

    // modify method
    void parseInput(fstream& inFile);
    void partition();

    // member functions about reporting
    void printSummary() const;
    void reportNet() const;
    void reportCell() const;
    void writeResult(fstream& outFile);

private:
    int                 _cutSize;       // cut size
    int                 _partSize[2];   // size (cell number) of partition A(0) and B(1)
    int                 _netNum;        // number of nets
    int                 _cellNum;       // number of cells
    int                 _maxPinNum;     // Pmax for building bucket list
    double              _bFactor;       // the balance factor to be met
    // Node*               _maxGainCell;   // pointer to max gain cell
    int                 _maxCellGain[2];   // entry for choosing cell with max gain
    vector<Net*>        _netArray;      // net array of the circuit
    vector<Cell*>       _cellArray;     // cell array of the circuit
    // map<int, Node*>     _bList[2];      // bucket list of partition A(0) and B(1)
    BucketList          _bList[2];      // bucket list
    map<string, int>    _netName2Id;    // mapping from net name to id
    map<string, int>    _cellName2Id;   // mapping from cell name to id

    // parameters that need to be reset for each fm iteration
    int                 _accGain;       // accumulative gain
    int                 _maxAccGain;    // maximum accumulative gain
    int                 _moveNum;       // number of cell movements
    int                 _iterNum;       // number of iterations
    int                 _bestMoveNum;   // store best number of movements
    int                 _unlockNum[2];  // number of unlocked cells
    vector<int>         _moveStack;     // history of cell movement
    int                 _max_extra_iters;

    int                 _verbose;       // 0 to print nothing, 1 for each iteration, 2 for each move
    clock_t             _start_time;

    // Clean up partitioner
    void clear();

    // PA1 add
    void initialize_partitions();
    void compute_net_part_count();
    void reset_all_parameters();
    void compute_cell_gain();
    void initialize_bucket_list();
    void fm_partition_iteration();
    void update_max_cell_gain();
    Cell* get_balanced_max_gain_cell();
    bool check_balance(Cell*);
    void update_gain(Cell*);
    void restore_best_move();
    void move_cell(Cell *cell, bool reverse=false);
    void lock_cell(Cell *cell);
    void estimate_cut_size();

    // sanity checks
    void check_net_part_count();
    void start_timing(){_start_time = clock();}
    double get_time() const {return (double)(clock() - _start_time) / CLOCKS_PER_SEC;}
};

#endif  // PARTITIONER_H
