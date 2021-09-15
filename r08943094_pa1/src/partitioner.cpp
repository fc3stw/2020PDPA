#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::partition()
{
    int verbosity = 0;
    int extra_iters = 0;
    if(_verbose>verbosity)  cout<<"Start Partitioning\n";
    start_timing();
    if(_verbose>verbosity)  cout<<"Looking for initial partition\n";
    // balance 2 partitions under constraint
    // find _maxPinNum for _bList size
    initialize_partitions();
    if(_verbose>verbosity)  cout<<"Initial partition found in "<<get_time()<<" sec\n";
    if(_verbose>verbosity)  cout<<"Estimate initial cut size\n";
    compute_net_part_count();
    if(_verbose>verbosity)  estimate_cut_size();
    if(_verbose>verbosity)  printSummary();
    if(_verbose>verbosity)  cout<<"Start optimization\n";
    // repeat FM algorithm until no further reduction in cutsize after the iteration is done
    _iterNum = 0;
    do{
        reset_all_parameters();
        // if(_verbose>verbosity)  cout<<"Parameters reset\n";
        compute_cell_gain();
        // if(_verbose>verbosity)  cout<<"Done computing initial cell gains\n";
        initialize_bucket_list();
        // if(_verbose>verbosity)  cout<<"Bucket list initialized\n";
        fm_partition_iteration();   // end if no unlocked cell left or all cells left cause unbalance (unmovable)
        // if(_verbose>verbosity)  cout<<"All movable cells moved\n";
        restore_best_move();
        // if(_verbose>verbosity)  cout<<"Best partition restored\n";
        // check_net_part_count();
        // if(_verbose>verbosity)  cout<<"net part count is correct\n";
        _iterNum++;
        if(_verbose>verbosity)  estimate_cut_size();
        if(_verbose>verbosity)  cout<<"Iteration "<<_iterNum<<": cut size = "<<_cutSize<<", max acc gain = "<<_maxAccGain<<" on move "<<_moveNum<<endl;
        if(_maxAccGain > 0){
            extra_iters = 0;
        }
        else{
            extra_iters++;
            // cout<<"extra iter = "<<extra_iters<<", max extra iter = "<<_max_extra_iters<<endl;
        }
    }while(extra_iters < _max_extra_iters);    // or maybe _bestMoveNum > 0?? No could lead to endless loop over the partitions without improvement
    estimate_cut_size();
    cout<<"Partitioning finished in "<<get_time()<<" sec\n";
}

void Partitioner::initialize_partitions()
{
    int half_cell_num = _cellNum / 2;
    _maxPinNum = 0;
    for(Cell *cell : _cellArray){
        if(cell->getNode()->getId() < half_cell_num){
            cell->move();
            _partSize[1]++;
        }
        else{
            _partSize[0]++;
        }
        // record max pin num for bucketlist
        _maxPinNum = max(_maxPinNum, cell->getPinNum());
    }
    _bList[0] = BucketList(_maxPinNum);
    _bList[1] = BucketList(_maxPinNum);
}

bool compare_cell_by_pins(const Cell *c1, const Cell *c2){return c1->getPinNum()<c2->getPinNum();}

// void Partitioner::initialize_partitions()
// {
//     // sort cells by num of pins
//     vector<Cell*> sort_cells_by_pins(_cellNum, nullptr);
//     copy(_cellArray.begin(), _cellArray.end(), sort_cells_by_pins.begin());
//     sort(sort_cells_by_pins.begin(), sort_cells_by_pins.end(), compare_cell_by_pins);
//     // move half of the cells of less pins to the other partition
//     int half_cell_num = _cellNum / 2;
//     for(int i = 0; i<half_cell_num; i++){
//         Cell *cell = sort_cells_by_pins[i];
//         cell->move();
//         _partSize[1]++;
//     }
//     _partSize[0] = _cellNum - _partSize[1];
//     // record max pin num for bucketlist
//     _maxPinNum = sort_cells_by_pins.back()->getPinNum();
//     _bList[0] = BucketList(_maxPinNum);
//     _bList[1] = BucketList(_maxPinNum);
// }

void Partitioner::compute_net_part_count()
{
    for(Cell *cell : _cellArray){
        int part = cell->getPart();
        for(int net_id : cell->getNetList()){
            _netArray[net_id]->incPartCount(part);
        }
    }
}

void Partitioner::reset_all_parameters()
{
    _accGain = 0;
    _maxAccGain = 0;
    _moveNum = 0;
    _bestMoveNum = 0;
    _unlockNum[0] = _partSize[0];
    _unlockNum[1] = _partSize[1];
    _moveStack.clear();
    // recalculate all net part count
    for(Net *net : _netArray){
        net->setPartCount(0, 0);
        net->setPartCount(1, 0);
        for(int cell_id : net->getCellList()){
            net->incPartCount(_cellArray[cell_id]->getPart());
        }
    }
}

void Partitioner::compute_cell_gain()
{
    // clear all cell gain
    for(Cell *cell : _cellArray){
        cell->setGain(0);
    }
    // compute initial cell gain
    for(Cell *cell : _cellArray){
        bool from = cell->getPart();
        bool to = !from;
        for(int net_id : cell->getNetList()){
            Net *net = _netArray[net_id];
            if(net->getPartCount(from)==1)  cell->incGain();
            if(net->getPartCount(to)==0) cell->decGain();
        }
    }
}

void Partitioner::initialize_bucket_list()
{
    _bList[0].clear();
    _bList[1].clear();
    for(Cell *cell : _cellArray){
        _bList[cell->getPart()].append(cell->getNode(), cell->getGain());
        assert((_cellArray[cell->getNode()->getId()]==cell) && "CELL ID MISMATCH CELL ARRAY INDEX");
        assert((_cellArray[_cellName2Id[cell->getName()]]==cell) && "CELL ID MISMATCH CELL ARRAY INDEX");
    }
}

void Partitioner::fm_partition_iteration()
{
    int verbosity = 1;
    int BC = _cellNum;  // balance coefficient
    // continue the iteration if there are still cells unlocked
    while(_bList[0].get_size() > 0 || _bList[1].get_size() > 0){
        // get cell with max gain under balance constraint
        // return nullptr if no cell to choose
        if(_verbose>verbosity)  cout<<"Update max cell gain\n";
        update_max_cell_gain();
        if(_verbose>verbosity)  cout<<"Looking for max gain cell\n";
        Cell *cell = get_balanced_max_gain_cell();
        if(cell==nullptr){
            // end the iteration if no cell to choose
            break;
        }
        if(_verbose>verbosity)  cout<<"Max gain cell under balance constraint acquired\n";
        // update step info
        move_cell(cell);
        lock_cell(cell);
        _accGain += cell->getGain();
        _moveNum++;
        int current_BC = abs(_partSize[0] - _partSize[1]);
        // if(_accGain > _maxAccGain){
        if(_accGain >= _maxAccGain && _moveNum!=_cellNum){
            BC = current_BC;
            _maxAccGain = _accGain;
            _bestMoveNum = _moveNum;
            if(_verbose>verbosity)  cout<<"move "<<_moveNum<<": base cell (id, gain) = ("<<cell->getNode()->getId()<<","<<cell->getGain()<<"), acc gain = "<<_accGain<<endl;
        }
        // else if(_accGain == _maxAccGain && current_BC<BC){
        //     BC = current_BC;
        //     _maxAccGain = _accGain;
        //     _bestMoveNum = _moveNum;
        //     if(_verbose>verbosity)  cout<<"move "<<_moveNum<<": base cell (id, gain) = ("<<cell->getNode()->getId()<<","<<cell->getGain()<<"), acc gain = "<<_accGain<<endl;
        // }
        if(_verbose>verbosity)  cout<<"Updating cell gain\n";
        update_gain(cell);
        if(_verbose>verbosity)  cout<<"Cell gain updated\n";
    }
}

void Partitioner::update_max_cell_gain()
{
    _maxCellGain[0] = _maxPinNum;
    while(_maxCellGain[0]>-_maxPinNum && !_bList[0].find_gain(_maxCellGain[0])){
        _maxCellGain[0]--;
    }
    _maxCellGain[1] = _maxPinNum;
    while(_maxCellGain[1]>-_maxPinNum && !_bList[1].find_gain(_maxCellGain[1])){
        _maxCellGain[1]--;
    }
}

Cell* Partitioner::get_balanced_max_gain_cell()
{
    int part = -1;
    Node *n[2];
    n[0] = _bList[0].get_node(_maxCellGain[0]);
    n[1] = _bList[1].get_node(_maxCellGain[1]);
    // check node existence
    if(n[0]==nullptr && n[1]==nullptr)  return nullptr;
    else if(n[0]==nullptr)  return _cellArray[n[1]->getId()];
    else if(n[1]==nullptr)  return _cellArray[n[0]->getId()];
    // check balance condition
    bool balanced0 = check_balance(_cellArray[n[0]->getId()]);
    bool balanced1 = check_balance(_cellArray[n[1]->getId()]);
    // decide the legal or better cell to move
    if(!balanced0 && !balanced1){
        assert(false && "Balancing error: cannot meet balance constraint");
    }
    else if(balanced0 && balanced1){
        part = (_maxCellGain[0]>=_maxCellGain[1])? 0 : 1;
        // if(_partSize[0]>_partSize[1])   part = 0;
        // else    part = 1;
    }
    else{
        part = balanced0? 0 : 1;
    }
    assert((part==0 || part==1) && "Error: cannot decide part while getting max gain cell");
    assert((n[part]!=nullptr) && "Error: acquire null cell");
    int cell_id = n[part]->getId();
    Cell *max_gain_cell = _cellArray[cell_id];
    return max_gain_cell;
}

bool Partitioner::check_balance(Cell *cell)
{
    bool part = cell->getPart();
    double ub = static_cast<double>(_cellNum)*(1. + _bFactor)/2.;
    double lb = static_cast<double>(_cellNum)*(1. - _bFactor)/2.;
    double new_from = _partSize[part]-1;
    double new_to = _partSize[!part]+1;
    return new_from > lb && new_to < ub;
}

void Partitioner::update_gain(Cell* cell)
{
    bool to_part = cell->getPart();
    bool from_part = !to_part;
    for(int net_id : cell->getNetList()){
        Net *net = _netArray[net_id];
        int from_size = net->getPartCount(from_part)+1;
        int to_size = net->getPartCount(to_part)-1;
        int only_from_cell = -1;
        int only_to_cell = -1;
        // before the move
        // cout<<"before move, to size = "<<to_size<<endl;
        if(to_size==0){
            // increment all the other cells
            for(int cell_id : net->getCellList()){
                Cell *c = _cellArray[cell_id];
                if(c->getLock())  continue;
                _bList[c->getPart()].remove(c->getNode(), c->getGain());
                c->incGain();
                _bList[c->getPart()].append(c->getNode(), c->getGain());
            }
        }
        else if(to_size==1){
            for(int cell_id : net->getCellList()){
                if(_cellArray[cell_id]->getPart()==to_part && !_cellArray[cell_id]->getLock()){
                    only_to_cell = cell_id;
                    break;
                }
            }
            if(only_to_cell!=-1){
                Cell *c = _cellArray[only_to_cell];
                _bList[c->getPart()].remove(c->getNode(), c->getGain());
                c->decGain();
                _bList[c->getPart()].append(c->getNode(), c->getGain());
            }
        }
        from_size--;
        to_size++;
        // after the move
        // cout<<"after move, from size = "<<from_size<<endl;
        if(from_size==0){
            // decrement all the other cells
            for(int cell_id : net->getCellList()){
                Cell *c = _cellArray[cell_id];
                if(c->getLock())  continue;
                _bList[c->getPart()].remove(c->getNode(), c->getGain());
                c->decGain();
                _bList[c->getPart()].append(c->getNode(), c->getGain());
            }
        }
        else if(from_size==1){
            for(int cell_id : net->getCellList()){
                if(_cellArray[cell_id]->getPart()==from_part && !_cellArray[cell_id]->getLock()){
                    only_from_cell = cell_id;
                    break;
                }
            }
            if(only_from_cell!=-1){
                Cell *c = _cellArray[only_from_cell];
                _bList[c->getPart()].remove(c->getNode(), c->getGain());
                c->incGain();
                _bList[c->getPart()].append(c->getNode(), c->getGain());
            }
        }
    }
}

void Partitioner::restore_best_move()
{
    // unlock all cells
    for(Cell *cell : _cellArray){
        cell->unlock();
    }
    // restore to the best move
    while(_moveNum != _bestMoveNum){
        Cell *cell = _cellArray[_moveStack.back()];
        _moveStack.pop_back();
        move_cell(cell, true);
        _moveNum--;
    }
}

void Partitioner::move_cell(Cell *cell, bool reverse)
{
    assert(!cell->getLock() && "cannot move locked cell");
    _unlockNum[cell->getPart()]--;
    _partSize[cell->getPart()]--;
    cell->move();
    _partSize[cell->getPart()]++;
    if(!reverse)    _moveStack.push_back(cell->getNode()->getId());
    bool to = cell->getPart();
    bool from = !to;
    for(int net_id : cell->getNetList()){
        Net *net = _netArray[net_id];
        net->decPartCount(from);
        net->incPartCount(to);
    }
}

void Partitioner::lock_cell(Cell *cell)
{
    cell->lock();
    bool part = !cell->getPart();
    _bList[part].remove(cell->getNode(), cell->getGain());
}

void Partitioner::estimate_cut_size()
{
    _cutSize = 0;
    for(Net *net : _netArray){
        if(net->getPartCount(0)>0 && net->getPartCount(1)>0) _cutSize++;
    }
}

void Partitioner::check_net_part_count()
{
    for(Net *net : _netArray){
        int part_cnt[2] = {};
        for(int cell_id : net->getCellList()){
            Cell *cell = _cellArray[cell_id];
            part_cnt[cell->getPart()]++;
        }
        assert((part_cnt[0]==net->getPartCount(0) && part_cnt[1]==net->getPartCount(1)) && "net part count mismatch");
    }
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
