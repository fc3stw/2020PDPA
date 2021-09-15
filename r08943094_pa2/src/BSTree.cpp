#include "BSTree.h"

using namespace std;

bool BST::remove(Node *node)
{
    if(node->left && node->right || node==get_root()) return false;
    Node *parent = node->parent;
    if(!node->left && !node->right){
        if(parent->left == node) parent->left = nullptr;
        else parent->right = nullptr;
    }
    else if(node->left){
        if(parent->left == node) parent->left = node->left;
        else if(parent->right == node) parent->right = node->left;
        else{
            cerr<<"Error: node parent does not point to node\n";
            exit(1);
        }
        node->left->parent = parent;
    }
    else{
        if(parent->left == node) parent->left = node->right;
        else if(parent->right == node) parent->right = node->right;
        else{
            cerr<<"Error: node parent does not point to node\n";
            exit(1);
        }
        node->right->parent = parent;
    }
    node->parent = nullptr;
    node->left = nullptr;
    node->right = nullptr;
    return true;
}

bool BST::append(Node *parent, Node *node, bool LR)
{
    if(parent==nullptr){
        set_root(node);
        return true;
    }
    if(LR){
        if(parent->right) return false;
        parent->right = node;
    }
    else{
        if(parent->left) return false;
        parent->left = node;
    }
    node->parent = parent;
    return true;
}

void BST::rotate(Node *node)
{
    node->orien = !node->orien;
}

void BST::swap(Node *n1, Node *n2)
{
    int tmp_id = n1->id;
    bool tmp_orien = n1->orien;
    n1->id = n2->id;
    n1->orien = n2->orien;
    n2->id = tmp_id;
    n2->orien = tmp_orien;
}

// replace an entire tree by another tree
void BST::replace(BST *tree)
{
    // delete all nodes
    delete_tree(get_root());
    // copy all nodes
    copy(tree->get_root(), nullptr);
    // annotate cost
    set_cost(tree->get_cost());
    // print(root);
}

void BST::delete_tree(Node *node)
{
    if(!node) return;
    delete_tree(node->left);
    delete_tree(node->right);
    delete node;
}

void BST::copy(Node *ori_node, Node *new_parent, bool LR)
{
    if(!ori_node) return;
    Node *new_node = new Node(ori_node);
    append(new_parent, new_node, LR);
    copy(ori_node->left, new_node, false);
    copy(ori_node->right, new_node, true);
}

void BST::print(Node *node)
{
    if(!node) return;
    cout<<"Node #"<<node->id<<":\n";
    if(node->parent) cout<<"parent= #"<<node->parent->id<<"\n";
    if(node->left) cout<<"left= #"<<node->left->id<<"\n";
    if(node->right) cout<<"right= #"<<node->right->id<<"\n";
    cout<<endl;
    print(node->left);
    print(node->right);
}

int Contour::find_blk_y(Block *b, int x)
{
    // block x has to be correct
    int blk_x1 = x;
    int blk_x2 = blk_x1 + b->get_width();
    // if there's a corner within the range of the block horizontally
    int max_y = 0;
    bool corner_in_hor_range = false;
    for(vector<XY>::iterator it = _hor.begin(); it!=_hor.end(); it++){
        XY corner = *it;
        if(corner.first < blk_x1) continue;
        if(corner.first > blk_x2) break;
        if(it!=_hor.end()-1 && corner.first==blk_x1 && corner.first==(it+1)->first && corner.second>(it+1)->second) continue;
        if(it!=_hor.begin() && corner.first==blk_x2 && corner.first==(it-1)->first && corner.second>(it-1)->second) break;
        max_y = max(max_y, corner.second);
        corner_in_hor_range = true;
    }
    // if there's no corner within the range of the block horizontally
    if(!corner_in_hor_range){
        for(XY corner : _hor){
            if(corner.first > blk_x1) break;
            max_y = corner.second;
        }
    }
    return max_y;
}

void Contour::update(Block *b)
{
    // new corners
    int x1, y1, x2, y2;
    tie(x1, y1) = b->get_low_xy();
    tie(x2, y2) = b->get_top_xy();
    // update horizontal corner
    // newly added horizontal corners
    vector<XY> new_corners;
    new_corners.emplace_back(x1,y1);
    new_corners.emplace_back(x1,y2);
    new_corners.emplace_back(x2,y2);
    new_corners.emplace_back(x2,y1);
    // merge with current corner list
    vector<XY>::iterator it = _hor.begin();
    while(it!=_hor.end()){
        XY xy = *it;
        // remove corners lying within the interval
        if(xy.first>x1 && xy.first<x2) it = _hor.erase(it)-1;
        if(xy.first>=x2) break;
        it++;
    }
    // insert an extra corner to connect with current contour
    if(it->second < new_corners.back().second){
        new_corners.emplace_back(x2,it->second);
    }
    _hor.insert(it, new_corners.begin(), new_corners.end());
    // remove redundant corners
    vector<XY>::iterator it_prev = _hor.begin();
    it = it_prev+1;
    while(it!=_hor.end()){
        XY xy = *it;
        XY xy_prev = *it_prev;
        if(xy.first != xy_prev.first){
            // if the track forms a vertical line, remove the points at the middle
            if(it - it_prev > 2){
                int dif = it-it_prev-2;
                it_prev++;
                while(dif--){
                    it_prev = _hor.erase(it_prev);
                }
                // if both end has the same points, remove both
                if(*it_prev==*(it_prev-1)){
                    it_prev = _hor.erase(it_prev-1);
                    it_prev = _hor.erase(it_prev);
                }
                it = it_prev;
            }
            else{
                it_prev = it;
            }
        }
        it++;
    }
    // update max height
    max_height = 0;
    for(XY c : _hor){
        max_height = max(max_height, c.second);
    }
    // WA: insert missing points
    it = _hor.begin();
    while(it != _hor.end()-1){
        if(it->first!=(it+1)->first && it->second!=(it+1)->second){
            if(it->second > (it+1)->second){
                it = _hor.insert(it+1, XY(it->first, (it+1)->second));
            }
            else{
                it = _hor.insert(it+1, XY((it+1)->first, it->second));
            }
        }
        it++;
    }
}

void Contour::check()
{
    // check horizontal
    vector<XY>::iterator it = _hor.begin()+1;
    while(it!=_hor.end()){
        XY xy = *it;
        if(xy.first < (it-1)->first){
            cerr<<"Error horizontal contour: next corner has larger x value than previous corner\n";
            print();
            exit(1);
        }
        if(xy.first != (it-1)->first && xy.second != (it-1)->second){
            cerr<<"Error horizontal contour: contour not continuous\n";
            print();
            exit(1);
        }
        it++;
    }
}