#include "floorplanner.h"

#define RAND unif(rand_gen)
#define RGB(r,g,b) (65536 * int(r) + 256 * int(g) + int(b))

using namespace std;

vector<int> Floorplanner::get_highest_placed_block_ids() const
{
    // traverse the b*-tree with BFS
    vector<int> block_ids;
    queue<Node*> q;
    q.push(current_bst->get_root());
    Node *node;
    while(!q.empty()){
        node = q.front();
        q.pop();
        int id = node->id;
        Block *block = block_list[id];
        // assumption: the blocks must not have right child
        if(block->get_top_xy().second==get_chip_height() && !node->right) block_ids.push_back(id);
        if(node->left)  q.push(node->left); 
        if(node->right)  q.push(node->right);
    }
    return block_ids;
}

void Floorplanner::place_block(Node *parent, Block *block, bool LR)
{
    int x,y;
    if(parent==nullptr){
        x = 0;
        y = 0;
        // cout<<"Place base block\n";
    }
    else if(LR){
        x = block_list[parent->id]->get_low_xy().first;
        y = contour.find_blk_y(block, x);
        // cout<<"Place block on the top\n";
    }
    else{
        x = block_list[parent->id]->get_top_xy().first;
        y = contour.find_blk_y(block, x);
        // cout<<"Place block on the right\n";
    }
    block->set_x(x);
    block->set_y(y);
    // block->print();
    contour.update(block);
    // contour.print();
    // cout<<endl;
}

// compute with in-order traversal for a compacted result
void Floorplanner::update_all_blocks(bool sampling_mode)
{
    contour.reset();
    current_bst->reset_num_nodes();
    chip_width = 0;
    chip_height = 0;
    // update coordinates for all blocks and maintain contour given current B*-tree topology
    // update chip width and height
    update_block_dfs(current_bst->get_root());
    contour.check();
    // update wirelength
    calculate_wl();
    // sample the result and update normalized cost
    if(sampling_mode) update_normalized_cost();
    // annotate normalized cost on B*-tree
    current_bst->set_cost(get_normalized_cost());
    if(current_bst->get_num_nodes() != num_blocks){
        cerr<<"Block number mismatch with node number\n";
        cerr<<"num nodes: "<<current_bst->get_num_nodes()<<endl;
        cerr<<"num blocks: "<<num_blocks<<endl;
        exit(1);
    }
}

void Floorplanner::update_block_dfs(Node *node)
{
    if(!node) return;
    current_bst->inc_node();
    Block *block = block_list[node->id];
    // check block orientation
    if(block->get_orien()!=node->orien) block->rotate();
    // set block coordinate
    if(node==current_bst->get_root()){
        place_block(nullptr, block);
    }
    else{
        // find parent and check LR
        Node *parent = node->parent;
        bool LR = (parent->right==node);
        place_block(parent, block, LR);
    }
    // update chip width and height
    int topx, topy;
    tie(topx, topy) = block->get_top_xy();
    chip_width = max(chip_width, topx);
    chip_height = max(chip_height, topy);
    if(node->left==node || node->right==node){
        cerr<<"Error: Node points to itself\n";
        exit(1);
    }
    update_block_dfs(node->left);
    update_block_dfs(node->right);
}

void Floorplanner::calculate_wl()
{
    WL = 0.;
    for(Net *net : net_list){
        double x1 = INT_MAX;
        double x2 = 0;
        double y1 = INT_MAX;
        double y2 = 0;
        for(Terminal *term : net->get_terminals()){
            double termx, termy;
            tie(termx, termy) = term->get_center_xy();
            x1 = min(x1, termx);
            x2 = max(x2, termx);
            y1 = min(y1, termy);
            y2 = max(y2, termy);
        }
        double wl = x2-x1 + y2-y1;
        WL += wl;
    }
}

void Floorplanner::backup_floorplan(bool force)
{
    if(force || is_legal() && current_bst->get_cost() < best_bst->get_cost()){
        best_bst->replace(current_bst);
        print_summary(0);
    }
}

void Floorplanner::restore_floorplan(BST *tree)
{
    current_bst->replace(tree);
    // restore link of node on blocks
    restore_block_node(current_bst->get_root());
    update_all_blocks();
}

void Floorplanner::restore_block_node(Node *node)
{
    if(!node) return;
    block_list[node->id]->set_node(node);
    restore_block_node(node->left);
    restore_block_node(node->right);
}

void Floorplanner::random_operation()
{
    int op = RAND*3;
    int b1_id = RAND*num_blocks;
    int b2_id = RAND*num_blocks;
    bool LR = RAND > 0.5;
    while(b2_id == b1_id){b2_id = RAND*num_blocks;}
    if(op==0){
        // remove b1
        Node *n1 = block_list[b1_id]->get_node();
        while(!current_bst->remove(n1)){
            b1_id = RAND*num_blocks;
            n1 = block_list[b1_id]->get_node();
        }
        // append b1 under b2
        Node *n2 = block_list[b2_id]->get_node();
        while(b1_id==b2_id || !current_bst->append(n2, n1, LR)){
            LR = RAND > 0.5;
            b2_id = RAND*num_blocks;
            n2 = block_list[b2_id]->get_node();
        }
    }
    else if(op==1){
        // rotate b1
        Node *n1 = block_list[b1_id]->get_node();
        current_bst->rotate(n1);
    }
    else if(op==2){
        // swap b1 and b2
        Node *n1 = block_list[b1_id]->get_node();
        Node *n2 = block_list[b2_id]->get_node();
        current_bst->swap(n1, n2);
    }
    // cout<<"op="<<op<<", b1="<<block_list[b1_id]->get_name()<<", b2="<<block_list[b2_id]->get_name()<<", LR="<<LR<<endl;
}

void Floorplanner::sample_normalized_cost()
{
    for(int i = 0; i<num_blocks*10; i++){
        random_operation();
        update_all_blocks(true);
    }
}

void Floorplanner::update_normalized_cost()
{
    norm_area *= num_samples;
    norm_wl *= num_samples;
    norm_area += get_area();
    norm_wl += get_wl();
    num_samples++;
    norm_area /= num_samples;
    norm_wl /= num_samples;
    // cout<<"sample #"<<num_samples<<", norm area = "<<norm_area<<", norm wl = "<<norm_wl<<endl;
}

bool compare_block_by_width(Block* b1, Block *b2)
{
    if(b1->get_width() > b2->get_width())   return true;
    else if(b1->get_width() < b2->get_width())   return false;
    else{
        if(b1->get_height() > b2->get_height())   return true;
        else return false;
    }
}

void Floorplanner::initial_floorplan(BlockList &unplaced_blocks)
{
    // all blocks have the same orientation
    // sort block by dominance relation
    sort(unplaced_blocks.begin(), unplaced_blocks.end(), compare_block_by_width);
    Block *block;
    Node *node;
    Node *rightmost;
    while(!unplaced_blocks.empty()){
        BlockList::iterator blkit = unplaced_blocks.begin();
        // start a new layer
        block = *blkit;
        node = new Node(block);
        if(current_bst->get_root()==nullptr){
            if(place_legal(nullptr, block, true)){
                place_block(nullptr, block);
                if(!current_bst->append(nullptr, block->get_node())){
                    cerr<<"Error: cannot append node\n";
                    exit(1);
                }
                rightmost = current_bst->get_rightmost();
                blkit = unplaced_blocks.erase(blkit);
            }
            else{
                block->rotate();
                place_block(nullptr, block);
                if(!current_bst->append(nullptr, block->get_node())){
                    cerr<<"Error: cannot append node\n";
                    exit(1);
                }
                rightmost = current_bst->get_rightmost();
                blkit = unplaced_blocks.erase(blkit);
            }
        }
        else{
            if(place_legal(rightmost, block, true)){
                place_block(rightmost, block, true);
                if(!current_bst->append(rightmost, block->get_node(), true)){
                    cerr<<"Error: cannot append node\n";
                    exit(1);
                }
                rightmost = node;
                blkit = unplaced_blocks.erase(blkit);
            }
            else{
                block->rotate();
                place_block(rightmost, block, true);
                if(!current_bst->append(rightmost, block->get_node(), true)){
                    cerr<<"Error: cannot append node\n";
                    exit(1);
                }
                rightmost = node;
                blkit = unplaced_blocks.erase(blkit);
            }
        }
        
        // place as many blocks as possible horizontally
        while(blkit!=unplaced_blocks.end()){
            block = *blkit;
            block->rotate();
            int block_y = contour.find_blk_y(block, block_list[node->id]->get_top_xy().first);
            if(place_legal(node, block, false) && contour.get_max_height()>=block_y+block->get_height()){
                // node->left = new Node(block);
                new Node(block);
                place_block(node, block, false);
                if(!current_bst->append(node, block->get_node(), false)){
                    cerr<<"Error: cannot append node\n";
                    exit(1);
                }
                node = node->left;
                blkit = unplaced_blocks.erase(blkit);
                continue;
            }
            else{
                // rotate back to keep dominance relation if still fail to place
                block->rotate();
                if(place_legal(node, block, false)){
                    // node->left = new Node(block);
                    new Node(block);
                    place_block(node, block, false);
                    if(!current_bst->append(node, block->get_node(), false)){
                        cerr<<"Error: cannot append node\n";
                        exit(1);
                    }
                    node = node->left;
                    blkit = unplaced_blocks.erase(blkit);
                    continue;
                }
            }
            blkit++;
        }
        // no more block can fill in the remaining horizontal space
        // the result would be a like merging several skewed trees
    }
    update_all_blocks(true);
}

// given an initial solution
void Floorplanner::SA()
{
    cout<<"Start SA optimization\n";
    int num_iter = 1;
    int num_operation = num_blocks*200;
    double T = 1.;
    double Tmin = 0.01;
    double r = 0.999;
    double delta;
    double prob;
    while(T > Tmin){
        // cout<<"Iter #"<<num_iter<<endl;
        restore_floorplan(best_bst);
        for(int i = 0; i<num_operation; i++){
            // cout<<"Op #"<<i+1<<"\r"<<flush;
            // backup current floorplan
            previous_bst->replace(current_bst);
            // randomly do an operation
            random_operation();
            update_all_blocks();
            delta = current_bst->get_cost() - previous_bst->get_cost();
            prob = min(1., exp(-delta/T));
            // restore undo the operation by probability
            if(RAND > prob){
                restore_floorplan(previous_bst);
            }
            // saved floorplan if it's the best
            backup_floorplan();
        }
        // cout<<endl;
        num_iter++;
        // reduce temperature
        T *= r;
    }
}

void Floorplanner::print_placed_blocks()
{
    cout<<"All placed blocks:\n";
    Node *node = current_bst->get_root();
    print_bst_dfs(node);
}

void Floorplanner::print_bst_dfs(Node *node)
{
    block_list[node->id]->print();
    // if(node->parent) cout<<"parent: "<<block_list[node->parent->id]->get_name()<<endl;
    if(node->left) print_bst_dfs(node->left);
    if(node->right) print_bst_dfs(node->right);
}

bool Floorplanner::check_outline()
{
    vector<int> ids;
    for(Block *block : block_list){
        int x2, y2;
        tie(x2, y2) = block->get_top_xy();
        if(x2 > outline_width || y2 > outline_height) ids.push_back(block->get_id());
    }
    for(int id : ids){
        cout<<"Block placed out of bounds\n";
        block_list[id]->print();
    }
    return ids.empty();
}

void Floorplanner::check_overlap()
{
    for(Block *b1 : block_list){
        for(Block *b2 : block_list){
            if(b1==b2) continue;
            // overlap if one of the corner of b1 lies within b2
            if(block_overlap(b1, b2))
            {
                cerr<<"Error: the following 2 blocks are overlapping\n";
                b1->print();
                b2->print();
            }
        }
    }
}

bool Floorplanner::block_overlap(Block *b1, Block *b2)
{
    // overlap if one of the corner of b1 lies within b2 or vise versa
    int b1x1, b1y1;
    int b1x2, b1y2;
    int b2x1, b2y1;
    int b2x2, b2y2;
    tie(b1x1, b1y1) =  b1->get_low_xy();
    tie(b1x2, b1y2) =  b1->get_top_xy();
    tie(b2x1, b2y1) =  b2->get_low_xy();
    tie(b2x2, b2y2) =  b2->get_top_xy();

    bool hor_overlap_1 = (b1x1 >= b2x1 && b1x1 < b2x2 || b1x2 > b2x1 && b1x2 <= b2x2);
    bool ver_overlap_1 = (b1y1 >= b2y1 && b1y1 < b2y2 || b1y2 > b2y1 && b1y2 <= b2y2);
    bool hor_overlap_2 = (b2x1 >= b1x1 && b2x1 < b1x2 || b2x2 > b1x1 && b2x2 <= b1x2);
    bool ver_overlap_2 = (b2y1 >= b1y1 && b2y1 < b1y2 || b2y2 > b1y1 && b2y2 <= b1y2);
    if(hor_overlap_1 && ver_overlap_1) return true;
    if(hor_overlap_2 && ver_overlap_2) return true;
    return false;
}

void Floorplanner::parse_block(fstream &block_file)
{
    int block_id = 0;
    int terminal_id = 0;
    string sbuf;
    while(block_file>>sbuf){
        if(sbuf == "Outline:"){
            string outline_w_in, outline_h_in;
            block_file>>outline_w_in;
            block_file>>outline_h_in;
            outline_width = stoi(outline_w_in);
            outline_height = stoi(outline_h_in);
        }
        else if(sbuf == "NumBlocks:"){
            block_file>>sbuf;
            num_blocks = stoi(sbuf);
        }
        else if(sbuf == "NumTerminals:"){
            block_file>>sbuf;
            num_terminals = stoi(sbuf);
        }
        else{
            string block_width;
            block_file>>block_width;
            if(block_width == "terminal"){
                if(terminalname2id.find(sbuf)!=terminalname2id.end()){
                    cerr<<"Terminal "<<sbuf<<" duplicates in .block file\n";
                    exit(1);
                }
                terminalname2id[sbuf] = terminal_id;
                string terminal_x, terminal_y;
                block_file>>terminal_x;
                block_file>>terminal_y;
                terminal_list.push_back(new Terminal(terminal_id, sbuf, stoi(terminal_x), stoi(terminal_y)));
                terminal_id++;
            }
            else{
                if(blkname2id.find(sbuf)!=blkname2id.end()){
                    cerr<<"Block "<<sbuf<<" duplicates in .block file\n";
                    exit(1);
                }
                blkname2id[sbuf] = block_id;
                string block_height;
                block_file>>block_height;
                block_list.push_back(new Block(block_id, sbuf, stoi(block_width), stoi(block_height)));
                block_id++;
            }
        }
    }
    if(block_id!=num_blocks || block_id!=block_list.size()){
        cerr<<"Wrong block number in .block file\n";
        exit(1);
    }
    if(terminal_id!=num_terminals || terminal_id!=terminal_list.size()){
        cerr<<"Wrong terminal number in .block file\n";
        exit(1);
    }
    cout<<"Outline dimension = ("<<outline_width<<","<<outline_height<<")\n";
    cout<<"Total "<<num_blocks<<" blocks and "<<num_terminals<<" terminals added successfully\n";
}

void Floorplanner::parse_net(fstream &net_file)
{
    int net_id = 0;
    string sbuf;
    while(net_file>>sbuf){
        if(sbuf == "NumNets:"){
            net_file>>sbuf;
            num_nets = stoi(sbuf);
        }
        else if(sbuf == "NetDegree:"){
            net_file>>sbuf;
            int degree = stoi(sbuf);
            Net *net = new Net(net_id);
            for(int i = 0; i<degree; i++){
                net_file>>sbuf;
                Name2Id::iterator blk_it = blkname2id.find(sbuf);
                Name2Id::iterator terminal_it = terminalname2id.find(sbuf);
                bool has_block = blk_it!=blkname2id.end();
                bool has_terminal = terminal_it!=terminalname2id.end();
                if(!has_block && !has_terminal){
                    cerr<<"Undifined block or terminal "<<sbuf<<" in .nets file\n";
                    exit(1);
                }
                else if(has_block){
                    net->add_terminal(block_list[blk_it->second]);
                }
                else if(has_terminal){
                    net->add_terminal(terminal_list[terminal_it->second]);
                }
            }
            net_list.push_back(net);
            net_id++;
        }
    }
    if(net_id!=num_nets || net_id!=net_list.size()){
        cerr<<"Wrong net number in .nets file\n";
        exit(1);
    }
    cout<<"Total "<<num_nets<<" nets added successfully\n";
    cout<<endl;
}

void Floorplanner::floorplan()
{
    BlockList unplaced_blocks(block_list.begin(), block_list.end());
    initial_floorplan(unplaced_blocks);
    cout<<"Initial floorplan:\n";
    backup_floorplan(true);
    sample_normalized_cost();
    restore_floorplan(best_bst);
    backup_floorplan(true);
    SA();
    restore_floorplan(best_bst);
    check_outline();
    check_overlap();
}

void Floorplanner::plot(string file_name)
{
	string plot_dir = "plot";
    string cmd;
    bool ret;
    cmd = "mkdir -p gnu";
    ret = system(cmd.c_str());
    cmd = "mkdir -p "+plot_dir;
    ret = system(cmd.c_str());
	string data_path = string("gnu/data")+string(".txt");
	string gnu_path = string("gnu/gnu");
    ofstream data_file(data_path);
    ofstream gnu_file(gnu_path);
    
    // plot blocks
    int i = 0;
    for(auto block: block_list){
        int x1, y1, x2, y2;
        tie(x1, y1) = block->get_low_xy();
        tie(x2, y2) = block->get_top_xy();
        data_file<<x1<<" "<<y1<<endl;
        data_file<<x2<<" "<<y2<<endl;
        gnu_file<<"set object " << i+1 << " rect from "<<x1<<","<<y1<<" to "<<x2<<","<<y2<<" fillcolor lt 2 linewidth 3"<<endl;
        gnu_file<<"set label \""<<block->get_name() <<"\" at "<<(x1+x2)/2<<","<<(y1+y2)/2<<" front center font \",40\""<<endl;
        i++;
    }
    // plot outline
    data_file<<outline_width<<" "<<outline_height<<endl;
    gnu_file<<"set arrow from 0,"<<outline_height<<" to "<<outline_width<<","<<outline_height<<" nohead lc 3 lw 5"<<endl;
    gnu_file<<"set arrow from "<<outline_width<<",0 to "<<outline_width<<","<<outline_height<<" nohead lc 3 lw 5"<<endl;
    gnu_file<<"set term png size 2000,2000"<<endl;
    gnu_file<<"set output \'"<<plot_dir<<"/"<<file_name<<"\'"<<endl;
    gnu_file<<"plot \'"<< data_path <<"\' using 1:2 with points"<<endl;
    data_file.close();
    gnu_file.close();
    string command ="gnuplot "+gnu_path; 
    
    ret = system(command.c_str());
    cout<<"Plot saved in "<<plot_dir<<"/"<<file_name<<endl;
}

void Floorplanner::write_result(fstream &output)
{
    output<<fixed;
    // final cost
    output<<get_cost()<<"\n";
    // total wirelength
    output<<get_wl()<<"\n";
    // chip area
    output<<get_area()<<"\n";
    output<<chip_width<<" "<<chip_height<<"\n";
    // runtime
    output<<get_time()<<"\n";
    // report macros
    for(Block *block : block_list){
        int x1,x2,y1,y2;
        tie(x1, y1) = block->get_low_xy();
        tie(x2, y2) = block->get_top_xy();
        output<<block->get_name()<<" "<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<"\n";
    }
}

void Floorplanner::print_summary(int verbose)
{
    if(verbose > 0){
        cout<<"\n###################### Floorplan Summary ######################\n";
        cout<<"Chip dimension: ("<<chip_width<<","<<chip_height<<")\n";
        cout<<"Chip area = "<<get_area()<<"\n";
        cout<<"Total wire length = "<<get_wl()<<"\n";
        cout<<"Normalized cost = "<<get_normalized_cost()<<"\n";
        print_placed_blocks();
        cout<<"Program runtime: "<<get_time()<<" sec\n";
        cout<<endl;
    }
    else{
        cout<<"Chip dimension: ("<<chip_width<<","<<chip_height<<")\t";
        cout<<"Chip area = "<<get_area()<<"\t";
        cout<<"Total wire length = "<<get_wl()<<"\t";
        cout<<"Normalized cost = "<<get_normalized_cost()<<"\t";
        cout<<"Program runtime: "<<get_time()<<" sec\n"<<flush;
    }
}

void Floorplanner::clear()
{
    for (size_t i = 0, end = block_list.size(); i < end; ++i) {
        delete block_list[i];
    }
    for (size_t i = 0, end = terminal_list.size(); i < end; ++i) {
        delete terminal_list[i];
    }
    for (size_t i = 0, end = net_list.size(); i < end; ++i) {
        delete net_list[i];
    }
    return;
}
