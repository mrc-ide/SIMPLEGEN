
#include "relatedness_Block2.h"
#include "misc_v14.h"

using namespace std;

//------------------------------------------------
// constructor
relatedness_Block2::relatedness_Block2(int haplo_ID, int contig, int generation,
                                       int left, int right, int time, vector<int> &descendant_IDs,
                                       map<pair<int, int>, relatedness_Block2>* block_map_ptr,
                                       std::vector<relatedness_Coaltracker>* coaltracker_ptr,
                                       Rcpp::List* store_output_ptr) {
  
  // copy over objects
  this->haplo_ID = haplo_ID;
  this->contig = contig;
  this->generation = generation;
  this->block_map_ptr = block_map_ptr;
  this->coaltracker_ptr = coaltracker_ptr;
  this->store_output_ptr = store_output_ptr;
  
  // add first block
  add_block(left, right, time, descendant_IDs);
}

//------------------------------------------------
// add a new block
void relatedness_Block2::add_block(int left, int right, int time, vector<int> &descendant_IDs) {
  this->left.push_back(left);
  this->right.push_back(right);
  this->time.push_back(time);
  this->descendant_IDs.push_back(descendant_IDs);
}

//------------------------------------------------
// look for overlap between this class and input block
void relatedness_Block2::get_overlap(int block_left, int block_right, int block_parent_haplo_ID,
                                    int block_coal_time) {
  
  if (block_left > block_right) {
    Rcpp::stop("get_overlap error, left-right not in order");
  }
  
  //print("get overlap. child:", haplo_ID, "parent:", block_parent_haplo_ID, ", from", block_left, "to", block_right);
  
  // loop through blocks and check for overlap
  for (int i = 0; i < left.size(); ++i) {
    
    // there are 6 possible configurations. Work through them in order
    // 1.    ---  |     |
    // 2.    -----|--   |
    // 3.    -----|-----|-----
    // 4.         | --- |
    // 5.         |   --|-----
    // 6.         |     |  ---
    
    if (block_right < left[i]) {  // config 1
      //print("config1");
      
      // no overlap
      
    } else if ((block_left < left[i]) && (block_right >= left[i]) && (block_right <= right[i])) {  // config 2
      //print("config2");
      
      // overlap interval is from left[i] to block_right
      map_parent(left[i], block_right, block_parent_haplo_ID, descendant_IDs[i], time[i], block_coal_time);
      
    } else if ((block_left < left[i]) && (block_right > right[i])) {  // config 3
      //print("config3");
      
      // overlap interval is from left[i] to right[i]
      map_parent(left[i], right[i], block_parent_haplo_ID, descendant_IDs[i], time[i], block_coal_time);
      
    } else if ((block_left >= left[i]) && (block_right <= right[i])) {  // config 4
      //print("config4");
      
      // overlap interval is from block_left to block_right
      map_parent(block_left, block_right, block_parent_haplo_ID, descendant_IDs[i], time[i], block_coal_time);
      
    } else if ((block_left >= left[i]) && (block_left <= right[i]) && (block_right > right[i])) {  // config 5
      //print("config5");
      
      // overlap interval is from block_left to right[i]
      map_parent(block_left, right[i], block_parent_haplo_ID, descendant_IDs[i], time[i], block_coal_time);
      
    }  else if (block_left > right[i]) {  // config 6
      //print("config6");
      
      // no overlap
      
    } else {
      Rcpp::stop("error in relatedness_Block2::get_overlap(), end of switch block");
    }
  }
  
}

//------------------------------------------------
// look for parent. Add to map if doesn't exist, otherwise look for coalescence
void relatedness_Block2::map_parent(int block_left, int block_right, int block_parent_haplo_ID,
                                   vector<int> &block_descendants, int block_descendant_time, int block_coal_time) {
  
  if (block_left > block_right) {
    Rcpp::stop("map_parent error, left-right not in order");
  }
  
  //print("map parent. child:", haplo_ID, "parent:", block_parent_haplo_ID, "descendant:", block_descendant_haplo_ID);
  
  // if parent doesn't exist then add to map
  // if parent does exist then check for coalescence within parent object
  if (block_map_ptr->count({block_parent_haplo_ID, contig}) == 0) {
    //print("new parent", block_parent_haplo_ID, ", from", block_left, "to", block_right);
    
    // create new parent
    relatedness_Block2 this_parent = relatedness_Block2(block_parent_haplo_ID, contig, generation + 1,
                                                        block_left, block_right, block_descendant_time, block_descendants,
                                                        block_map_ptr, coaltracker_ptr, store_output_ptr);
    
    // add to map
    block_map_ptr->insert({{block_parent_haplo_ID, contig}, this_parent});
    
  } else {
    //print("existing parent");
    
    // check coalescence in existing parent
    block_map_ptr->at({block_parent_haplo_ID, contig}).check_coalescence(block_left, block_right, block_descendant_time, block_descendants, block_coal_time);
    
  }
  
  //print("end map_parent");
}

//------------------------------------------------
// check for and implement coalescence
void relatedness_Block2::check_coalescence(int block_left, int block_right, int block_descendant_time, vector<int> &block_descendants,
                                           int block_coal_time) {
  
  if (block_left > block_right) {
    Rcpp::stop("check_coalescence error, left-right not in order");
  }
  
  //print("check coalescence within parent:", haplo_ID);
  //print(block_descendant_time, block_coal_time);
  //print_vector(time);
  
  // check blocks in order
  if (!is_sorted(left.begin(), left.end())) {
    Rcpp::stop("blocks out of order");
  }
  
  // get number of blocks as currently stands
  int n_blocks = left.size();
  
  // loop through blocks and check for overlap
  bool leftover_block = true;
  for (int i = 0; i < n_blocks; ++i) {
    //print("  parent block", i, ": from", left[i], "to", right[i]);
    //print("  new block: from", block_left, "to", block_right);
    
    // there are 6 possible configurations. Work through them in order
    // 1.    ---  |     |
    // 2.    -----|--   |
    // 3.    -----|-----|-----
    // 4.         | --- |
    // 5.         |   --|-----
    // 6.         |     |  ---
    
    if (block_right < left[i]) {  // config 1
      //print("config1");
      
      // no overlap with first parental block. Add new block to this parent
      // without coalescence
      add_block(block_left, block_right, block_descendant_time, block_descendants);
      leftover_block = false;
      break;
      
    } else if ((block_left < left[i]) && (block_right >= left[i]) && (block_right <= right[i])) {  // config 2
      //print("config2");
      
      // block independent from block_left to left[i]-1
      // coalesce from left[i] to block_right
      coalescence(i, left[i], block_right, block_descendant_time, block_descendants, block_coal_time);
      add_block(block_left, left[i] - 1, block_descendant_time, block_descendants);
      if (block_right == right[i]) {
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      } else {
        add_block(block_right + 1, right[i], block_descendant_time, descendant_IDs[i]);
        right[i] = block_right;
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      }
      
      leftover_block = false;
      break;
      
    } else if ((block_left < left[i]) && (block_right > right[i])) {  // config 3
      //print("config3");
      
      // block independent from block_left to left[i]-1
      // coalesce from left[i] to right[i]
      // block continues from right[i]+1 to block_right
      coalescence(i, left[i], right[i], block_descendant_time, block_descendants, block_coal_time);
      add_block(block_left, left[i] - 1, block_descendant_time, block_descendants);
      push_back_multiple(descendant_IDs[i], block_descendants);
      time[i] = block_coal_time;
      block_left = right[i] + 1;
      
    } else if ((block_left >= left[i]) && (block_right <= right[i])) {  // config 4
      //print("config4");
      
      // coalesce from block_left to block_right
      coalescence(i, block_left, block_right, block_descendant_time, block_descendants, block_coal_time);
      if ((block_left == left[i]) && (block_right == right[i])) {
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      } else if (block_left == left[i]) {
        add_block(block_right + 1, right[i], block_descendant_time, block_descendants);
        right[i] = block_right;
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      } else if (block_right == right[i]) {
        add_block(left[i], block_left - 1, block_descendant_time, block_descendants);
        left[i] = block_left;
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      } else {
        add_block(left[i], block_left - 1, block_descendant_time, block_descendants);
        add_block(block_right + 1, right[i], block_descendant_time, block_descendants);
        left[i] = block_left;
        right[i] = block_right;
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      }
      
      leftover_block = false;
      break;
      
    } else if ((block_left >= left[i]) && (block_left <= right[i]) && (block_right > right[i])) {  // config 5
      //print("config5");
      
      // coalesce from block_left to right[i]
      // block continues from right[i]+1 to block_right
      coalescence(i, block_left, right[i], block_descendant_time, block_descendants, block_coal_time);
      if (block_left == left[i]) {
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      } else {
        add_block(left[i], block_left - 1, block_descendant_time, block_descendants);
        left[i] = block_left;
        push_back_multiple(descendant_IDs[i], block_descendants);
        time[i] = block_coal_time;
      }
      
      block_left = right[i] + 1;
      continue;
      
      
    }  else if (block_left > right[i]) {  // config 6
      //print("config6");
      
      // no overlap. Continue to next
      
    } else {
      Rcpp::stop("error in relatedness_Block2::check_coalescence(), end of switch block");
    }
    
  }  // end i loop
  
  // create leftover block
  if (leftover_block) {
    add_block(block_left, block_right, block_descendant_time, block_descendants);
  }
  
  // get blocks in inreasing order of genomic position
  get_blocks_in_order();
  
}

//------------------------------------------------
// implement coalescence
void relatedness_Block2::coalescence(int parent_index, int block_left, int block_right, int block_descendant_time,
                                     std::vector<int> &block_descendants, int block_coal_time) {
  
  
  //print("coalescence");
  //print("Parent chunk: from", left[parent_index], "to", right[parent_index]);
  //print("Input block : from", block_left, "to", block_right);
  //print(block_descendant_time, block_coal_time);
  
  //print("coaltracker before:");
  //coaltracker_ptr->at(contig).print_status();
  
  // update coalescent tracker
  coaltracker_ptr->at(contig).coalescence(block_left, block_right);
  
  //print("coaltracker after:");
  //coaltracker_ptr->at(contig).print_status();
  
  // add to record of coalescent events
  int t1 = (time[parent_index] < block_coal_time) ? time[parent_index] : block_coal_time;
  int t2 = (time[parent_index] < block_coal_time) ? block_coal_time : time[parent_index];
  
  Rcpp::List tmp1;
  tmp1["contig"] = contig + 1;
  tmp1["left"] = block_left;
  tmp1["right"] = block_right;
  tmp1["time_old"] = t1;
  tmp1["time_new"] = t2;
  tmp1["descendants"] = descendant_IDs[parent_index];
  store_output_ptr->push_back(tmp1);
  
  t1 = (block_descendant_time < block_coal_time) ? block_descendant_time : block_coal_time;
  t2 = (block_descendant_time < block_coal_time) ? block_coal_time : block_descendant_time;
  
  Rcpp::List tmp2;
  tmp2["contig"] = contig + 1;
  tmp2["left"] = block_left;
  tmp2["right"] = block_right;
  tmp2["time_old"] = t1;
  tmp2["time_new"] = t2;
  tmp2["descendants"] = block_descendants;
  store_output_ptr->push_back(tmp2);
}

//------------------------------------------------
// reorder blocks (if needed) to ensure they are ordered from left to right by
// start position
void relatedness_Block2::get_blocks_in_order() {
  
  // return if single element
  if (left.size() == 1) {
    return;
  }
  
  // return if already sorted
  if (is_sorted(left.begin(), left.end())) {
    return;
  }
  
  // get vectors in order
  vector<int> block_order = get_order(left);
  apply_order(left, block_order);
  apply_order(right, block_order);
  apply_order(time, block_order);
  apply_order(descendant_IDs, block_order);
  
}

//------------------------------------------------
// print all blocks
void relatedness_Block2::print_status() {
  print("haplo_ID:", haplo_ID, ", contig:", contig);
  for (int i = 0; i < left.size(); ++i) {
    print(left[i], right[i]);
  }
}
