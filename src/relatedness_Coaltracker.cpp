
#include "relatedness_Coaltracker.h"
#include "misc_v14.h"

using namespace std;

//------------------------------------------------
// initialise
void relatedness_Coaltracker::init(int contig, int left, int right, int lineages_remaining) {
  this->contig = contig;
  add_block(left, right, lineages_remaining);
}

//------------------------------------------------
// add new block to end
void relatedness_Coaltracker::add_block(int left, int right, int lineages_remaining) {
  this->left.push_back(left);
  this->right.push_back(right);
  this->lineages_remaining.push_back(lineages_remaining);
}

//------------------------------------------------
// find which sections still have more than one lineage remaining. Return in
// res_left and res_right
void relatedness_Coaltracker::get_overlap(int block_left, int block_right, vector<int> &res_left, vector<int> &res_right) {
  
  bool loop_open = false;
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
      
      // no overlap. Should not be possible to reach this stage as should have
      // broken loop in previous block
      print("args:", block_left, block_right);
      print_status();
      Rcpp::stop("coaltracker::get_overlap() should never reach config1");
      
    } else if ((block_left < left[i]) && (block_right >= left[i]) && (block_right <= right[i])) {  // config 2
      //print("config2");
      
      // block independent from block_left to left[i]-1
      // overlap from left[i] to block_right
      if (loop_open) {
        if (lineages_remaining[i] == 1) {
          res_right.push_back(left[i] - 1);
        } else {
          res_right.push_back(block_right);
        }
      } else {
        if (lineages_remaining[i] == 1) {
          // do nothing
        } else {
          res_left.push_back(left[i]);
          res_right.push_back(block_right);
        }
      }
      break;
      
    } else if ((block_left < left[i]) && (block_right > right[i])) {  // config 3
      //print("config3");
      
      // block independent from block_left to left[i]-1
      // overlap from left[i] to right[i]
      // block continues from right[i]+1 to block_right
      if (loop_open) {
        if (lineages_remaining[i] == 1) {
          res_right.push_back(left[i] - 1);
          loop_open = false;
        } else {
          // do nothing
        }
      } else {
        if (lineages_remaining[i] == 1) {
          // do nothing
        } else {
          res_left.push_back(left[i]);
          loop_open = true;
        }
      }
      
    } else if ((block_left >= left[i]) && (block_right <= right[i])) {  // config 4
      //print("config4");
      
      // overlap from block_left to block_right
      if (lineages_remaining[i] == 1) {
        // do nothing
      } else {
        res_left.push_back(block_left);
        res_right.push_back(block_right);
      }
      break;
      
    } else if ((block_left >= left[i]) && (block_left <= right[i]) && (block_right > right[i])) {  // config 5
      //print("config5");
      
      // overlap from block_left to right[i]
      // block continues from right[i]+1 to block_right
      if (lineages_remaining[i] == 1) {
        // do nothing
      } else {
        res_left.push_back(block_left);
        loop_open = true;
      }
      
    }  else if (block_left > right[i]) {  // config 6
      //print("config6");
      
      // no overlap
      
    } else {
      Rcpp::stop("error in relatedness_Coaltracker::get_overlap(), end of switch block");
    }
    
  }  // end loop through blocks
  
}

//------------------------------------------------
// implement coalescence
void relatedness_Coaltracker::coalescence(int block_left, int block_right) {
  
  //print("coaltracker: coalescence");
  //print("block_left:", block_left, "block_right:", block_right);
  
  int n_blocks = left.size();
  for (int i = 0; i < n_blocks; ++i) {
    //print("i =", i, "of", n_blocks - 1);
    
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
      
      // block independent from block_left to left[i]-1
      // coalesce from left[i] to block_right
      if (block_right == right[i]) {
        lineages_remaining[i]--;
      } else {
        add_block(block_right + 1, right[i], lineages_remaining[i]);
        right[i] = block_right;
        lineages_remaining[i]--;
      }
      
    } else if ((block_left < left[i]) && (block_right > right[i])) {  // config 3
      //print("config3");
      
      // block independent from block_left to left[i]-1
      // coalesce from left[i] to right[i]
      // block continues from right[i]+1 to block_right
      lineages_remaining[i]--;
      
    } else if ((block_left >= left[i]) && (block_right <= right[i])) {  // config 4
      //print("config4");
      
      // coalesce from block_left to block_right
      if ((block_left == left[i]) && (block_right == right[i])) {
        lineages_remaining[i]--;
      } else if (block_left == left[i]) {
        add_block(block_right + 1, right[i], lineages_remaining[i]);
        right[i] = block_right;
        lineages_remaining[i]--;
      } else if (block_right == right[i]) {
        add_block(left[i], block_left - 1, lineages_remaining[i]);
        left[i] = block_left;
        lineages_remaining[i]--;
      } else {
        add_block(left[i], block_left - 1, lineages_remaining[i]);
        add_block(block_right + 1, right[i], lineages_remaining[i]);
        left[i] = block_left;
        right[i] = block_right;
        lineages_remaining[i]--;
      }
      
    } else if ((block_left >= left[i]) && (block_left <= right[i]) && (block_right > right[i])) {  // config 5
      //print("config5");
      
      // coalesce from block_left to right[i]
      // block continues from right[i]+1 to block_right
      if (block_left == left[i]) {
        lineages_remaining[i]--;
      } else {
        add_block(left[i], block_left - 1, lineages_remaining[i]);
        left[i] = block_left;
        lineages_remaining[i]--;
      }
      
    }  else if (block_left > right[i]) {  // config 6
      //print("config6");
      
      // no overlap. Continue to next
      
    } else {
      Rcpp::stop("error in relatedness_Coaltracker::coalescence(), end of switch block");
    }
    
  }  // end loop through blocks
  
  // tidy up blocks
  tidy_blocks();
  
}

//------------------------------------------------
// remove redundancy by joining blocks together as needed
void relatedness_Coaltracker::tidy_blocks() {
  
  // nothing to tidy if single block
  if (left.size() == 1) {
    return;
  }
  
  // get blocks in order
  if (!is_sorted(left.begin(), left.end())) {
    vector<int> block_order = get_order(left);
    apply_order(left, block_order);
    apply_order(right, block_order);
    apply_order(lineages_remaining, block_order);
  }
  
  // loop through and replace blocks
  int i = 1;
  while (i < left.size()) {
    if (lineages_remaining[i-1] == lineages_remaining[i]) {
      right[i-1] = right[i];
      drop_block(i);
    } else {
      i++;
    }
  }
  
}

//------------------------------------------------
// remove ith block
void relatedness_Coaltracker::drop_block(int i) {
  left.erase(left.begin() + i);
  right.erase(right.begin() + i);
  lineages_remaining.erase(lineages_remaining.begin() + i);
}

//------------------------------------------------
// print status
void relatedness_Coaltracker::print_status() {
  print("coaltracker, contig:", contig);
  for (int i = 0; i < left.size(); ++i) {
    print(left[i], right[i], lineages_remaining[i]);
  }
}

