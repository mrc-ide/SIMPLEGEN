
#include <set>
#include "prune_main.h"

using namespace std;

//------------------------------------------------
cpp11::list prune_transmission_record(cpp11::list args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract input args
  string transmission_record_location = cpp_to_string(args["transmission_record_location"]);
  vector<int> focal_IDs_vec = cpp_to_vector_int(args["inoc_IDs"]);
  set<int> focal_IDs(focal_IDs_vec.begin(), focal_IDs_vec.end());
  bool silent = cpp_to_bool(args["silent"]);
  
  // open filestream to transmission record and check that opened
  if (!silent) {
    print("Opening filestream to transmission record");
  }
  ifstream infile;
  infile.open(transmission_record_location);
  if (!infile.is_open()) {
    cpp11::stop("unable to open filestream at specified location. Check the path exists and that you have read access");
  }
  
  if (!silent) {
    print("Pruning transmission record");
  }
  
  // initialise vectors for storing transmission record values
  vector<int> time;
  vector<int> event;
  vector<int> human_ID;
  vector<int> mosquito_ID;
  vector<int> child_infection_ID;
  vector<vector<int>> parent_infection_ID;
  vector<int> deme;
  
  // read in transmission record to memory
  string line;
  int i = 0;
  while (getline(infile, line)) {
    
    // skip first line
    if (i == 0) {
      i++;
      continue;
    }
    
    
    
    // read line elements
    istringstream s(line);
    string field;
    int j = 0;
    while (getline(s, field,',')) {
      if (j == 0) {
        time.push_back(stoi(field));
      } else if (j == 1) {
        if (field == "H") {
          event.push_back(0);
        } else if (field == "M") {
          event.push_back(1);
        } else {
          cpp11::stop("unrecognised letter in transmission record (not H or M)");
        }
      } else if (j == 2) {
        human_ID.push_back(stoi(field));
      } else if (j == 3) {
        mosquito_ID.push_back(stoi(field));
      } else if (j == 4) {
        child_infection_ID.push_back(stoi(field));
      } else if (j == 5) {
        // read in elements delimited by ;
        stringstream ss(field);
        string ss_segment;
        vector<int> tmp;
        while (std::getline(ss, ss_segment, ';')) {
          tmp.push_back(stoi(ss_segment));
        }
        parent_infection_ID.push_back(tmp);
      } else if (j == 6) {
        deme.push_back(stoi(field));
      }
      j++;
    }
    i++;
  }
  
  // initialise vectors for storing pruned values
  vector<int> time_pruned;
  vector<int> event_pruned;
  vector<int> human_ID_pruned;
  vector<int> mosquito_ID_pruned;
  vector<int> child_infection_ID_pruned;
  vector<vector<int>> parent_infection_ID_pruned;
  vector<int> deme_pruned;
  
  // search backwards through transmission record
  for (int i = time.size() - 1; i >= 0; --i) {
    if (focal_IDs.count(child_infection_ID[i])) {
      
      // add to pruned record
      time_pruned.push_back(time[i]);
      event_pruned.push_back(event[i]);
      human_ID_pruned.push_back(human_ID[i]);
      mosquito_ID_pruned.push_back(mosquito_ID[i]);
      child_infection_ID_pruned.push_back(child_infection_ID[i]);
      parent_infection_ID_pruned.push_back(parent_infection_ID[i]);
      deme_pruned.push_back(deme[i]);
      
      // drop child ID from focal set and add parental IDs
      focal_IDs.erase(child_infection_ID[i]);
      for (int j = 0; j < parent_infection_ID[i].size(); ++j) {
        focal_IDs.insert(parent_infection_ID[i][j]);
      }
    }
  }
  
  // close filestreams
  if (!silent) {
    print("Closing filestreams");
  }
  infile.close();
  
  // end timer
  if (!silent) {
    chrono_timer(t1);
  }
  
  // write output list
  //cpp11::writable::integers_matrix<> parent_infection_ID_pruned_mat = copy_mat(parent_infection_ID_pruned);
  
  cpp11::writable::list parent_infection_ID_pruned_list;
  for (size_t i = 0; i < parent_infection_ID_pruned.size(); ++i) {
    parent_infection_ID_pruned_list.push_back({cpp11::literals::operator""_nm("z", 0) = parent_infection_ID_pruned[i]});
  }
  
  cpp11::writable::list tmp_output;
  tmp_output.push_back({cpp11::literals::operator""_nm("time", 0) = time_pruned});
  tmp_output.push_back({cpp11::literals::operator""_nm("event", 0) = event_pruned});
  tmp_output.push_back({cpp11::literals::operator""_nm("human_ID", 0) = human_ID_pruned});
  tmp_output.push_back({cpp11::literals::operator""_nm("mosquito_ID", 0) = mosquito_ID_pruned});
  tmp_output.push_back({cpp11::literals::operator""_nm("child_infection_ID", 0) = child_infection_ID_pruned});
  tmp_output.push_back({cpp11::literals::operator""_nm("parent_infection_ID", 0) = parent_infection_ID_pruned_list});
  tmp_output.push_back({cpp11::literals::operator""_nm("deme", 0) = deme_pruned});
  return tmp_output;
}
