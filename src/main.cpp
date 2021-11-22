
#include "main.h"
#include "Parameters.h"
#include "Dispatcher.h"
#include "Tree_node.h"
#include "probability_v11.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
Rcpp::List indiv_sim_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define parameters object and load values from arguments
  Parameters params;
  params.load_params(args);
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.init(params);
  dispatcher.run_simulation(args_functions, args_progress);
  
  // end timer
  chrono_timer(t1);
  
  return Rcpp::List::create(Rcpp::Named("daily_output") = dispatcher.daily_output,
                            Rcpp::Named("sweep_output") = dispatcher.sweep_output,
                            Rcpp::Named("survey_output") = dispatcher.surveys_indlevel_output,
                            Rcpp::Named("survey_output_infection_IDs") = dispatcher.surveys_indlevel_output_infection_IDs);
  
}
#endif

//------------------------------------------------
// read in transmission record and prune to only those events that run up to the
// final sample
void prune_transmission_record_cpp(Rcpp::List args) {
  
  // extract input args
  string transmission_record_location = rcpp_to_string(args["transmission_record_location"]);
  string pruned_record_location = rcpp_to_string(args["pruned_record_location"]);
  vector<int> focal_IDs_vec = rcpp_to_vector_int(args["infection_IDs"]);
  set<int> focal_IDs(focal_IDs_vec.begin(), focal_IDs_vec.end());
  bool silent = rcpp_to_bool(args["silent"]);
  
  // open filestream to transmission record and check that opened
  if (!silent) {
    print("Opening filestream to transmission record");
  }
  ifstream infile;
  infile.open(transmission_record_location);
  if (!infile.is_open()) {
    Rcpp::stop("unable to open filestream at specified location. Check the path exists and that you have read access");
  }
  
  // open filestream to pruned record and check that opened
  if (!silent) {
    print("Opening filestream to pruned transmission record");
  }
  ofstream outfile;
  outfile.open(pruned_record_location);
  if (!outfile.is_open()) {
    Rcpp::stop("unable to open filestream at specified location. Check the path exists and that you have write access");
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
        event.push_back(stoi(field));
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
  
  // write pruned record to file
  outfile << "time,event,human_ID,mosquito_ID,child_infection_ID,parent_infection_ID,deme\n";
  for (int i = time_pruned.size() - 1; i >= 0; --i) {
    outfile << time_pruned[i] << "," << event_pruned[i] << "," << human_ID_pruned[i] << "," << mosquito_ID_pruned[i] << "," << child_infection_ID_pruned[i] << "," << parent_infection_ID_pruned[i][0];
    for (int j = 1; j < parent_infection_ID_pruned[i].size(); ++j) {
      outfile << ";" << parent_infection_ID_pruned[i][j];
    }
    outfile << "," << deme_pruned[i] << "\n";
  }
  
  // close filestreams
  if (!silent) {
    print("Closing filestreams");
  }
  infile.close();
  outfile.close();
  
}

//------------------------------------------------
// read in a vector of inoc_IDs. These are put in ascending order, and searched
// for using the information in first_IDs to find the time at which each inoc_ID
// is first created. Finally, inoc_IDs are added to the population array at
// these time points. IDs are only added if they are not already present in the
// population array at that time point.
void add_to_pop(vector<int> &inoc_IDs, const vector<pair<int, int>> &first_IDs,
                vector<vector<int>> &pop, bool print_pop) {
  
  // sort inoc_IDs
  sort(inoc_IDs.begin(), inoc_IDs.end());
  
  // push back IDs to population at correct time
  int j = 0;
  bool break_i = false;
  for (int i = 0; i < (int(first_IDs.size()) - 1); ++i) {
    if (break_i) {
      break;
    }
    while (inoc_IDs[j] < first_IDs[i+1].second) {
      int t = first_IDs[i].first;
      if (find(pop[t].begin(), pop[t].end(), inoc_IDs[j]) == pop[t].end()) {
        pop[t].push_back(inoc_IDs[j]);
      }
      j++;
      if (j > (int(inoc_IDs.size()) - 1)) {
        break_i = true;
        break;
      }
    }
  }
  
  // any leftover IDs must be beyond the end of first_IDs, and so are added to
  // final timepoint
  for (int i = j; i < int(inoc_IDs.size()); ++i) {
    int t = first_IDs[first_IDs.size()-1].first;
    if (find(pop[t].begin(), pop[t].end(), inoc_IDs[i]) == pop[t].end()) {
      pop[t].push_back(inoc_IDs[i]);
    }
  }
  
  // optionally print pop array
  if (print_pop) {
    for (int t = 0; t < int(pop.size()); ++t) {
      Rcpp::Rcout << "t = " << t << ": ";
      for (int i = 0; i < int(pop[t].size()); ++i) {
        Rcpp::Rcout << pop[t][i] << " ";
      }
      Rcpp::Rcout << "\n";
    }
  }
  
}

//------------------------------------------------
// read in a pruned transmission record and simulate relatedness intervals
// forwards-in-time through this tree
Rcpp::List sim_relatedness_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract input args
  string pruned_record_location = rcpp_to_string(args["pruned_record_location"]);
  double r = rcpp_to_double(args["r"]);
  double alpha = rcpp_to_double(args["alpha"]);
  vector<double> oocyst_distribution = rcpp_to_vector_double(args["oocyst_distribution"]);
  vector<double> hepatocyte_distribution = rcpp_to_vector_double(args["hepatocyte_distribution"]);
  vector<int> contig_lengths = rcpp_to_vector_int(args["contig_lengths"]);
  bool silent = rcpp_to_bool(args["silent"]);
  
  // open filestream to pruned record and check that opened
  if (!silent) {
    print("Opening filestream to pruned transmission record");
  }
  ifstream infile;
  infile.open(pruned_record_location);
  if (!infile.is_open()) {
    Rcpp::stop("unable to open filestream at specified location. Check the path exists, and that you have read access");
  }
  
  // create samplers for drawing number of oocysts and hepatocytes
  Sampler sampler_oocyst(oocyst_distribution, 1000);
  Sampler sampler_hepatocyte(hepatocyte_distribution, 1000);
  
  // read pruned record into an array. For each row, first value gives innoc_ID,
  // second value gives time of creation, and any subsequent values give
  // parental inoc_IDs that are ancestral
  vector<vector<int>> pruned_array;
  vector<int> this_line;
  string line, block, element;
  int t = 0;
  while (getline(infile, line)) {
    istringstream iss(line);
    while (getline(iss, block, ';')) {
      bool new_block = true;
      istringstream iss2(block);
      while (getline(iss2, element, ' ')) {
        
        if (new_block) {
          this_line.clear();
          int ID_key = stoi(element);
          this_line.push_back(ID_key);
          this_line.push_back(t);
          new_block = false;
        } else {
          int ID_val = stoi(element);
          this_line.push_back(ID_val);
        }
        
      }
      pruned_array.push_back(this_line);
    }
    t++;
  }
  
  // create a map, where each element represents an inoculation in the pruned
  // transmission tree
  map<int, Tree_node> inoc_tree;
  
  // populate map
  int lineage_ID = 1;
  for (int i = 0; i < int(pruned_array.size()); ++i) {
    
    // create node for this inoc_ID
    int inoc_ID = pruned_array[i][0];
    int t = pruned_array[i][1];
    inoc_tree[inoc_ID] = Tree_node(t, contig_lengths, sampler_oocyst, sampler_hepatocyte, inoc_tree);
    
    // draw lineages de novo or by recombination from ancestral inoculations
    if (pruned_array[i].size() == 2) {
      inoc_tree[inoc_ID].draw_lineages_denovo(lineage_ID, alpha);
    } else {
      inoc_tree[inoc_ID].draw_lineages_recombine(lineage_ID, pruned_array[i], r, alpha);
    }
  }
  
  // get into more convenient output format
  Rcpp::List ret;
  for (int i = 0; i < int(pruned_array.size()); ++i) {
    int inoc_ID = pruned_array[i][0];
    
    // get descriptive details of this node
    Rcpp::List ret_details;
    ret_details["inoc_ID"] = inoc_ID;
    ret_details["time"] = inoc_tree[inoc_ID].t;
    ret_details["lineage_IDs"] = inoc_tree[inoc_ID].lineage_ID_vec;
    ret_details["lineage_densities"] = inoc_tree[inoc_ID].lineage_density;
    
    // get lineage intervals
    Rcpp::List ret_lineages;
    for (int j = 0; j < inoc_tree[inoc_ID].n_lineages; ++j) {
      ret_lineages.push_back(inoc_tree[inoc_ID].intervals[j]);
    }
    
    // push to return list
    ret.push_back(Rcpp::List::create(Rcpp::Named("details") = ret_details,
                                     Rcpp::Named("lineages") = ret_lineages));
  }
  
  // close filestreams
  if (!silent) {
    print("Closing filestream");
  }
  infile.close();
  
  // end timer
  chrono_timer(t1);
  
  
  return ret;
}
