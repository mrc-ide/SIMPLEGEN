
#include "main.h"
#include "Parameters.h"
#include "Dispatcher.h"
#include "Haplo_node.h"
#include "probability_v10.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
Rcpp::List indiv_sim_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define parameters object and load values
  Parameters params;
  params.load_epi_params(args);
  params.load_deme_params(args);
  params.load_migration_params(args);
  params.load_demog_params(args);
  params.load_sampling_params(args);
  params.load_run_params(args);
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.init();
  dispatcher.run_simulation(args_functions, args_progress);
  
  // end timer
  chrono_timer(t1);
  
  return Rcpp::List::create(Rcpp::Named("daily_values") = dispatcher.daily_values,
                            Rcpp::Named("age_distributions") = dispatcher.age_distributions,
                            Rcpp::Named("sample_details") = dispatcher.sample_details);
}
#endif

//------------------------------------------------
// read in transmission record, and prune to only those events that run up to
// the final sample
#ifdef RCPP_ACTIVE
void prune_transmission_record_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract input args
  string transmission_record_location = rcpp_to_string(args["transmission_record_location"]);
  string pruned_record_location = rcpp_to_string(args["pruned_record_location"]);
  bool silent = rcpp_to_bool(args["silent"]);
  vector<int> inoc_IDs = rcpp_to_vector_int(args["inoc_IDs"]);
  
  // open filestream to transmission record and check that opened
  if (!silent) {
    print("Opening filestream to transmission record");
  }
  ifstream infile;
  infile.open(transmission_record_location);
  if (!infile.is_open()) {
    Rcpp::stop("unable to open filestream at specified location. Check the path exists, and that you have read access");
  }
  
  // open filestream to pruned record and check that opened
  if (!silent) {
    print("Opening filestream to pruned transmission record");
  }
  ofstream outfile;
  outfile.open(pruned_record_location);
  if (!outfile.is_open()) {
    Rcpp::stop("unable to open filestream at specified location. Check the path exists, and that you have write access");
  }
  
  
  // start by getting positions of all line breaks, and the first ID in every
  // line
  vector<int> line_starts(1);
  vector<pair<int, int>> first_IDs;
  stringstream ss;
  int ID;
  bool new_line = true;
  int i = 0;
  int t = 0;
  char c;
  while (infile.get(c)) {
    i++;
    if (c == '\n') {
      line_starts.push_back(i);
      new_line = true;
      t++;
    }
    if (new_line) {
      if (c == ' ' || c == ';') {
        ss >> ID;
        first_IDs.push_back(make_pair(t, ID));
        ss.clear();
        ss.str(std::string());
        new_line = false;
      } else {
        ss << c;
      }
    }
  }
  
  // count the number of generations. Note that the line_starts vector contains
  // (max_time+1) elements, as the final element represents the newlinw at the
  // end of the file
  int max_time = int(line_starts.size()) - 1;
  
  // create a vector of vectors for storing IDs we are interested in at each
  // time point. Initialise with inoc_IDs
  vector<vector<int>> pop(max_time);
  add_to_pop(inoc_IDs, first_IDs, pop);
  
  // create a map for storing final results. First value in the vector will
  // store the time
  map<int, vector<int>> pruned;
  
  // loop through the file again, this time reading lines in reverse order
  infile.clear();
  infile.seekg(0, ios::beg);
  ss.clear();
  ss.str(std::string());
  for (int t = (max_time-1); t >= 0; --t) {
    
    // some housekeeping
    int this_line_start = line_starts[t];
    int this_line_end = line_starts[t+1] - 1;
    int popsize = int(pop[t].size());
    inoc_IDs.clear();
    
    // go to start of line t
    infile.seekg(this_line_start, ios::beg);
    
    // loop through characters until reach end of line
    bool save_IDs_on = false;
    bool new_block = true;
    int ID_key;
    for (int i = this_line_start; i < this_line_end; ++i) {
      infile.get(c);
      
      // process character, either storing or discarding
      if (new_block) {
        if (c == ' ' || c == ';') {
          ss >> ID_key;
          ss.clear();
          ss.str(std::string());
          if (c == ' ') {
            new_block = false;
          }
          
          // if ID_key is in pop[t] then store time, and activate ID saving for
          // remainder of this block
          save_IDs_on = false;
          for (int j = 0; j < popsize; ++j) {
            if (pop[t][j] == ID_key) {
              pruned[ID_key].push_back(t);
              save_IDs_on = true;
              break;
            }
          }
          
        } else {
          ss << c;
        }
      } else {
        if (save_IDs_on) {
          if (c == ' ' || c == ';') {
            ss >> ID;
            pruned[ID_key].push_back(ID);
            ss.clear();
            ss.str(std::string());
            inoc_IDs.push_back(ID);
            if (c == ';') {
              new_block = true;
            }
          } else {
            ss << c;
          }
        } else if (c == ';') {
          ss.clear();
          ss.str(std::string());
          new_block = true;
        }
      }
      
      
    }  // end i loop
    
    // add new inoc_IDs to pop
    add_to_pop(inoc_IDs, first_IDs, pop);
    
  }  // end t loop
  
  // write pruned record to file
  t = 0;
  for (const auto & x : pruned) {
    while (x.second[0] > t) {
      outfile << "\n";
      t++;
    }
    outfile << x.first;
    for (int i = 1; i < int(x.second.size()); ++i) {
      outfile << " " << x.second[i];
    }
    outfile << ";";
  }
  
  // close filestreams
  if (!silent) {
    print("Closing filestreams");
  }
  infile.close();
  outfile.close();
  
  // end timer
  chrono_timer(t1);
  
}
#endif

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
// read in a pruned transmission record, and simulate relatedness intervals
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
  
  // create a map, where each element represents a node in the pruned
  // transmission tree
  map<int, Haplo_node> haplo_tree;
  
  // populate map by parsing puned transmission record
  vector<int> inoc_ID_vec;
  string line, block, element;
  int haplo_ID = 0;
  int t = 0;
  while (getline(infile, line)) {
    istringstream iss(line);
    while (getline(iss, block, ';')) {
      
      // parse block into inoc IDs
      istringstream iss2(block);
      inoc_ID_vec.clear();
      while (getline(iss2, element, ' ')) {
        int inoc_ID = stoi(element);
        inoc_ID_vec.push_back(inoc_ID);
      }
      
      // initialise node using the key ID
      haplo_tree[inoc_ID_vec[0]] = Haplo_node(t, contig_lengths, sampler_oocyst, sampler_hepatocyte, haplo_tree);
      
      // create haplotypes within this node either de novo or from ancestral
      // inoculations
      if (inoc_ID_vec.size() == 1) {
        haplo_tree[inoc_ID_vec[0]].draw_haplotypes_denovo(haplo_ID, alpha);
      } else {
        haplo_tree[inoc_ID_vec[0]].draw_haplotypes_recombine(haplo_ID, inoc_ID_vec, r, alpha);
      }
      
    }
    t++;
  }
  
  // DEBUG - DELETE THE FOLLOWING LINES ONCE COMPLETE
  return Rcpp::List::create(Rcpp::Named("details") = -9);
  
  
  /*
  // get into more convenient output format
  Rcpp::List ret;
  for (int i = 0; i < int(pruned_array.size()); ++i) {
    int ID_key = pruned_array[i][0];
    
    // get descriptive details of this node
    Rcpp::List ret_details;
    ret_details["inoc_ID"] = ID_key;
    ret_details["time"] = tree[ID_key].t;
    ret_details["haplo_IDs"] = tree[ID_key].haplo_ID_vec;
    ret_details["haplo_densities"] = tree[ID_key].haplo_density;
    
    // get haplotype intervals
    Rcpp::List ret_haplotypes;
    for (int j = 0; j < tree[ID_key].n_haplotypes; ++j) {
      ret_haplotypes.push_back(tree[ID_key].intervals[j]);
    }
    
    // push to return list
    ret.push_back(Rcpp::List::create(Rcpp::Named("details") = ret_details,
                                     Rcpp::Named("haplotypes") = ret_haplotypes));
  }
  
  // close filestreams
  if (!silent) {
    print("Closing filestream");
  }
  infile.close();
  
  // end timer
  chrono_timer(t1);
  
  return ret;
  */
}

