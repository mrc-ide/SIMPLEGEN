
#include "main.h"
#include "Parameters.h"
#include "Dispatcher.h"
#include "probability_v14.h"
#include "genmodel_Host.h"
#include "genmodel_Mosquito.h"
#include "genmodel_Infection.h"
#include "relatedness_Haplo.h"

#include <chrono>
#include <map>

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
Rcpp::List prune_transmission_record_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract input args
  string transmission_record_location = rcpp_to_string(args["transmission_record_location"]);
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
  
  // close filestreams
  if (!silent) {
    print("Closing filestreams");
  }
  infile.close();
  
  // end timer
  if (!silent) {
    chrono_timer(t1);
  }
  
  // return as list
  return Rcpp::List::create(Rcpp::Named("time") = time_pruned,
                            Rcpp::Named("event") = event_pruned,
                            Rcpp::Named("human_ID") = human_ID_pruned,
                            Rcpp::Named("mosquito_ID") = mosquito_ID_pruned,
                            Rcpp::Named("child_infection_ID") = child_infection_ID_pruned,
                            Rcpp::Named("parent_infection_ID") = parent_infection_ID_pruned,
                            Rcpp::Named("deme") = deme_pruned);
}

//------------------------------------------------
// read in a pruned transmission record and simulate tree connecting haplotypes
// from genetic model
Rcpp::List sim_haplotype_tree_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract input args
  vector<int> sample_human_IDs = rcpp_to_vector_int(args["sample_human_IDs"]);
  vector<vector<int>> sample_infection_IDs = rcpp_to_matrix_int(args["sample_infection_IDs"]);
  
  vector<double> oocyst_distribution = rcpp_to_vector_double(args["oocyst_distribution"]);
  vector<double> hepatocyte_distribution = rcpp_to_vector_double(args["hepatocyte_distribution"]);
  double alpha = rcpp_to_double(args["alpha"]);
  bool defined_densities = rcpp_to_bool(args["defined_densities"]);
  bool defined_deme = rcpp_to_bool(args["defined_deme"]);
  bool silent = rcpp_to_bool(args["silent"]);
  
  Rcpp::List pruned_record = args["pruned_record"];
  vector<int> time = rcpp_to_vector_int(pruned_record["time"]);
  vector<int> event = rcpp_to_vector_int(pruned_record["event"]);
  vector<int> human_ID = rcpp_to_vector_int(pruned_record["human_ID"]);
  vector<int> mosquito_ID = rcpp_to_vector_int(pruned_record["mosquito_ID"]);
  vector<int> child_infection_ID = rcpp_to_vector_int(pruned_record["child_infection_ID"]);
  vector<vector<int>> parent_infection_ID = rcpp_to_matrix_int(pruned_record["parent_infection_ID"]);
  vector<vector<double>> parent_infection_density;
  if (defined_densities) {
    parent_infection_density = rcpp_to_matrix_double(pruned_record["parent_infection_density"]);
  }
  vector<vector<int>> deme;
  if (defined_deme) {
    deme = rcpp_to_matrix_int(pruned_record["deme"]);
  }
  
  // create samplers for drawing number of oocysts and hepatocytes
  Sampler sampler_oocyst(oocyst_distribution, 1000);
  Sampler sampler_hepatocyte(hepatocyte_distribution, 1000);
  
  // create dynamic populations of host and mosquito objects
  map<int, genmodel_Host> host_pop;
  map<int, genmodel_Mosquito> mosquito_pop;
  
  // initialise next haplo ID
  int next_haplo_ID = 1;
  
  // objects for storing results
  vector<int> output_record_row;
  vector<int> output_child_haplo_ID;
  vector<vector<int>> output_parent_haplo_ID;
  
  if (!silent) {
    print("Simulating haplotype tree");
  }
  
  // loop through pruned transmission record
  for (int i = 0; i < time.size(); ++i) {
    
    if (event[i] == 1) {
      // human receives parasites from mosquito
      
      // add host to population if doesn't exist
      if (host_pop.count(human_ID[i]) == 0) {
        host_pop[human_ID[i]] = genmodel_Host();
      }
      
      // pass in haplos from mosquito, or alternatively create de novo
      if (mosquito_ID[i] == -1) {
        host_pop[human_ID[i]].denovo_infection(i, child_infection_ID[i], next_haplo_ID, alpha, output_record_row,
                                               output_child_haplo_ID, output_parent_haplo_ID);
      } else {
        
        // check that mosquito exists in population
        if (mosquito_pop.count(mosquito_ID[i]) == 0) {
          Rcpp::stop("human host receives infectious bite from mosquito that has not yet been infected");
        }
        
        // pass in haplos from mosquito
        host_pop[human_ID[i]].infection(i, child_infection_ID[i], sampler_hepatocyte, next_haplo_ID, alpha,
                                        mosquito_pop[mosquito_ID[i]], parent_infection_ID[i],
                                        output_record_row, output_child_haplo_ID, output_parent_haplo_ID);
      }
      
      //host_pop[human_ID[i]].print_status();
      
    } else if (event[i] == 2) {
      // mosquito receives parasites from human
      
      // mosquito should not already exist in population as no superinfection
      // allowed in this model
      if (mosquito_pop.count(mosquito_ID[i]) != 0) {
        Rcpp::stop("same mosquito infected multiple times, which is not allowed in no-mosquito-superinfection model");
      }
      
      // corresponding host should exist in population
      if (host_pop.count(human_ID[i]) == 0) {
        Rcpp::stop("mosquito receives infectious bite from host who has not themselves been infected");
      }
      
      // add new mosquito to population
      mosquito_pop[mosquito_ID[i]] = genmodel_Mosquito();
      
      // get infection densities
      vector<double> this_infection_density;
      if (defined_densities) {
        this_infection_density = parent_infection_density[i];
      } else {
        this_infection_density = vector<double>(parent_infection_ID[i].size(), 1.0);
      }
      
      // infect mosquito from host
      mosquito_pop[mosquito_ID[i]].infection(i, child_infection_ID[i], sampler_oocyst,
                                             host_pop[human_ID[i]], parent_infection_ID[i], this_infection_density);
      
      //mosquito_pop[mosquito_ID[i]].print_status();
    }
    
  }  // end loop through transmission record
  
  // get haplo IDs in final sample
  vector<vector<vector<int>>> sample_haplo_IDs(sample_human_IDs.size());
  vector<vector<vector<double>>> sample_haplo_densities(sample_human_IDs.size());
  for (int i = 0; i < sample_human_IDs.size(); ++i) {
    
    // check that host exists in population
    if (host_pop.count(sample_human_IDs[i]) == 0) {
      Rcpp::stop("Host requested from sample details is not impled by pruned transmission record");
    }
    genmodel_Host this_host = host_pop[sample_human_IDs[i]];
    
    // initialise output objects
    sample_haplo_IDs[i] = vector<vector<int>>(sample_infection_IDs[i].size());
    sample_haplo_densities[i] = vector<vector<double>>(sample_infection_IDs[i].size());
    
    // loop through infection IDs
    for (int j = 0; j < sample_infection_IDs[i].size(); ++j) {
      int this_infection_ID = sample_infection_IDs[i][j];
      
      // check that this infection ID present in host
      if (this_host.infections_map.count(this_infection_ID) == 0) {
        Rcpp::stop("Infection ID requested from sample details not present in host");
      }
      
      // get haplos and densities
      genmodel_Infection this_infection = this_host.infections_map[this_infection_ID];
      sample_haplo_IDs[i][j] = this_infection.haplo_IDs;
      sample_haplo_densities[i][j] = this_infection.haplo_densities;
      
    }
  }
  
  
  // end timer
  if (!silent) {
    chrono_timer(t1);
  }
  
  return Rcpp::List::create(Rcpp::Named("record_row") = output_record_row,
                            Rcpp::Named("child_haplo_ID") = output_child_haplo_ID,
                            Rcpp::Named("parent_haplo_ID") = output_parent_haplo_ID,
                            Rcpp::Named("sample_haplo_IDs") = sample_haplo_IDs,
                            Rcpp::Named("sample_haplo_densities") = sample_haplo_densities);
}

//------------------------------------------------
// get relatedness between haplotypes by considering all possible paths through
// haplotype tree
Rcpp::List get_haplotype_relatedness_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract input args
  vector<int> time = rcpp_to_vector_int(args["time"]);
  vector<int> child_haplo_ID = rcpp_to_vector_int(args["child_haplo_ID"]);
  vector<vector<int>> parent_haplo_ID = rcpp_to_matrix_int(args["parent_haplo_ID"]);
  vector<int> target_haplo_IDs = rcpp_to_vector_int(args["target_haplo_IDs"]);
  int generations = rcpp_to_double(args["generations"]);
  bool silent = rcpp_to_bool(args["silent"]);
  
  // inialise map for storing a flexible population of haplotypes
  map<int, relatedness_Haplo> haplo_pop;
  
  // populate map with target haplo IDs
  for (int i = 0; i < target_haplo_IDs.size(); ++i) {
    haplo_pop[target_haplo_IDs[i]] = relatedness_Haplo(target_haplo_IDs[i]);
  }
  
  if (!silent) {
    print("Traversing haplotype tree");
  }
  
  // walk back up the tree
  for (int i = time.size() - 1; i >= 0; --i) {
    
    // skip if not a haplo we are looking for
    if (haplo_pop.count(child_haplo_ID[i]) == 0) {
      continue;
    }
    
    // skip if any descendent of this haplo exceeds generation limit
    if (max(haplo_pop[child_haplo_ID[i]].generations) >= generations) {
      continue;
    }
    
    // if single parent (i.e. human sampling of mosquito recombinant product)
    int this_child = child_haplo_ID[i];
    if (parent_haplo_ID[i].size() == 1) {
      int this_parent = parent_haplo_ID[i][0];
      
      // skip if no real parent (child haplo created de novo)
      if (this_parent == -1) {
        continue;
      }
      
      // create parental haplo from child if it does not exist in
      // population, otherwise merge child with existing
      if (haplo_pop.count(this_parent) == 0) {
        haplo_pop[this_parent] = haplo_pop[this_child];
      } else {
        haplo_pop[this_parent].merge(haplo_pop[this_child]);
      }
      
      // drop child from pop
      haplo_pop.erase(this_child);
      
    } else {
    // if two parents (i.e. recombinant sampling of parental gametocytes)
      
      // if parents are clonal
      if (parent_haplo_ID[i][0] == parent_haplo_ID[i][1]) {
        int this_parent = parent_haplo_ID[i][0];
        
        // increment generations but not recombinations
        haplo_pop[this_child].increment_generations();
        
        // create parental haplo from child if it does not exist in
        // population, otherwise merge child with existing
        if (haplo_pop.count(this_parent) == 0) {
          haplo_pop[this_parent] = haplo_pop[this_child];
        } else {
          haplo_pop[this_parent].merge(haplo_pop[this_child]);
        }
        
        // drop child from pop
        haplo_pop.erase(this_child);
        
      } else {
      // not clonal parents
      
        int parent1 = parent_haplo_ID[i][0];
        int parent2 = parent_haplo_ID[i][1];
        
        // increment generations and recombinations
        haplo_pop[this_child].increment_generations();
        haplo_pop[this_child].increment_recombinations();
        
        // for each parent, create parental haplo from child if it does not
        // exist in population, otherwise merge child with existing
        if (haplo_pop.count(parent1) == 0) {
          haplo_pop[parent1] = haplo_pop[this_child];
        } else {
          haplo_pop[parent1].merge(haplo_pop[this_child]);
        }
        if (haplo_pop.count(parent2) == 0) {
          haplo_pop[parent2] = haplo_pop[this_child];
        } else {
          haplo_pop[parent2].merge(haplo_pop[this_child]);
        }
        
        // drop child from pop
        haplo_pop.erase(this_child);
        
      }
      
    }
    
    
  }  // end loop up tree
  
  // print haplo pop
  //for (auto &x : haplo_pop) {
  //  print("Haplo:", x.first);
  //  x.second.print_status();
  //}
  
  //---------------------
  // Calculate relatedness between target haplotypes
  
  // initialise objects for storing pairwise relatedness
  int n_target_haplo_IDs = target_haplo_IDs.size();
  vector<double> pairwise_relatedness(0.5 * n_target_haplo_IDs * (n_target_haplo_IDs - 1));
  
  // loop through all pairwise combos of target haplos
  int i2 = 0;
  for (int i = 0; i < (n_target_haplo_IDs - 1); ++i) {
    for (int j = (i + 1); j < n_target_haplo_IDs; ++j) {
      
      // sum relatedness over all members of haplo_pop that have targets i and j
      // as descendants
      for (auto &x : haplo_pop) {
        bool descendant_i = false;
        bool descendant_j = false;
        for (int k = 0; k < x.second.descendant_IDs.size(); ++k) {
          if (x.second.descendant_IDs[k] == target_haplo_IDs[i]) {
            descendant_i = true;
          }
          if (x.second.descendant_IDs[k] == target_haplo_IDs[j]) {
            descendant_j = true;
          }
          if (descendant_i && descendant_j) {
            break;
          }
        }
        if (!(descendant_i && descendant_j)) {
          continue;
        }
        
        // sum relatedness to each descendant
        double relatedness_i = 0.0;
        double relatedness_j = 0.0;
        for (int k = 0; k < x.second.descendant_IDs.size(); ++k) {
          if (x.second.descendant_IDs[k] == target_haplo_IDs[i]) {
            relatedness_i += pow(0.5, x.second.recombinations[k]);
          }
          if (x.second.descendant_IDs[k] == target_haplo_IDs[j]) {
            relatedness_j += pow(0.5, x.second.recombinations[k]);
          }
        }
        
        // multiply to get relatedness i to j and add to running sum
        pairwise_relatedness[i2] += relatedness_i * relatedness_j;
        
      }  // end loop through haplo_pop
      
      i2++;
    }
  }
  
  // end timer
  if (!silent) {
    chrono_timer(t1);
  }
  
  return Rcpp::List::create(Rcpp::Named("pairwise_relatedness") = pairwise_relatedness);
}

//------------------------------------------------
// simulate recombination blocks from haplotype tree
Rcpp::List sim_block_tree_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract input args
  double r = rcpp_to_double(args["r"]);
  vector<int> contig_lengths = rcpp_to_vector_int(args["contig_lengths"]);
  int n_contigs = contig_lengths.size();
  vector<int> time = rcpp_to_vector_int(args["time"]);
  vector<int> child_haplo_ID = rcpp_to_vector_int(args["child_haplo_ID"]);
  vector<vector<int>> parent_haplo_ID = rcpp_to_matrix_int(args["parent_haplo_ID"]);
  bool silent = rcpp_to_bool(args["silent"]);
  
  // define objects for storing output
  vector<int> output_time;
  vector<int> output_child_haplo_ID;
  vector<int> output_contig;
  vector<int> output_start;
  vector<int> output_end;
  vector<int> output_parent_haplo_ID;
  
  // step through haplotype tree
  for (int i = 0; i < time.size(); ++i) {
    
    // establish if clonal
    bool clonal = true;
    if (parent_haplo_ID[i].size() == 2) {
      if (parent_haplo_ID[i][0] != parent_haplo_ID[i][1]) {
        clonal = false;
      }
    }
    
    // if clonal then same ancestor over all contigs
    if (clonal) {
      vector<int> this_time(n_contigs, time[i]);
      vector<int> this_child(n_contigs, child_haplo_ID[i]);
      vector<int> this_contig = seq_int(1, n_contigs);
      vector<int> this_start(n_contigs, 1);
      vector<int> this_parent(n_contigs, parent_haplo_ID[i][0]);
      
      push_back_multiple(output_time, this_time);
      push_back_multiple(output_child_haplo_ID, this_child);
      push_back_multiple(output_contig, this_contig);
      push_back_multiple(output_start, this_start);
      push_back_multiple(output_end, contig_lengths);
      push_back_multiple(output_parent_haplo_ID, this_parent);
    }
    
    // if not clonal then simulate recombination
    if (!clonal) {
      
      // loop through contigs
      vector<int> this_start;
      vector<int> this_end;
      vector<int> this_contig;
      vector<int> this_parent;
      for (int j = 0; j < n_contigs; ++j) {
        
        // draw starting parent
        int parent_index = rbernoulli1(0.5);
        
        // draw blocks until reach end of contig
        int block_start = 1;
        while (block_start <= contig_lengths[j]) {
          
          // draw end position
          int block_end = block_start + round(rexp1(r));
          if (block_end > contig_lengths[j]) {
            block_end = contig_lengths[j];
          }
          
          // add to output
          this_start.push_back(block_start);
          this_end.push_back(block_end);
          this_contig.push_back(j + 1);
          this_parent.push_back(parent_haplo_ID[i][parent_index]);
          
          // swap parent
          parent_index = 1 - parent_index;
          
          // move block forward
          block_start = block_end + 1;
        }
        
      }  // end loop through contigs
      
      // fill in other output vectors
      int n_blocks = this_start.size();
      vector<int> this_time(n_blocks, time[i]);
      vector<int> this_child(n_blocks, child_haplo_ID[i]);
      
      push_back_multiple(output_time, this_time);
      push_back_multiple(output_child_haplo_ID, this_child);
      push_back_multiple(output_contig, this_contig);
      push_back_multiple(output_start, this_start);
      push_back_multiple(output_end, this_end);
      push_back_multiple(output_parent_haplo_ID, this_parent);
    }
    
  }  // end loop through haplotype tree
  
  // end timer
  if (!silent) {
    chrono_timer(t1);
  }
  
  return Rcpp::List::create(Rcpp::Named("time") = output_time,
                            Rcpp::Named("child_haplo_ID") = output_child_haplo_ID,
                            Rcpp::Named("contig") = output_contig,
                            Rcpp::Named("start") = output_start,
                            Rcpp::Named("end") = output_end,
                            Rcpp::Named("parent_haplo_ID") = output_parent_haplo_ID);
}
