
#include "Mosquito_pop.hpp"
#include "Host_pop.hpp"

using namespace std;

//------------------------------------------------
// initialise
void Mosquito_pop::init(Parameters &params_, int M_) {
  
  // store pointer to parameters
  params = &params_;
  
  // counts of mosquito types
  M = M_;
  Sv = M;
  Ev = 0;
  Iv = 0;
  
  // ring buffer of deaths in latent stage
  ringbuffer_index = 0;
  ringbuffer = vector<int>(params->v);
  
  // vector of mosquitoes
  memory_increment = 100;
  mosq_vec_n = 0;
  mosq_vec = vector<Mosquito>(memory_increment);
  
}

//------------------------------------------------
// step ringbuffer forwards and apply all scheduled deaths
void Mosquito_pop::apply_ringbuffer_death() {
  
  // step ringbuffer index forwards
  ringbuffer_index++;
  if (ringbuffer_index > (params->v - 1)) {
    ringbuffer_index = 0;
  }
  
  // move deaths from Ev to Sv
  Sv += ringbuffer[ringbuffer_index];
  Ev -= ringbuffer[ringbuffer_index];
  
  // reset buffer ready to hold new values
  ringbuffer[ringbuffer_index] = 0;
}

//------------------------------------------------
// loop through all mosquitoes and apply life events
void Mosquito_pop::update_pop(int t) {
  
  int i = 0;
  while (i < mosq_vec_n) {
    
    // transition from Ev to Iv
    if (!mosq_vec[i].infectious_on) {
      if (mosq_vec[i].time_infectious == t) {
        mosq_vec[i].Ev_to_Iv();
        Ev--;
        Iv++;
      }
    }
    
    // death. Rather than popping from the vector, swap the last mosquito in at
    // this position then decrease the maximum vector index by one. This
    // achieves the same result but avoids reallocating memory. Note that deaths
    // will always be in the Iv state, because mosquitoes that die in the latent
    // state are not present in this vector
    if (mosq_vec[i].time_death == t) {
      mosq_vec[i].init_from_mosq(mosq_vec[mosq_vec_n - 1]);
      mosq_vec_n--;
      Iv--;
      Sv++;
    } else {
      ++i;
    }
    
  }
}

//------------------------------------------------
// for every host, draw all mosquitoes who are infected by this individual (only
// if host is in active sexual state). Split into infections that do/don't make
// it through EIP and update all objects as needed
void Mosquito_pop::draw_new_infections(Host_pop &host_pop, int t, int k) {
  
  // TODO - more efficient to draw non-surviving in aggregate? i.e. get total
  // infectivity over all hosts, then draw from binomial taking account prob
  // survival and distribute these values over ringbuffer. The second binomial
  // draw of the *remaining* Sv for those surviving. Would need to walk back
  // through hosts and allocate these in proportion to infectivity. Could do
  // this using the inverse CDF method, i.e. make N sorted uniform draws from 0
  // to total_infectivity, walk through pop taking cumulative sum of infectivity
  // and infect whenever draw falls inside this section
  
  // loop through all hosts in this deme
  for (int i2 = 0; i2 < host_pop.H[k]; ++i2) {
    int this_host = host_pop.host_index[k][i2];
    
    // get infectivity of this host, and calculate the probability of a
    // single susceptible mosquito being infected by this single host
    double infectivity = host_pop.host_vec[this_host].get_infectivity(t);
    if (infectivity == 0) {
      continue;
    }
    double prob_viable_infection = params->a * infectivity / double(host_pop.H[k]);
    
    // binomial draw of number of mosquitoes infected by this host, and decrease
    // Sv accordingly
    int new_infections = rbinom1(Sv, prob_viable_infection);
    if (new_infections == 0) {
      continue;
    }
    Sv -= new_infections;
    Ev += new_infections;
    
    // draw whether mosquitoes die in Ev state, and if so add to ringbuffer
    int surviving_infections = new_infections;
    int j = ringbuffer_index;
    for (int i = 0; i < params->v; ++i) {
      j++;
      if (j > (params->v - 1)) {
        j = 0;
      }
      int n_die_today = rbinom1(surviving_infections, 1 - params->p);
      ringbuffer[j] += n_die_today;
      surviving_infections -= n_die_today;
      if (surviving_infections == 0) {
        break;
      }
    }
    if (surviving_infections == 0) {
      continue;
    }
    
    // increase size of mosq_vec if needed
    if ((mosq_vec_n + surviving_infections) > mosq_vec.size()) {
      mosq_vec.resize(mosq_vec.size() + memory_increment);
    }
    
    // surviving infections are added to mosq_vec, and time of death is drawn
    // from the end of the EIP onwards
    for (int i = 0; i < surviving_infections; ++i) {
      int t_death = t + params->v + rgeom1(1.0 - params->p) + 1;
      mosq_vec[mosq_vec_n + i].init_from_host(*params, t, t_death, host_pop.host_vec[this_host]);
    }
    mosq_vec_n += surviving_infections;
  }
}

//------------------------------------------------
// print summary
void Mosquito_pop::print_summary() {
  
  // count number of Ev and Iv individuals in mosq_vec
  int Ev_in_mosq_vec = 0;
  int Iv_in_mosq_vec = 0;
  for (int i = 0; i < mosq_vec_n; ++i) {
    if (mosq_vec[i].infectious_on) {
      Iv_in_mosq_vec++;
    } else {
      Ev_in_mosq_vec++;
    }
  }
  
  // print values
  print("Sv:", Sv, ", Ev:", Ev, ", Iv", Iv);
  print("ringbuffer_index:", ringbuffer_index);
  print("ringbuffer:");
  print_vector(ringbuffer);
  print("sum(ringbuffer):", sum(ringbuffer));
  print("mosq_vec_n:", mosq_vec_n);
  print("mosq_vec.size():", mosq_vec.size());
  print("Ev in mosq_vec:", Ev_in_mosq_vec);
  print("Iv in mosq_vec:", Iv_in_mosq_vec);
  /*
  for (int i = 0; i < mosq_vec_n; ++i) {
    std::cout << mosq_vec[i].time_death << " ";
  }
  print("");
  */
  // check that values correspond
  bool pass_checks = true;
  if ((Sv + Ev + Iv) != M) {
    cpp11::stop("Error: (Sv + Ev + Iv) != M");
  }
  if ((sum(ringbuffer) + Ev_in_mosq_vec) != Ev) {
    cpp11::stop("Error: (sum(ringbuffer) + Ev_in_mosq_vec) != Ev");
  }
  if (Iv_in_mosq_vec != Iv) {
    cpp11::stop("Error: Iv_in_mosq_vec != Iv");
  }
  if (mosq_vec_n > mosq_vec.size()) {
    cpp11::stop("mosq_vec_n > mosq_vec.size()");
  }
  if (pass_checks) {
    print("passed all checks");
  }
}
