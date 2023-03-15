
#include "sim_indiv.hpp"

using namespace std;

//------------------------------------------------
void sim_indiv(cpp11::list args_progress) {
  
  //vector<double> x = cpp11::as_cpp<vector<double>>(args_list["x"]);
  //print(x[0], x[1]);
  
  // copy update_progress function from package
  cpp11::function update_progress = cpp11::package("SIMPLEGEN")["update_progress"];
  
  
  double z = 100.0;
  for (int i = 0; i < 1e8; ++i) {
    
    if ((i % (int)1e6) == 0) {
      update_progress(args_progress, "pb_sim", i, 1e8, true);
    }
    
    z = pow(z, 0.9);
  }
  print(z);
  
  //update_progress(args_progress, "pb_sim", 0, 100, true);
  
  //cpp11::stop("sub_stop");
}
