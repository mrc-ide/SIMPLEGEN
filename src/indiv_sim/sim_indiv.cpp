
#include "sim_indiv.hpp"

using namespace std;

//------------------------------------------------
void sim_indiv(cpp11::list args_progress) {
  
  //vector<double> x = cpp11::as_cpp<vector<double>>(args_list["x"]);
  //print(x[0], x[1]);
  
  cpp11::function update_progress = cpp11::package("SIMPLEGEN")["update_progress"];
  update_progress(args_progress, "pb_sim", 0, 100, true);
  
  //cpp11::stop("sub_stop");
}
