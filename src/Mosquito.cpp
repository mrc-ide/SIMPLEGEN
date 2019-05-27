
#include "Mosquito.h"
#include "misc_v6.h"
#include "probability_v2.h"

using namespace std;

//------------------------------------------------
// constructor for mosquito class
Mosquito::Mosquito(int host_ID, int infection_time) :
  host_ID(host_ID),
  infection_time(infection_time) {}
