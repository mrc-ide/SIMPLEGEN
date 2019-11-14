
#pragma once

#include <vector>

//------------------------------------------------
// class defining mosquito
class Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  std::vector<int> inoc_ID;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito() {};
  Mosquito(std::vector<int> &inoc_ID);
  
};
