
#include "relatedness_Haplo.h"
#include "misc_v12.h"

using namespace std;

//------------------------------------------------
// constructor
relatedness_Haplo::relatedness_Haplo(int haplo_ID) {
  initialise(haplo_ID);
}

//------------------------------------------------
// initialise (descended from itself at generation 0)
void relatedness_Haplo::initialise(int haplo_ID) {
  descendant_IDs = vector<int>(1, haplo_ID);
  generations = vector<int>(1, 0);
  recombinations = vector<int>(1, 0);
}

//------------------------------------------------
// merge in values from another haplo
void relatedness_Haplo::merge(relatedness_Haplo m) {
  push_back_multiple(descendant_IDs, m.descendant_IDs);
  push_back_multiple(generations, m.generations);
  push_back_multiple(recombinations, m.recombinations);
}

//------------------------------------------------
// increment generation counter for all descendents
void relatedness_Haplo::increment_generations() {
  for (int i = 0; i < generations.size(); ++i) {
    generations[i]++;
  }
}

//------------------------------------------------
// increment recombination counter for all descendents
void relatedness_Haplo::increment_recombinations() {
  for (int i = 0; i < recombinations.size(); ++i) {
    recombinations[i]++;
  }
}

//------------------------------------------------
// print status
void relatedness_Haplo::print_status() {
  print("descendant_IDs:");
  print_vector(descendant_IDs);
  print("generations:");
  print_vector(generations);
  print("recombinations:");
  print_vector(recombinations);
}
