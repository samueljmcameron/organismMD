#ifndef CELLSPLIT_COMPUTE_HPP
#define CELLSPLIT_COMPUTE_HPP


#include <set>

#include "pointers.hpp"

namespace CellSplit_NS {

class Compute : protected Pointers
{
public:
  Compute(CellSplit *);

  std::vector<double> array;
  std::string name;

  virtual void init(const std::vector<std::string> &);
  virtual void in_fourier() = 0;
  virtual void end_of_step() = 0;

  void start_of_step();
  
  inline static std::vector<std::string> NAMES;
  bool per_atom;
  bool scalar;
  bool vector;


  bool this_step;
  std::set<std::string> dump_callers;

  int numberofcomponents; // number of components in the array (e.g. 1 for scalar, 3 for vector, etc.)
  
};

}
#endif
