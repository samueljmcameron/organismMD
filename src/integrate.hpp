#ifndef CELLSPLIT_INTEGRATE_HPP
#define CELLSPLIT_INTEGRATE_HPP
#include <cstdint>
#include "pointers.hpp"

namespace CellSplit_NS {

class Integrate : protected Pointers {
public:
  Integrate(CellSplit *);

  ~Integrate();

  void setup();
  void run();

  int64_t nsteps; // total number of steps
  int64_t firststep; // first step
  int64_t timestep; // current step
  double dt; // time interval between steps

};
  
}
#endif
