#ifndef CELLSPLIT_INPUT_HPP
#define CELLSPLIT_INPUT_HPP

#include "pointers.hpp"

namespace CellSplit_NS {

class Input : protected Pointers {
public:
  Input(CellSplit *);

  ~Input();

  void read();


};
  
}
#endif
