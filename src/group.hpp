#ifndef CELLSPLIT_GROUP_HPP
#define CELLSPLIT_GROUP_HPP

#include <vector>
#include <string>

#include "pointers.hpp"

namespace CellSplit_NS {

class Group : protected Pointers {
public:
  Group(CellSplit *);
  
  void create_all();
  void create_group(const std::vector<std::string> &);

  // iterate from start_indices[i] to end_indices[i]
  //  - length of these two vectors will be 1 if atom group,
  //  and >= 1 if molecule group

  
  std::vector<int> start_indices; // first index in chunk
  std::vector<int> end_indices; // last index + 1 in chunk

  std::string name,style;

  inline static std::vector<std::string> NAMES;
  
private:
  // group properties, all groups must be in chunks

  void group_atoms(int,int);
  void group_molecule(int);

  
};


}
#endif
