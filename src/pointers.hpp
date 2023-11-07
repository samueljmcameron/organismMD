#ifndef CELLSPLIT_POINTERS_HPP
#define CELLSPLIT_POINTERS_HPP


#include <mpi.h>
#include "cellsplit.hpp"

namespace CellSplit_NS {
  

class Pointers {
public:
  
  Pointers(CellSplit *ptr) : cellsplit(ptr), input(ptr->input),commbrick(ptr->commbrick),
			     domain(ptr->domain),world(ptr->world),inputfile(ptr->inputfile),
			     atoms(ptr->atoms),groups(ptr->groups),computes(ptr->computes),
			     dumps(ptr->dumps),varmap(ptr->varmap),
			     fixes(ptr->fixes),integrate(ptr->integrate) {};
  
protected:


  CellSplit *cellsplit;

  
  std::unique_ptr<Input> &input;
  std::unique_ptr<Domain> &domain;
  std::unique_ptr<CommBrick> &commbrick;
  std::unique_ptr<Atom> &atoms;
  std::unique_ptr<Integrate> &integrate;
  std::vector<std::unique_ptr<Group>> &groups;
  std::vector<std::unique_ptr<Fix>> &fixes;
  std::vector<std::unique_ptr<Compute>> &computes;
  std::vector<std::unique_ptr<Dump>> &dumps;
  
  
  std::unique_ptr<std::fstream> &inputfile;
  std::unique_ptr<std::map<std::string,std::string>> &varmap;
  
  MPI_Comm &world;
};
}


#endif
