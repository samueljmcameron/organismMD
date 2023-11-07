#ifndef CELLSPLIT_CELLSPLIT_HPP
#define CELLSPLIT_CELLSPLIT_HPP


#include <vector>
#include <memory>
#include <mpi.h>
#include <fstream>
#include <map>

namespace CellSplit_NS {


class CellSplit
{
  void create();
public:
  
  CellSplit(MPI_Comm,int, char **);
  ~CellSplit();

  std::unique_ptr<class Atom> atoms;
  std::unique_ptr<class Domain> domain;
  std::unique_ptr<class Input> input;
  std::unique_ptr<class CommBrick> commbrick;  
  std::unique_ptr<class Integrate> integrate;
  std::vector<std::unique_ptr<class Group>> groups;
  std::vector<std::unique_ptr<class Fix>> fixes;
  std::vector<std::unique_ptr<class Compute>> computes;
  std::vector<std::unique_ptr<class Dump>> dumps;
  
  
  MPI_Comm world;
  
  std::unique_ptr<std::fstream> inputfile;
  
  std::unique_ptr<std::map<std::string,std::string>> varmap;
  
};
}
#endif
