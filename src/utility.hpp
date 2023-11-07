#ifndef CELLSPLIT_UTILITY_HPP
#define CELLSPLIT_UTILITY_HPP

#include <cstdint>
#include <mpi.h>
#include <string>
#include <vector>
#include <map>

namespace CellSplit_NS {

  namespace utility {
    std::vector<std::string> split_line(std::string&);
    void replacePercentages(std::string &, int);
    
    void convertVariables(std::string &,
			  std::map<std::string, std::string> const&);

  
    void check_MPI_duplicates(const std::vector<int> &,MPI_Comm ,int ,int ,
			      std::string);

    int make_unique_seed(int,const MPI_Comm &, int, int);

  }

union ubuf {
  double d;
  int64_t i;
  ubuf(const double &arg) : d(arg) {}
  ubuf(const int64_t &arg) : i(arg) {}
  ubuf(const int &arg) : i(arg) {}
};
}
#endif
