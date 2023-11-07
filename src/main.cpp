
#include <iostream>
#include <mpi.h>
#include "cellsplit.hpp"
#include "input.hpp"


int main(int argc, char **argv)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  int ierr = MPI_Init(NULL,NULL);

  {
    auto cellsplit = new CellSplit_NS::CellSplit(comm,argc,argv);
    cellsplit->input->read();

    
    delete cellsplit;
  }
  
  
  ierr = MPI_Finalize();

  return 0;
    
}
