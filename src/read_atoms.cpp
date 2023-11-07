#include <iostream>
#include <algorithm>
#include <set>

#include "utility.hpp"
#include "atom.hpp"

#include "read_atoms.hpp"
 

#include "comm_brick.hpp"

using namespace CellSplit_NS;

ReadAtoms::ReadAtoms(CellSplit *cellsplit) : Pointers(cellsplit) {};

int ReadAtoms::read_file(const std::string & fname_in)
{


  fname = fname_in;
  int errflag = SUCCESS;
  int total_errflag;

  // open the file to be read
  
  datafile.open(fname);
  
  if (datafile.fail()) {
    std::cerr << "could not open file " << fname << std::endl;
    errflag = NOFILE;
  }

  if (errflag != SUCCESS) errflag = 1;
  else errflag = 0;

  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_MAX,world);
  
  if (total_errflag) return NOFILE;

  // create the atoms from the file info
  
  errflag = create_atoms(atoms->xs.cols());
  
  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);

  if (total_errflag) return FORMAT_ERROR;
  
  datafile.close();
  datafile.clear();
  
  return SUCCESS;
  
}



int ReadAtoms::create_atoms(int totalatoms)
/* Construct atoms from the lines of the data file. */
{

  // allocate unique_ptr for atoms



  int iatom = 0;


  std::string line;

  std::vector<std::string> v_line;

  int created_atoms;

  int errflag = 0;
  
  while (std::getline(datafile,line)) {
      
    if (line == "" || line[0] == '#') continue;

    v_line = utility::split_line(line);
    

    created_atoms = atoms->add_atom(v_line,iatom,errflag);
    if (errflag) break;
	  
    iatom += created_atoms;
    
  }
  int totalerr;
  MPI_Allreduce(&errflag,&totalerr,1,MPI_INT,MPI_MAX,world);
  if (totalerr)
    throw std::runtime_error("Failed to create atoms.");

  
  std::cout << "read all atoms on processor " << commbrick->me << std::endl;

  if (iatom != totalatoms) {
    std::cerr << "Incorrect number of atoms specified in " << fname << std::endl;
    std::cerr << "atoms in file = " << iatom << ", total atoms = " << totalatoms << std::endl;
    return FORMAT_ERROR;
  }

  return SUCCESS;
}
