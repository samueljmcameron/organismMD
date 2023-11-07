#include "integrate.hpp"

#include "domain.hpp"
#include "comm_brick.hpp"
#include "atom.hpp"
#include "fix.hpp"
#include "dump.hpp"
#include "compute.hpp"

#include <mpi.h>
#include <iostream>
#include <string>
#include <set>
#include <chrono>

using namespace CellSplit_NS;

Integrate::Integrate(CellSplit *cellsplit) : Pointers(cellsplit) {
  nsteps = 0;
  firststep = 0;
  timestep = 0;
  dt = 0.0;
  
}

Integrate::~Integrate() = default;

void Integrate::setup()
{  

  if (atoms->ntypes >= 0) {
    // need to apply pbc to atoms in case they are outside of simulation domain
    domain->pbc();
    
    atoms->Fs.setZero();
  }
  
  int errflag = 0;
  int total_errflag;
  std::string ewhat = "";


  for (auto &dump : dumps) {
    try {
      dump->setup();
    } catch (const std::runtime_error &e) {
      errflag = 1;
      ewhat = e.what();
    }
    MPI_Allreduce(&errflag, &total_errflag,1,MPI_INT,MPI_SUM,world);

    if (total_errflag) {
      throw std::runtime_error(ewhat.c_str());
    }


    dump->write_collection_header();
  }


  for (auto &fix : fixes)
    fix->reset_dt();

  for (auto &fix : fixes)
    fix->setup();

  // note that this will dump out zero vectors for almost everything as nothing has been
  // computed yet
  for (auto &dump :dumps)
    if (timestep % dump->every == 0) {
      if (commbrick->me == 0)
	std::cout << "Saving on step " << timestep << std::endl;
      
      
      dump->write_collection_middle();
      
      
    }
  
}

void Integrate::run()
{

  if (commbrick->me == 0) 
    std::cout << "Running simulation of solution." << std::endl;
  
  
  int errflag = 0;
  int total_errflag = 0;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


  for (int i = 0; i < nsteps; i++) {
  
    timestep ++;


    /*
    for (auto &compute: computes)
      compute->start_of_step();
    */
    for (auto &fix: fixes)
      fix->start_of_step();

    for (auto &dump : dumps)
      dump->start_of_step();
    
    
    domain->pbc();
    
    atoms->Fs.setZero();

    for (auto &fix : fixes)
      fix->initial_integrate();

    for (auto &fix: fixes)
      fix->post_force();


    /*
    for (auto &compute : computes)
      compute->end_of_step();
    
    */
    for (auto &fix : fixes)
      fix->end_of_step();
    
    for (auto &dump :dumps)
      if (timestep % dump->every == 0) {
	if (commbrick->me == 0)
	  std::cout << "Saving on step " << timestep << std::endl;


	dump->write_collection_middle();

      
      }
    
    
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  for (auto &dump :dumps)
    dump->write_collection_footer();
  
  if (commbrick->me == 0) {
    std::cout << "Run time = "
	      << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	      << "seconds." << std::endl;  
  }


}

