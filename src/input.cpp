#include "input.hpp"

#include "utility.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "atom.hpp"
#include "integrate.hpp"
#include "fix_organism_birthdeath.hpp"
#include "group.hpp"
#include "dump.hpp"

#include <string>
#include <set>
#include <chrono>
#include <iostream>

using namespace CellSplit_NS;

Input::Input(CellSplit *cellsplit) : Pointers(cellsplit) {}

Input::~Input() = default;

void Input::read()
{


  if (commbrick->me == 0) {
    std::cout << "Reading input file..." << std::endl << std::endl << std::endl;
  }


  
  std::vector<std::string> v_line;
  std::string line, firstword;
  
  while(std::getline(*inputfile,line)) {
    
    utility::convertVariables(line,*varmap);
    utility::replacePercentages(line,commbrick->me);

    v_line = utility::split_line(line);

    if (v_line.size() == 0) continue;
    
    firstword = v_line.at(0);
    v_line.erase(v_line.begin());

    
    if (firstword == "dimensions")
      domain->dimensions = std::stoi(v_line.at(0));

    else if (firstword == "domain")
      domain->set_box(v_line);
    else if (firstword == "atom_style") 
      atoms->setup(v_line);
    else if (firstword == "atom_populate") {

      std::cout << "got here on processor " << commbrick->me << std::endl;
      if  (!atoms->atomset)
	throw std::invalid_argument("CAnnot use command atom_populate before atom_style is set.");
      atoms->populate(v_line);
      std::cout << "populated on processor " << commbrick->me << std::endl;
      // call this to ensure everything is filled in that needs to be filled in.
      atoms->check_tags_and_types();
      groups.push_back(std::make_unique<Group>(cellsplit));
      groups.back()->create_all();
    } else if (firstword == "fix") {

      firstword = v_line.at(0);
      v_line.erase(v_line.begin());
      if (firstword == "organism/birthdeath") {
	if (atoms->ntypes <= 0)
	  throw std::runtime_error("Fix requires atoms to be created.");
	
	
	fixes.push_back(std::make_unique<FixOrganismBirthDeath>(cellsplit));
      
      }
      fixes.back()->init(v_line);
    } else if (firstword == "dt")
      integrate->dt = std::stod(v_line.at(0));
    else if (firstword == "dump") {
      dumps.push_back(std::make_unique<Dump>(cellsplit));
      dumps.back()->init(v_line);
    }
    else if (firstword == "run") {
      if (integrate->dt <= 0.0)
	throw std::runtime_error("Cannot run simulation without setting dt > 0.");

      integrate->nsteps = std::stoll(v_line.at(0));
      integrate->setup();
      integrate->run();
    }

  }

  /*
  MPI_Barrier(world);
  for (int proc = 0; proc < commbrick->nprocs; proc++) {
    if (proc == commbrick->me) {
      std::cout << "Processor " << proc << " domains:" <<  std::endl;
      for (int i = 0; i < 3; i++)
	std::cout << "(" << domain->sublo[i] << "," << domain->subhi[i] << ")" << std::endl;
      std::cout << "Processor " << proc << " atoms:" <<  std::endl; 
      for (int n = 0; n < atoms->xs.cols(); n++) {
	std::cout << "( ";
	for (int i = 0; i < atoms->xs.rows(); i++) 
	  std::cout << atoms->xs(i,n) << " ";
	std::cout << ")" << std::endl;
      }
    }
    MPI_Barrier(world);
  }
  */
  
  return;
}
