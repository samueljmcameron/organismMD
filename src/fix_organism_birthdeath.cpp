#include <algorithm>
#include <iostream>
#include <random>

#include "utility.hpp"

#include "atom.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "integrate.hpp"

#include "fix_organism_birthdeath.hpp"

using namespace CellSplit_NS;

FixOrganismBirthDeath::FixOrganismBirthDeath(CellSplit *cellsplit)
  : Fix(cellsplit), dist(0,1), gen(0) {
};



void FixOrganismBirthDeath::init(const std::vector<std::string> &v_line)
{

  Fix::init(v_line);

  find_group(v_line.at(1));


  // v_line contains a seed in the argument, so need to extract that and initialise
  //  a new RNG with it.
  // need to reconvert to a single string (removing the fix name and the group name)
  std::string line = "";

  std::vector<std::string> new_v_line(v_line);

  int seed = std::stoi(new_v_line.at(2));

  seed = utility::make_unique_seed(seed,world,commbrick->me,commbrick->nprocs);
  
  gen.seed(seed);

  birthrate = std::stod(new_v_line.at(3));
  deathrate = std::stod(new_v_line.at(4));
  shift = std::stod(new_v_line.at(5));
  cleanevery = std::stoi(new_v_line.at(6));



  

}

  
void FixOrganismBirthDeath::setup()
{

}


void FixOrganismBirthDeath::post_force()
{
  delete_dead_atoms();
  double dt = integrate->dt;

  count_vitals();  // create local_alive_list and compute nalive, ndead

  double birthrn, deathrn;
  int i;



  std::vector<int> new_atoms;

  for (int ia = 0; ia < localalive; ia++) {
    i = local_alive_list[ia];
    
    birthrn = dist(gen);

    deathrn = dist(gen);
    
    if (birthrn  < birthrate*dt &&  ! (deathrn < deathrate*nalive*dt)) {

      bool procreated = birth_from_dead(i);

      if (! procreated) {
	new_atoms.push_back(i);
      }
      
      
      
    } else if (deathrn < deathrate*nalive*dt && ! (birthrn  < birthrate*dt)) {

      atoms->alive[i] = 0;

      
    }
    
  }

  if (new_atoms.size() + atoms->nlocal > atoms->xs.cols())
    atoms->resize(atoms->nlocal+new_atoms.size());


  int jatom = atoms->nlocal;
  for (auto iatom : new_atoms) {
    atoms->create_atom(atoms->types[iatom],&atoms->xs(0,iatom));
    procreate(iatom,jatom++);
  }
  if (jatom != atoms->nlocal)
    std::cout << "somtthing is gone wrong?" << std::endl;

  atoms->tag_extend();
  atoms->check_tags_and_types();

  
}



void FixOrganismBirthDeath::delete_dead_atoms()
{

  if (integrate->timestep % cleanevery != 0)
    return;

  int nlocal = atoms->nlocal;

  std::vector<int> dlist(nlocal);
  
  for (int i = 0; i < nlocal; i++) {
    if (! atoms->alive[i])
      dlist[i] = 1;
    else
      dlist[i] = 0;
  }

  int i = 0;

  while (i < nlocal) {
    if (dlist[i]) {
      atoms->copy(nlocal - 1, i);
      dlist[i] = dlist[nlocal - 1];
      nlocal--;
    } else
      i++;
  }

  int natoms_previous = atoms->nlocal;
  atoms->nlocal = nlocal;
  
  int ndelete = natoms_previous - atoms->nlocal;

}



/* Atom i must be alive */
bool FixOrganismBirthDeath::birth_from_dead(int i)
{

  bool procreated = false;
  
  for (int j = 0; j < atoms->nlocal; j++) {

    if (i == j) continue;

    if (!atoms->alive[j]) {
      procreate(i,j);
      procreated = true;
      break;
    }
  }


  return procreated;

}

void FixOrganismBirthDeath::procreate(int i, int j)
{
  double polar,azim;
  double cx,cy,cz,dx,dy,dz;

  cx = atoms->xs(0,i);
  cy = atoms->xs(1,i);
  cz = atoms->xs(2,i);
  
  if (domain->dimensions == 3) {
    polar = M_PI*dist(gen);
    azim = 2*M_PI*dist(gen);
    dx = shift*cos(azim)*sin(polar)/2.;
    dy = shift*sin(azim)*sin(polar)/2.;
    dz = shift*cos(polar)/2;
    
  } else {
    azim = 2*M_PI*dist(gen);
    dx = shift*cos(azim)/2.;
    dy = shift*sin(azim)/2.;
    dz = 0.0;
  }

  atoms->alive[j] = 1;

  atoms->xs(0,i) = cx-dx;
  atoms->xs(1,i) = cy-dy;
  atoms->xs(2,i) = cz-dz;

  atoms->xs(0,j) = cx+dx;
  atoms->xs(1,j) = cy+dy;
  atoms->xs(2,j) = cz+dz;
  
  return;
  
}


/* ---------------------------------------------------------------------- */
/* Store atoms which are alive in local_alive_list and count total
   number of alive atoms and dead atoms across all processors. */
/* ---------------------------------------------------------------------- */

void FixOrganismBirthDeath::count_vitals()
{
  
  localalive = 0;
  localdead = 0;
  if (local_alive_list.size() < atoms->nlocal) {
    local_alive_list.resize(atoms->nlocal);
  }
  
  
  for (int i = 0; i < atoms->nlocal; i++) {
    if (atoms->alive[i]) {
      local_alive_list[localalive++] = i;
    } else {
      localdead++;
    }
  }

  MPI_Allreduce(&localalive, &nalive, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&localdead, &ndead, 1, MPI_INT, MPI_SUM, world);

  return;
}
