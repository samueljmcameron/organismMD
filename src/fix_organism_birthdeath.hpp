
#ifndef CELLSPLIT_FIXATOM_DRAG_HPP
#define CELLSPLIT_FIXATOM_DRAG_HPP

#include <random>

#include "fix.hpp"

#include <memory>

namespace CellSplit_NS {

class FixOrganismBirthDeath : public Fix {
public:
  FixOrganismBirthDeath(CellSplit *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;
  virtual void initial_integrate() override {};
  virtual void post_force() override;
  virtual void pre_final_integrate() override {};
  virtual void final_integrate() override {};
  virtual void post_final_integrate() override {};

  virtual void end_of_step() override {};

private:

  std::uniform_real_distribution<> dist;
  std::mt19937 gen;

  std::vector<int> local_alive_list;
  bool birth_from_dead(int);

  void delete_dead_atoms();

  void procreate(int , int );
  void count_vitals();

  double birthrate, deathrate,shift;
  int localalive,localdead,nalive,ndead;
  int cleanevery;
};

}

#endif
