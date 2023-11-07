#ifndef CELLSPLIT_DUMP_HPP
#define CELLSPLIT_DUMP_HPP

#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Core>
#include <cstdint>

#include "pointers.hpp"

namespace CellSplit_NS {

class Dump : protected Pointers
{
public:
  Dump(CellSplit *);
  ~Dump();

  void init(const std::vector<std::string> &);

  void setup();
  void start_of_step();
  void require_calculations();
  void write_collection_header();
  void write_collection_footer();
  void write_collection_middle();
  
  inline static std::vector<std::string> NAMES;

  int64_t every;

  std::string name;
  

protected:
  std::string base_name; // base name of file without extension, e.g. vtkfiles/field_p0
  std::string nopath_base_name; // name of file without path or extension, e.g. field_p0
  std::string collection_name; // collection name, .e.g vtkfiles/field_p0.pvd
  std::string instance_name; // base name plus time and extension, e.g. vtkfiles/field_p0_0.vti
  std::string pinstance_name; // base name plus time and extension, e.g. vtkfiles/field_0.pvti
  std::string pnopath_instance_name; // above but without path, e.g. field_0.pvti
  std::string fext; // file extension of per time files (e.g. ".vtp" for atom)
  
  std::string dump_type; // either atom, or lammpsatom


  std::vector<std::string> attributes;

private:
  void create_instance_name();

  void write_atom_timestep();
  void write_grid_timestep();

  void write_atom_timestep_pvtp();
  void write_grid_timestep_pvti();
  
  void create_lammps_style_file();
  void write_lammps_dump();
  void process_attribute_name(std::fstream &, const std::string &,
			      bool for_pvtp=false);
  void write_ascii_data(std::fstream &, const std::string &,
			Eigen::Ref<Eigen::Matrix3Xd>,bool for_pvtp);

  template <typename T>
  void write_ascii_data(std::fstream &, const std::string &,
			std::vector<T>,int,bool for_pvtp);
  
  std::string which_atoms;

  unsigned int bytelength;

  int precision;
  std::string fmt;
  
};

}
#endif
