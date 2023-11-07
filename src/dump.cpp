#include <iomanip>
#include <iostream>
#include "dump.hpp"
#include "comm_brick.hpp"

#include "domain.hpp"
#include "fix.hpp"
#include "atom.hpp"
#include "compute.hpp"
#include "integrate.hpp"

using namespace CellSplit_NS;

/* dump file structure is going to be that each timestep is written to an individual (VTK)
   file, and then there is a .pvd file to display the timesteps collectively. Files are
   always going to be outputted in parallel. */
Dump::Dump(CellSplit *phafd) : Pointers(phafd) {
  precision = -1;

};

Dump::~Dump() {
}


void Dump::init(const std::vector<std::string> &v_line ) {


  name = v_line.at(0);

  for (auto &dumpname : NAMES)
    if (dumpname == name)
      throw std::runtime_error("Error: Duplicate of dump " + name + std::string("."));

  NAMES.push_back(name);
  
  dump_type = v_line.at(1);

  if (dump_type.rfind("atom",0) == 0) {
    which_atoms = dump_type.substr(5);
    dump_type = "atom";
  }

  if (dump_type == "atom") {
    fext = "vtp";
  } else if (dump_type == "lammpsatom") {
    fext = "lmmpstrj";
  } else
    throw std::runtime_error("Invalid dump type (must be atom).");
  
  base_name = v_line.at(2);
  collection_name = base_name + std::string(".pvd");

  
  nopath_base_name = base_name;
  
  size_t endslash = base_name.find_last_of("\\/");

  if (endslash != std::string::npos) {
    nopath_base_name = base_name.substr(endslash+1);
  }



  fmt = "lf";
  every = std::stoll(v_line.at(3));

  auto it = v_line.begin()+4;
  while (it != v_line.end()) {
    
    if (*it == "attributes") {
      ++it;
      int numattributes = std::stoi(*it);
      ++it;
      for (int i = 0; i < numattributes; i++) {
	if (it == v_line.end())
	  throw std::invalid_argument("Incorrect number of attributes in dump file.");
	attributes.push_back(*it);
	++it;
      }

      if (attributes.size() != numattributes)
	throw std::invalid_argument("Duplicate attribute values in dump file.");
    } else if (*it == "decimalplaces") {
      ++it;
      precision = std::stoi(*it);
      ++it;
      if (precision <= 0)
	throw std::invalid_argument("Dump precision must be greater than or equal to zero.");
    } else if (*it == "fmt") {
      ++it;
      fmt = *it;
      ++it;
      if (fmt != "lf" && fmt != "e")
	throw std::invalid_argument("Dump fmt must be either lf or e.");
    } else
      throw std::invalid_argument("Invalid dump command.");
  }

  if (attributes.size() == 0)
    throw std::invalid_argument("Need attributes to write dump files.");

  
}

void Dump::setup()
{

  for (auto &word : attributes) {
    if (word.rfind("c_",0) == 0) {
      
      
      std::string cid = word.substr(2);
      
      int index = 0;
      
      for (auto &name : Compute::NAMES) {

	if (cid == name) {
	  break;
	}
	index += 1;
      }
      if (index == computes.size())
	throw std::runtime_error("Cannot write dump, compute ID " + cid
				 + std::string("doesn't exist."));

      computes.at(index)->dump_callers.insert(name);
      
    } else if (word.rfind("f_",0) == 0) {
      
      std::string fid = word.substr(2);
      
      int index = 0;
      
      for (auto &name : Fix::NAMES) {
	if (fid == name) {
	  break;
	}
	index += 1;
      }
      
      
      if (index == fixes.size())
	throw std::runtime_error("Cannot write dump, fix ID " + fid
				 + std::string("doesn't exist."));

      fixes.at(index)->dump_callers.insert(name);
    }
  }

  if (dump_type=="lammpsatom") {
    create_lammps_style_file();
  }
}

void Dump::start_of_step()
{

  if (integrate->timestep % every == 0) {
    for (auto &compute : computes)
      if (compute->dump_callers.find(name) != compute->dump_callers.end())
	compute->this_step = true;
    
    for (auto & fix : fixes)
      if (fix->dump_callers.find(name) != fix->dump_callers.end())
	fix->this_step = true;
  }
}

void Dump::require_calculations()
{
  for (auto &compute : computes)
    if (compute->dump_callers.find(name) != compute->dump_callers.end())
      compute->this_step = true;
  
  for (auto & fix : fixes)
    if (fix->dump_callers.find(name) != fix->dump_callers.end())
      fix->this_step = true;

}
  
/* called at start of simulation */
void Dump::write_collection_header()
{

  if (commbrick->me == 0 && dump_type != "lammpsatom") {
    auto myfile = std::ofstream(collection_name);
    if (not myfile.is_open()) {
      throw std::runtime_error(std::string("Cannot open file ") + collection_name);
    }
    
    myfile << "<?xml version=\"1.0\"?>" << std::endl
	   << "<VTKFile type=\"Collection\"  version=\"1.0\""
	   << " byte_order=\"LittleEndian\">" << std::endl
	   << "<Collection>" << std::endl;
    
    
    myfile.close();
  }


}



/* called at end of simulation */
void Dump::write_collection_footer() 
{

  if (commbrick->me == 0 && dump_type != "lammpsatom") {
    auto myfile = std::fstream(collection_name, std::ios_base::app);
    if (not myfile.is_open()) {
      throw std::runtime_error(std::string("Cannot open file ") + collection_name);
    }
    
    myfile << "</Collection>" << std::endl
	   << "</VTKFile>";
    
    myfile.close();
  }
}


/* called whenever new file is dumped. */
void Dump::write_collection_middle()
{


  if (dump_type != "lammpsatom") {
    double time = integrate->timestep*integrate->dt;
    create_instance_name();

    if (commbrick->me == 0) {
      auto myfile = std::fstream(collection_name,std::ios_base::app);
      if (not myfile.is_open()) {
	throw std::runtime_error(std::string("Cannot open file ") + collection_name);
      }
      
      myfile << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\""
	     << " file=\"" << pnopath_instance_name << "\"/>" << std::endl;
      
      myfile.close();
    }
  }

  

  if (dump_type == "atom") {
    write_atom_timestep();
    write_atom_timestep_pvtp();
  } else if (dump_type == "lammpsatom") {
    write_lammps_dump();
  }
}





void Dump::create_instance_name()
{

  instance_name = base_name + std::string("_p") + std::to_string(commbrick->me)
    + std::string("_") + std::to_string(integrate->timestep) + std::string(".") + fext;
  
  pinstance_name = base_name +  std::string("_") + std::to_string(integrate->timestep)
    + std::string(".p") + fext;
  
  
  
  
  pnopath_instance_name = nopath_base_name + std::string("_")
    + std::to_string(integrate->timestep) + std::string(".p") + fext;
 

  
}

void Dump::create_lammps_style_file()
{
  instance_name = base_name + std::string("_p") + std::to_string(commbrick->me)
    +std::string(".") + fext;
  
  auto myfile = std::fstream(instance_name, std::ios::out);  

  myfile.close();


}



void Dump::write_atom_timestep()
/*============================================================================*/
/*
  Write atom data to a vtk (and paraview) compatible file
  (extension .vtp).

*/
/*============================================================================*/

{
  
  auto myfile = std::fstream(instance_name, std::ios::out);

  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + instance_name);
  }

  if (precision > 0 && fmt == "lf") {
    myfile << std::fixed << std::setprecision(precision);
  } else if (precision > 0 && fmt == "e") {
    myfile << std::scientific << std::setprecision(precision);
  } else if (fmt == "e") {
    myfile << std::scientific;
  }

  int numpoints;
  if (which_atoms == "all")
    numpoints = atoms->xs.cols();
  else if (which_atoms == "local")
    numpoints = atoms->nlocal;
  else if (which_atoms == "ghost")
    numpoints = atoms->nghost;
  else
    throw std::runtime_error("need atoms to be dumped to be all, local, or ghost");
    
  
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"PolyData\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<PolyData>" << std::endl
	 << "<Piece NumberOfPoints=\"" << numpoints
	 << "\" NumberOfLines=\"1\">" << std::endl;
  


  // first find and save the position data, as it must always be present in the atom dump.
  auto it = std::find(attributes.begin(),attributes.end(),"x");
  if (it == attributes.end())
    throw std::runtime_error("Atom dump must contain positions.");

  myfile << "<Points>" << std::endl;
  write_ascii_data(myfile,"x",atoms->xs,false);
  myfile << "</Points>" << std::endl;

  myfile << "<PointData>" << std::endl;


  // don't re-save position data
  for (auto &word : attributes)
    if (word != "x") 
      process_attribute_name(myfile,word);

  myfile << "</PointData>" << std::endl;
  myfile << "</Piece>" << std::endl
	 << "</PolyData>" << std::endl
	 << "</VTKFile>" << std::endl;    
  myfile.close();

  
}


void Dump::write_atom_timestep_pvtp()
/*============================================================================*/
/*
  Write atom data to a vtk (and paraview) compatible file
  (extension .pvtp).

*/
/*============================================================================*/

{

  if (commbrick->me == 0) {
    auto myfile = std::fstream(pinstance_name, std::ios::out);
    
    if (not myfile.is_open()) {
      throw std::runtime_error(std::string("Cannot open file ") + instance_name);
    }
    
    int numpoints;
    if (which_atoms == "all")
      numpoints = atoms->xs.cols();
    else if (which_atoms == "local")
      numpoints = atoms->nlocal;
    else if (which_atoms == "ghost")
      numpoints = atoms->nghost;
    else
      throw std::runtime_error("need atoms to be dumped to be all, local, or ghost");
    
    
    myfile << "<?xml version=\"1.0\"?>" << std::endl
	   << "<VTKFile type=\"PPolyData\"  version=\"1.0\""
	   << " byte_order=\"LittleEndian\">" << std::endl
	   << "<PPolyData GhostLevel=\"0\">" << std::endl
	   << "<PPoints>" << std::endl;
    
    write_ascii_data(myfile,"x",atoms->xs,true);
    myfile << "</PPoints>" << std::endl;
    
    myfile << "<PPointData>" << std::endl;

    for (auto &word : attributes)
      if (word != "x")
	process_attribute_name(myfile,word,true);
    
    myfile << "</PPointData>" << std::endl;
    for (int proc = 0; proc < commbrick->nprocs; proc++) {
      std::string sourcename =  nopath_base_name  + std::string("_p") + std::to_string(proc)
	+ std::string("_") + std::to_string(integrate->timestep) + std::string(".") + fext;

      myfile << "<Piece Source=\"" + sourcename + std::string("\"/>") << std::endl;
    }
    
    
    myfile << "</PPolyData>" << std::endl
	   << "</VTKFile>" << std::endl;    
    myfile.close();
    
  }
}


void Dump::write_lammps_dump()
/*============================================================================*/
/*
  Write atom data to a similar format to a lammps dump file.

*/
/*============================================================================*/

{
  
  auto myfile = std::fstream(instance_name, std::ios::out|std::ios::app);

  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + instance_name);
  }

  if (precision > 0 && fmt == "lf") {
    myfile << std::fixed << std::setprecision(precision);
  } else if (precision > 0 && fmt == "e") {
    myfile << std::scientific << std::setprecision(precision);
  } else if (fmt == "e") {
    myfile << std::scientific;
  }
  
  int numpoints = atoms->nlocal;
    


  myfile << "ITEM: TIMESTEP" << std::endl << integrate->timestep << std::endl;
  myfile << "ITEM: NUMBER OF ATOMS" << std::endl << numpoints << std::endl;
  myfile << "ITEM: BOX BOUNDS pp pp pp" << std::endl 
	 << domain->boxlo[0] << " " << domain->boxhi[0] << std::endl
	 << domain->boxlo[1] << " " << domain->boxhi[1] << std::endl
	 << domain->boxlo[2] << " " << domain->boxhi[2] << std::endl;

  myfile << "ITEM: ATOMS ";
  for (auto &word : attributes)
    if (word == "x")
      myfile << "x y z ";
    else if (word == "F")
      myfile << "Fx Fy Fz ";
    else if (word == "ux")
      myfile << "ux uy uz ";
    else {
      int nc = 1;
      if (word.rfind("c_",0) == 0) {
	
	
	std::string cid = word.substr(2);
	
	int index = 0;
	
	for (auto &name : Compute::NAMES) {
	  
	  if (cid == name) {
	    break;
	  }
	  index += 1;
	}
	
	if (index == computes.size())
	  throw std::runtime_error("Cannot write dump, compute ID " + cid
				   + std::string("doesn't exist."));
	
	nc = computes.at(index)->numberofcomponents;
	
      } else if (word.rfind("f_",0) == 0) {
	
	std::string fid = word.substr(2);
	
	int index = 0;
	
	for (auto &name : Fix::NAMES) {
	  if (fid == name) {
	    break;
	  }
	  index += 1;
	}
	
	
	if (index == fixes.size())
	  throw std::runtime_error("Cannot write dump, fix ID " + fid
				   + std::string("doesn't exist."));
	
	nc = fixes.at(index)->numberofcomponents;
	
      }

      if (nc > 1) 

	for (int component = 0; component < nc; component++)
	  myfile << word << "_" << component << " ";

      else
	myfile << word << " ";

    }

  myfile << std::endl;

  for (int point = 0; point < numpoints; point++) {
    
    for (auto &word : attributes) {
      
      if (word.rfind("c_",0) == 0) {
	
	
	std::string cid = word.substr(2);
	
	int index = 0;
	
	for (auto &name : Compute::NAMES) {
	  
	  if (cid == name) {
	    break;
	  }
	  index += 1;
	}
	
	if (index == computes.size())
	  throw std::runtime_error("Cannot write dump, compute ID " + cid
				   + std::string("doesn't exist."));
	
	int nc = computes.at(index)->numberofcomponents;
	for (int component = 0; component < nc; component ++ )
	  myfile << computes.at(index)->array[point*nc+component] << " ";
	
      } else if (word.rfind("f_",0) == 0) {
	
	std::string fid = word.substr(2);
	
	int index = 0;
	
	for (auto &name : Fix::NAMES) {
	  if (fid == name) {
	    break;
	  }
	  index += 1;
	}
	
	
	if (index == fixes.size())
	  throw std::runtime_error("Cannot write dump, fix ID " + fid
				   + std::string("doesn't exist."));
	
	int nc = fixes.at(index)->numberofcomponents;
	for (int component = 0; component < nc; component ++ )
	  myfile << fixes.at(index)->array[point*nc+component] << " ";
	
      } else {
	
	if (word == "x")
	  myfile << atoms->xs.col(point)(0) << " "
		 << atoms->xs.col(point)(1) << " " 
		 << atoms->xs.col(point)(2) << " " ;
	
	else if (word == "F")
	  myfile << atoms->Fs.col(point)(0) << " "
		 << atoms->Fs.col(point)(1) << " " 
		 << atoms->Fs.col(point)(2) << " " ;
	
	else if (word == "ux")
	  myfile << atoms->uxs.col(point)(0) << " "
		 << atoms->uxs.col(point)(1) << " " 
		 << atoms->uxs.col(point)(2) << " " ;
	
	
	else if (word == "ix")
	  myfile << atoms->images[point] << " ";
	
	else if (word == "ID")
	  myfile << atoms->tags[point] << " ";
	
	else if (word == "alive")
	  myfile << atoms->alive[point] << " ";
	
	else if (word == "type")
	  myfile << atoms->types[point] << " ";
	
	else if (word == "label")
	  myfile << atoms->labels[point] << " ";
	
	else 
	  throw std::runtime_error("Invalid per atom quantity in dump file.");

      }
    }
    myfile << std::endl;
  }

  myfile.close();

  
}


void Dump::process_attribute_name(std::fstream &myfile,const std::string &word,
				  bool for_pvtp)
{

  
  if (word.rfind("c_",0) == 0) {


    std::string cid = word.substr(2);
    
    int index = 0;

    for (auto &name : Compute::NAMES) {

      if (cid == name) {
	break;
      }
      index += 1;
    }

    if (index == computes.size())
      throw std::runtime_error("Cannot write dump, compute ID " + cid
			       + std::string("doesn't exist."));


    if (dump_type == "atom") {
      if (!computes.at(index)->per_atom)
	throw std::runtime_error("Cannot write dump, compute ID " + cid
				 + std::string("is not a per_atom quantity."));
      write_ascii_data(myfile,word,computes.at(index)->array,
		       computes.at(index)->numberofcomponents,for_pvtp);
    } else
      throw std::runtime_error("Something wrong, should not get here. ");

    
    
  } else if (word.rfind("f_",0) == 0) {

    std::string fid = word.substr(2);

    int index = 0;
    
    for (auto &name : Fix::NAMES) {
      if (fid == name) {
	break;
      }
      index += 1;
    }

    
    if (index == fixes.size())
      throw std::runtime_error("Cannot write dump, fix ID " + fid
			       + std::string("doesn't exist."));



    if (dump_type == "atom") {
      if (!fixes.at(index)->per_atom)
	throw std::runtime_error("Cannot write dump, fix ID " + fid
				 + std::string("is not a per_atom quantity."));
      
      write_ascii_data(myfile,word,fixes.at(index)->array,
		       fixes.at(index)->numberofcomponents,for_pvtp);
      
    } else
      throw std::runtime_error("Something wrong, should not get here. ");

    


    
  } else if (dump_type == "atom") {

    if (word == "F")
      write_ascii_data(myfile,word,atoms->Fs,for_pvtp);
    else if (word == "ux")
      write_ascii_data(myfile,word,atoms->uxs,for_pvtp);
    else if (word == "ix")
      write_ascii_data(myfile,word,atoms->images,1,for_pvtp);
    else if (word == "ID")
      write_ascii_data(myfile,word,atoms->tags,1,for_pvtp);
    else if (word == "alive")
      write_ascii_data(myfile,word,atoms->alive,1,for_pvtp);
    else if (word == "type")
      write_ascii_data(myfile,word,atoms->types,1,for_pvtp);
    else if (word == "label")
      write_ascii_data(myfile,word,atoms->labels,1,for_pvtp);
    
    else
      throw std::runtime_error("Invalid per atom quantity in dump file.");

  } else {
    throw std::runtime_error("Something wrong, should not get here. ");
  }
  return;
}




void Dump::write_ascii_data(std::fstream &myfile,const std::string &arrname,
			    Eigen::Ref<Eigen::Matrix3Xd> array,
			    bool for_pvtp)
{

  if (!for_pvtp) {
    myfile << "<DataArray Name=\"" << arrname << "\" "
	   << "type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    
    if (which_atoms == "all") {
      for (int i = 0; i < array.cols(); i++) {
	myfile << array.col(i)(0) << " " << array.col(i)(1) << " "
	       << array.col(i)(2) << " ";
	
      }
    } else if (which_atoms == "local") {
      for (int i = 0; i < atoms->nlocal; i++) {
	if (atoms->labels[i] == CellSplit_NS::Atom::OWNED)
	  myfile << array.col(i)(0) << " " << array.col(i)(1) << " "
		 << array.col(i)(2) << " ";
	else
	  throw std::runtime_error("Somehow an un-owned atom snuck into the owned? shouldn't happen.");
      }
    } else if (which_atoms == "ghost") {
      for (int i = atoms->nlocal; i < atoms->nlocal+atoms->ngathered; i++) {
	if (atoms->labels[i] == CellSplit_NS::Atom::GHOST)
	  myfile << array.col(i)(0) << " " << array.col(i)(1) << " "
		 << array.col(i)(2) << " ";
      } 
    }
    
    myfile << std::endl << "</DataArray>" << std::endl;
  } else {
    myfile << "<PDataArray Name=\"" << arrname << "\" "
	   << "type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>" << std::endl;

  }
}


template <typename T>
void Dump::write_ascii_data(std::fstream &myfile,const std::string &arrname,
			    std::vector<T> array,int numberofcomponents,
			    bool for_pvtp)
{
  std::string dtype;
  if (typeid(T) == typeid(double))
    dtype = "Float64";
  else if (typeid(T) == typeid(int))
    dtype = "Int64";
  else
    throw std::runtime_error("Cannot use non double non int type for writing ascii data.");
  

  if (!for_pvtp) {
    myfile << "<DataArray Name=\"" << arrname << "\" "
	   << "type=\"" << dtype << "\" NumberOfComponents=\""
	   << numberofcomponents << "\" format=\"ascii\">" << std::endl;
    
    
    if (which_atoms == "all") {
      for (int i = 0; i < array.size(); i++) {
	myfile << array[i] << " ";
	
      }
    } else if (which_atoms == "local") {
      
      int index = 0;
      
      for (int i = 0; i < atoms->nlocal; i++) {
	if (atoms->labels[i] == CellSplit_NS::Atom::OWNED)
	  for (int components = 0; components < numberofcomponents; components++)
	    myfile << array[index++] << " ";
	else
	  throw std::runtime_error("Somehow an un-owned atom snuck into the owned? shouldn't happen.");
	
      }
    } else if (which_atoms == "ghost") {
      
      int index = 0;
      
      for (int i = atoms->nlocal; i < atoms->nlocal+atoms->ngathered; i++) {
	if (atoms->labels[i] == CellSplit_NS::Atom::GHOST)
	  for (int components = 0; components < numberofcomponents; components++)
	    myfile << array[index++] << " ";
	else
	  index += numberofcomponents;
	
      } 
    }
    
    myfile << std::endl << "</DataArray>" << std::endl;
  } else {
    myfile << "<PDataArray Name=\"" << arrname << "\" "
	   << "type=\"" << dtype << "\" NumberOfComponents=\""
	   << numberofcomponents << "\" format=\"ascii\"/>" << std::endl;

  }
  

}
