#include <iostream>
#include <algorithm>
#include <set>

#include "atom.hpp"
#include "utility.hpp"
//#include "read_dump.hpp"
#include "read_atoms.hpp"
#include "domain.hpp"
#include "comm_brick.hpp"
#include "group.hpp"


using namespace CellSplit_NS;

Atom::Atom(CellSplit *cellsplit)
  : Pointers(cellsplit),size_forward(3), size_reverse(3),size_border(6) {
  ntypes = -1;
  atomset = false;

};


void Atom::setup(const std::vector<std::string> & v_line) {

  organism_flag = false;

  std::vector<int> Natomsperproc(commbrick->nprocs,0);

  
  int nargs = v_line.size();  
  int iarg = 0;

  while (v_line.at(iarg) != "natoms") {
    if (v_line.at(iarg) == "organism") {
      organism_flag = true;

    } else
      throw std::runtime_error("Unrecognised atom style.");


    iarg += 1;
    
  }

  iarg += 1;

  std::string tmpstr;

  int proc = 0;

  // get number of atoms per processor
  while (iarg < nargs) {

    tmpstr = v_line.at(iarg);
    std::string::size_type pos = tmpstr.find("*");
    if (pos != std::string::npos) {
      int atoms_per_proc = std::stoi(tmpstr.substr(0,pos));
      int nps = std::stoi(tmpstr.substr(pos+1));

      if (proc+nps > commbrick->nprocs)
	throw std::runtime_error("More atom groupings than there are processors.");
      for (int i = 0; i < nps; i++) {
	Natomsperproc.at(proc++) = atoms_per_proc;

      }
    } else
      Natomsperproc.at(proc++) = std::stoi(tmpstr);
    iarg += 1;
  }


  int Natoms = Natomsperproc.at(commbrick->me);

  ngathered = 0;
  nlocal = 0;
  nghost = 0;

  
  resize(Natoms);
  xs.setZero();
  Fs.setZero();
  uxs.setZero();
  
  atomset = true;
  
}

void Atom::resize(int Natoms) {

  xs.conservativeResize(3,Natoms);
  uxs.conservativeResize(3,Natoms);
  Fs.conservativeResize(3,Natoms);

  tags.resize(Natoms);
  types.resize(Natoms);
  labels.resize(Natoms);
  images.resize(Natoms);

  if (organism_flag) {
    alive.resize(Natoms);
  }

}

void Atom::populate(const std::vector<std::string> & v_line) {

  if (v_line.at(0) == "read_atoms") {
     ReadAtoms read_atoms(cellsplit);
      
     int errflag = read_atoms.read_file(v_line.at(1));
     if (errflag != ReadAtoms::SUCCESS) errflag = 1;
     else errflag = 0;
     int total_errflag;
     MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,world);
     if (total_errflag)
       throw std::runtime_error("Could not read data atom file. ");

  }/* else if (v_line.at(0) == "read_dump") {

    ReadDump read_dump(cellsplit);

    std::vector<std::string> new_v_line(v_line);

    new_v_line.at(0) = "atom";
    read_dump.init(new_v_line);
    read_dump.process_attributes();

    } */
  else
    throw std::runtime_error("Invalid atom_populate command.");

}

void Atom::check_tags_and_types()
{

  
  std::set<int> tset(tags.begin(),std::next(tags.begin(), nlocal));


  int errflag = 0;
  int totalerr;

  /*
  if (tags.size() != nlocal) {
    std::cerr << "number of atom allocated on processor " << commbrick->me
	      << " doesn't match number of atoms created." << std::endl;
    errflag = 1;
  }
  */
  
  MPI_Allreduce(&errflag,&totalerr,1,MPI_INT,MPI_MAX,world);
  if (totalerr)
    throw std::runtime_error("Incorrect atoms allocated on processor(s).");



  if (tset.size() != nlocal) {
    std::cerr << "duplicate IDs on processor " << commbrick->me << "." << std::endl;
    errflag = 1;
  }

  MPI_Allreduce(&errflag,&totalerr,1,MPI_INT,MPI_MAX,world);
  if (totalerr)
    throw std::runtime_error("duplicate IDs found on processor(s).");



  // check if any IDs are the same across all processors

  std::vector<int> tmptags(tags.begin(),std::next(tags.begin(), nlocal));
  utility::check_MPI_duplicates(tmptags,world,commbrick->me,commbrick->nprocs,"IDs");



  
  
  int typemax;

  if (types.size() == 0)
    typemax = -1;
  else
    typemax = *(std::max_element(types.begin(),types.end()));



  int all_typemax; 

  MPI_Allreduce(&typemax,&all_typemax,1,MPI_INT,MPI_MAX,world);

  ntypes = all_typemax + 1;




}


void Atom::tag_extend()
{
  int maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = (maxtag > tags[i] ? maxtag : tags[i]);
  int maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_INT,MPI_MAX,world);

  int notag = 0;
  for (int i = 0; i < nlocal; i++) if (tags[i] == -1) notag ++;


  int notag_sum;
  MPI_Scan(&notag,&notag_sum,1,MPI_INT,MPI_SUM,world);

  int itag = maxtag_all + notag_sum - notag+1;
  for (int i = 0; i < nlocal; i++) if (tags[i] == -1) tags[i] = itag++;

}

void Atom::copy(int i, int j)
{

  tags[j] = -1;
  if (organism_flag) 
    alive[j] = alive[i];


  types[j] = types[i];
  xs(0,j) = xs(0,i);
  xs(1,j) = xs(1,i);
  xs(2,j) = xs(2,i);

  labels[j] = labels[i];
  images[j] = images[i];

  return;
  
}

void Atom::create_atom(int itype,double *coord) {

  int iatom = nlocal;
  tags[iatom] = -1;
  if (organism_flag) 
    alive[iatom] = 0;
  types[iatom] = itype;
  xs(0,iatom) = coord[0];
  xs(1,iatom) = coord[1];
  xs(2,iatom) = coord[2];

  labels[iatom] = Atom::OWNED;
  images[iatom] = domain->set_image();

  nlocal += 1;

}

int Atom::add_atom(std::vector<std::string> v_line, int iatom, int &errflag)
{
  

  int l_index = 0;

  tags[iatom] = std::stoi(v_line.at(l_index++));

  if (organism_flag) {
    std::string word = v_line.at(l_index++);
    if (word == "alive")
      alive[iatom] = 1;
    else if (word == "dead")
      alive[iatom] = 0;
    else {
      std::cerr << "need to be dead or alive " << std::endl;
      errflag = 1;
    }
  }

  
  types[iatom] = std::stoi(v_line.at(l_index++));
  xs(0,iatom) = std::stod(v_line.at(l_index++));
  xs(1,iatom) = std::stod(v_line.at(l_index++));
  xs(2,iatom) = std::stod(v_line.at(l_index++));
  
  labels[iatom] = Atom::OWNED;
  images[iatom] = domain->set_image();

  nlocal += 1;
  
  return 1;
  

}

/*
int Atom::unpack_reverse(int n,const std::vector<int> & sendlist, double *buf)
{

  int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    Fs(0,j) += buf[m++];
    Fs(1,j) += buf[m++];
    Fs(2,j) += buf[m++];
  }

  return m;
  
}

int Atom::pack_comm(int n,const std::vector<int> & sendlist, double *buf,
		    const std::array<double,3> & pbc)
{

  int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j) + pbc[0];
    buf[m++] = xs(1,j) + pbc[1];
    buf[m++] = xs(2,j) + pbc[2];
  }

  return m;
  
}

int Atom::pack_comm_in_z(int n,const std::vector<int> & sendlist, double * buf,
			 const std::vector<double> &pbcz )
{
 int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j);
    buf[m++] = xs(1,j);
    buf[m++] = xs(2,j) + pbcz[i];
  }
  return m;

}

int Atom::pack_border_in_z(int nsend, const std::vector<int> & sendlist,double *buf,
			   const std::vector<double> & pbcz,
			   const std::vector<int> & labelz)
{

  int j;

  int m = 0;

  for (int i = 0; i < nsend; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j);
    buf[m++] = xs(1,j);
    buf[m++] = xs(2,j) + pbcz[i];
    buf[m++] = ubuf(tags[j]).d;
    buf[m++] = ubuf(types[j]).d;
    buf[m++] = ubuf(labelz[i]).d;

  }

  return m;
}


int Atom::pack_border(int nsend, const std::vector<int> & sendlist,double *buf,
		      const std::array<double,3> & pbc)
{

  int j;

  int m = 0;

  for (int i = 0; i < nsend; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j) + pbc[0];
    buf[m++] = xs(1,j) + pbc[1];
    buf[m++] = xs(2,j) + pbc[2];
    buf[m++] = ubuf(tags[j]).d;
    buf[m++] = ubuf(types[j]).d;
    buf[m++] = ubuf(Atom::GHOST).d;
    
  }

  return m;
}



void Atom::unpack_border(int nrecv,int first, const double *buf)
{
  int m = 0;
  int last = first + nrecv;

  xs.conservativeResize(Eigen::NoChange,last);
  Fs.conservativeResize(Eigen::NoChange,last);
  tags.resize(last);
  types.resize(last);
  labels.resize(last);

  
  for (int i = first; i < last; i++) {

    xs(0,i) = buf[m++];
    xs(1,i) = buf[m++];
    xs(2,i) = buf[m++];
    tags[i] = (int) ubuf(buf[m++]).i;
    types[i] = (int) ubuf(buf[m++]).i;
    labels[i] = (int) ubuf(buf[m++]).i;
    
  }

  
  return;
}
*/
