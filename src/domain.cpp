
#include "domain.hpp"
#include "atom.hpp"
#include "comm_brick.hpp"

#define IMGMASK 1023
#define IMGMAX 512
#define IMGBITS 10
#define IMG2BITS 20


using namespace CellSplit_NS;

Domain::Domain(CellSplit *cellsplit)
  : Pointers(cellsplit),dimensions(-1)
{

  boxset = false;
  subboxset = false;
  
  
};

Domain::~Domain() = default;


void Domain::set_box(const std::vector<std::string> &v_line)
{

  if (dimensions == -1)
    throw std::invalid_argument("Error: Need to set dimension before setting box.");
  
  int iargs = 0;





  bool period_flag = false;
  bool origin_flag = false;
  
  while (iargs < v_line.size()) {
    
    if (v_line[iargs] == "boxdims") {
      period_flag = true;
      iargs ++;
      period[0] = std::stod(v_line[iargs++]);
      period[1] = std::stod(v_line[iargs++]);
      period[2] = std::stod(v_line[iargs++]);
    } else if (v_line[iargs] == "boxorigin") {
      origin_flag = true;
      iargs ++;
      boxlo[0] = std::stod(v_line[iargs++]);
      boxlo[1] = std::stod(v_line[iargs++]);
      boxlo[2] = std::stod(v_line[iargs++]);
      
    } else {
      throw std::invalid_argument("Bad arguments for domain_setup.");
      
    }
  }

  if (!period_flag)
    throw std::runtime_error("need to set box dimensions.");

  
  if (!origin_flag)
    throw std::runtime_error("need to set box origin.");

  
  for (int i = 0; i < 3; i++)  
    boxhi[i] = boxlo[i] + period[i];

  boxset = true;
  set_subbox();
}

void Domain::set_subbox()
{

  if (commbrick->nprocs == 1) {
    sublo[0] = boxlo[0];
    sublo[1] = boxlo[1];
    sublo[2] = boxlo[2];
  
    subhi[0] = boxhi[0];
    subhi[1] = boxhi[1];
    subhi[2] = boxhi[2];

  } else if (commbrick->nprocs < 4) {

    double interval = period[0]/commbrick->nprocs;
    
    sublo[0] = boxlo[0] + commbrick->me*interval;
    sublo[1] = boxlo[1];
    sublo[2] = boxlo[2];
  
    subhi[0] = boxlo[0] + (commbrick->me+1)*interval;
    subhi[1] = boxhi[1];
    subhi[2] = boxhi[2];


  } else if (commbrick->nprocs == 4) {

    double x_interval = period[0]/2;
    double y_interval = period[1]/2;
    int xint = commbrick->me % 2;
    int yint = commbrick->me/2;

    sublo[0] = boxlo[0] + xint*x_interval;
    sublo[1] = boxlo[1] + yint*y_interval;
    sublo[2] = boxlo[2];
  
    subhi[0] = boxlo[0] + (xint+1)*x_interval;
    subhi[1] = boxlo[1] + (yint+1)*y_interval;
    subhi[2] = boxhi[2];

  } else {
    throw std::runtime_error("Can currently only handle max 4 processors.");
  }
  
  subboxset = true;
  return;
}


void Domain::pbc() const
{

  int idim,otherdims;
  for (int i = 0; i < atoms->nlocal; i++) {
    if (atoms->xs(0,i) < boxlo[0]) {
      atoms->xs(0,i) += period[0];
      idim = atoms->images[i] & IMGMASK;
      otherdims = atoms->images[i] ^ idim;
      idim --;
      idim &= IMGMASK;
      atoms->images[i] = otherdims | idim;
    }
    if (atoms->xs(0,i) >= boxhi[0]) {
      atoms->xs(0,i) -= period[0];
      atoms->xs(0,i) = (atoms->xs(0,i) > boxlo[0] ? atoms->xs(0,i) : boxlo[0]);
      idim = atoms->images[i] & IMGMASK;
      otherdims = atoms->images[i] ^ idim;
      idim ++;
      idim &= IMGMASK;
      atoms->images[i] = otherdims | idim;
    }
    if (atoms->xs(1,i) < boxlo[1]) {
      atoms->xs(1,i) += period[1];
      idim = (atoms->images[i] >> IMGBITS) & IMGMASK;
      otherdims = atoms->images[i] ^ (idim << IMGBITS);
      idim--;
      idim &= IMGMASK;
      atoms->images[i] = otherdims | (idim << IMGBITS);
    }
    if (atoms->xs(1,i) >= boxhi[1]) {
      atoms->xs(1,i) -= period[1];
      atoms->xs(1,i) = (atoms->xs(1,i) > boxlo[1] ? atoms->xs(1,i) : boxlo[1]);
      idim = (atoms->images[i]  >> IMGBITS) & IMGMASK;
      otherdims = atoms->images[i] ^ (idim << IMGBITS);
      idim ++;
      idim &= IMGMASK;
      atoms->images[i] = otherdims | (idim << IMGBITS);
    }
    if (atoms->xs(2,i) < boxlo[2]) {
      atoms->xs(2,i) += period[2];
      idim = atoms->images[i] >> IMG2BITS;
      otherdims = atoms->images[i] ^ (idim << IMG2BITS);
      idim--;
      idim &= IMGMASK;
      atoms->images[i] = otherdims | (idim << IMG2BITS);
    }
    if (atoms->xs(2,i) >= boxhi[2]) {
      atoms->xs(2,i) -= period[2];
      atoms->xs(2,i) = (atoms->xs(2,i) > boxlo[2] ? atoms->xs(2,i) : boxlo[2]);
      idim = atoms->images[i]  >> IMG2BITS;
      otherdims = atoms->images[i] ^ (idim << IMG2BITS);
      idim ++;
      idim &= IMGMASK;
      atoms->images[i] = otherdims | (idim << IMG2BITS);
    }
    // store atom in unwrapped coords (needed for polymer updates).
    unmap(atoms->uxs.col(i),atoms->xs.col(i),atoms->images[i]);
  }

  return;
}


void Domain::unmap(Eigen::Ref<Eigen::Vector3d> unwrapped,
		   const Eigen::Ref<const Eigen::Vector3d> & wrapped,
		   int image) const
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  unwrapped(0) = wrapped(0) + xbox*period[0];
  unwrapped(1) = wrapped(1) + ybox*period[1];
  unwrapped(2) = wrapped(2) + zbox*period[2];

}


void Domain::map(Eigen::Ref<Eigen::Vector3d> wrapped,
		 const Eigen::Ref<const Eigen::Vector3d> & unwrapped,
		 int image) const
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  wrapped(0) = unwrapped(0) - xbox*period[0];
  wrapped(1) = unwrapped(1) - ybox*period[1];
  wrapped(2) = unwrapped(2) - zbox*period[2];

}


int Domain::set_image() const
{
  int i1 = IMGMAX << IMG2BITS;
  int i2 = IMGMAX << IMGBITS;
  int i3 = IMGMAX;

  return i1 | i2 | i3;

}
