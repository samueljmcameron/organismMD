#ifndef CellSplit_DOMAIN_HPP
#define CellSplit_DOMAIN_HPP


#include <array>
#include <Eigen/Core>

#include "pointers.hpp"

namespace CellSplit_NS {

  
class Domain : protected Pointers
{
public:
  Domain(CellSplit *);
  ~Domain();
  
  void set_box(const std::vector<std::string> &);

  void pbc () const;
  int set_image() const;

  bool boxset, subboxset;


  int dimensions;
  
  std::array<double,3> period,boxlo,boxhi;
  std::array<double,3> sublo,subhi;

  void map(Eigen::Ref<Eigen::Vector3d>,
   	   const Eigen::Ref<const Eigen::Vector3d> &,int ) const;

  
private:
  void set_subbox();
  void unmap(Eigen::Ref<Eigen::Vector3d>,
  	     const Eigen::Ref<const Eigen::Vector3d> &,int image) const;

};
}
#endif
