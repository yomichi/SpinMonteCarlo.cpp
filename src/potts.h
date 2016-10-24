#ifndef POTTS_H
#define POTTS_H

#include <vector>
#include <cmath>
#include "model.h"
#include "util/union_find.hpp"

#include <boost/math/special_functions/expm1.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace potts{

template <class Lattice, class RNG01>
class Potts : public Model<Lattice, RNG01>{
public:
  Potts(int q, int L):
    Model<Lattice, RNG01>(L),
    q_(q),
    nsites_(this->lat_.num_sites()),nbonds_(this->lat_.num_bonds()),
    NoverQ_(1.0*nsites_/q_),
    spins_(nsites_,0),
    nzero_(nsites_),
    ene_(-nbonds_){}
  
  double mag() const{ return nzero_ - NoverQ_;}
  double ene() const{ return ene_;}

  void update(double beta, RNG01 &rng){ SW_update(beta, rng);};
  void local_update(double beta, RNG01 &rng);
  void SW_update(double beta, RNG01 &rng);

private:
  int q_;
  int nsites_;
  int nbonds_;
  double NoverQ_;

  std::vector<int> spins_;
  int nzero_; // number of s=0 spins
  double ene_;
};

template <class Lattice, class RNG01>
void Potts<Lattice, RNG01>::local_update(double beta, RNG01 &rng)
{
  for(int site = 0; site < nsites_; ++site){
    const int now  = spins_[site];
    const int cand = q_ * rng();
    int de = 0;
    foreach(int n, this->lat_.neighbors(site)){
      de += now  == spins_[n] ? 1 : 0;
      de -= cand == spins_[n] ? 1 : 0;
    }
    const double p = std::exp(-beta*de);
    if(rng() < p){
      spins_[site] = cand;
      ene_ += de;
      nzero_ -= now == 0 ? 1:0;
      nzero_ += cand == 0 ? 1:0;
    }
  }
}

template <class Lattice, class RNG01>
void Potts<Lattice, RNG01>::SW_update(double beta, RNG01 &rng)
{
  const double p = -boost::math::expm1(-beta);
  union_find::Nodes nodes(nsites_);

  for(int bond = 0; bond < nbonds_; ++bond){
    const int src = this->lat_.source(bond);
    const int tgt = this->lat_.source(bond);
    if(spins_[src] == spins_[tgt] && rng() < p){
      union_find::unify(nodes, src, tgt);
    }
  }
  const int nc = union_find::clusterize(nodes);
  std::vector<int> cluster_size(nc);
  std::vector<int> cluster_spin(nc);
  for(int cl=0; cl<nc; ++cl){
    cluster_spin[cl] = q_ * rng();
  }
  for(int site=0; site<nsites_; ++site){
    const int id = union_find::cluster_id(nodes, site);
    ++cluster_size[id];
    spins_[site] = cluster_spin[id];
  }
}

} // end of namespace potts

#endif
