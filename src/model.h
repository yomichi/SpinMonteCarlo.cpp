#ifndef MODEL_H
#define MODEL_H

namespace potts{

template <class Lattice, class RNG01>
class Model{
public:
  Model(int L):lat_(L){}
  virtual ~Model(){}

  virtual double mag() const=0;
  virtual double ene() const=0;

  virtual void update(double beta, RNG01 &rng)=0;
  virtual void local_update(double beta, RNG01 &rng);
  virtual void SW_update(double beta, RNG01 &rng);
protected:
  Lattice lat_;
};

} // end of namespace potts

#endif
