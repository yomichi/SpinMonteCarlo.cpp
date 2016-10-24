#include "potts.h"
#include "ising.h"
#include "util/observable.hpp"
#include "util/squarelattice.hpp"
#include <iostream>
#include <boost/random.hpp>
#include <ctime>
#include <mpi.h>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

int main(int argc, char **argv)
{
  typedef util::SquareLattice Lattice;
  typedef boost::variate_generator<boost::mt19937, boost::uniform_real<> > RNG01;

  RNG01 rnd(boost::mt19937(static_cast<uint32_t>(std::time(0))), boost::uniform_real<>(0.0, 1.0));

  using namespace boost;
  using namespace boost::program_options;
  options_description opt("options");
  opt.add_options()
    ("help", "show this message")
    ("q,q", value<int>()->default_value(2), "number of state of a spin")
    ("L,L", value<int>()->default_value(10), "length of lattice")
    ("beta-min", value<double>()->default_value(0.1), "minimum beta")
    ("beta-max", value<double>()->default_value(1.0), "maximum beta")
    ("beta-num", value<int>()->default_value(10), "number of beta")
    ("thermalization", value<int>()->default_value(8), "number of exchange to thermalize")
    ("mcs", value<int>()->default_value(64), "number of exchange to measure")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, opt), vm);
  notify(vm);

  if( vm.count("help") ){
    std::cout << opt << std::endl;
    return 0;
  }

  const int q = vm["q"].as<int>();
  const int L = vm["L"].as<int>();
  const double bmin = vm["beta-min"].as<double>();
  const double bmax = vm["beta-max"].as<double>();
  const int nbeta = vm["beta-num"].as<int>();
  const int therm = vm["thermalization"].as<int>();
  const int MCS = vm["mcs"].as<int>();

  const double dbeta = (bmax-bmin)/(nbeta-1);

  potts::Potts<Lattice, RNG01> model(q, L);

  for(int ibeta=0; ibeta<nbeta; ++ibeta){
    const double beta = bmin + dbeta*ibeta;
    util::Observable obs_time;
    util::Observable obs_speed;
    
    for(int mcs=0; mcs < therm+MCS; ++mcs){
      boost::timer::cpu_timer tm;
      model.SW_update(beta, rnd);
      const double sec = tm.elapsed().wall * 1.0e-9;
      const double speed = 1.0/sec;
      obs_time << sec;
      obs_speed << speed;
    }
    std::cout << beta
       << " " << obs_time.mean() << " " << obs_time.error()
       << " " << obs_speed.mean() << " " << obs_speed.error()
       << std::endl;
  }

  return 0;
}
