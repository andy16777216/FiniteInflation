#include <iostream>
//#include <boost/math/special_functions/bessel.hpp>
#include <complex.h>
#include <math.h>
#include <complex_bessel.h>

int main(int argc, char *argv[])
{
  using namespace std;

  const double Pi = M_PI;

  cout << "Hello, world!" << endl;

  double kmin = 1.E-6;
  double kmax = 1.;

  //std::complex<double> D10(0,-sqrt(Pi/4));
  //std::complex<double> D20(0,-sqrt(Pi/4));


  //double D10 = 1.0*atof(argv[1]);
  //double D20 = 1.0e-3*atof(argv[2]);
  //double logkc = -1.0*atof(argv[1]);

//cout << boost::math::cyl_bessel_j(1, 3);
//cout << boost::math::cyl_neumann(3, 4);

  int order = 2;
  complex<double> z(1.0,0.0);
  cout << sp_bessel::besselJ(order, z) << endl;



return 0;
}
