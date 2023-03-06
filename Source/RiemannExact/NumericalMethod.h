#include "EulerEOS.h"

class NumericalMethod
{
public:

  NumericalMethod(double, int);

  Eigen::ArrayXXf PressureApprox(Eigen::ArrayXXf, Eigen::ArrayXXf);

  Eigen::ArrayXXf SoundSpeed(Eigen::ArrayXXf);

  Eigen::ArrayXXf InitialPressure(Eigen::ArrayXXf);

  std::array<double,2> PressureFlux(double, Eigen::ArrayXXf, double);

  Eigen::ArrayXXf NewtonRaphson(Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);

private:

  double gamma;
  int nCells;
  EulerEOS eos;

  
 
};
