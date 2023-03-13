#include "common.h"

class NumericalMethod
{
public:
  
  NumericalMethod(double, int);

  Eigen::ArrayXXf SamplePattern(double, double, double, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);

  Eigen::ArrayXXf leftShock(double, double , double , Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);

  Eigen::ArrayXXf leftRar(double, double , double, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);

  Eigen::ArrayXXf rightShock(double, double , double , Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);

  Eigen::ArrayXXf rightRar(double, double , double, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);
  
  double PressureApprox(Eigen::ArrayXXf, Eigen::ArrayXXf);

  Eigen::ArrayXXf SoundSpeed(Eigen::ArrayXXf, Eigen::ArrayXXf);

  std::array<double,2> PressureFlux(double, Eigen::ArrayXXf, double);

  std::array<double,2> NewtonRaphson(double, Eigen::ArrayXXf, Eigen::ArrayXXf, Eigen::ArrayXXf);

private:
  double gamma;
  int nCells;
 
};
