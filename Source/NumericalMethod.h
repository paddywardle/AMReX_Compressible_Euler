#include "EulerEOS.h"
#include "Limiters.h"
#include <vector>

class NumericalMethod
{
public:

  NumericalMethod(double);

  double reconstruction_uL(double, double, double, Limiters);

  double reconstruction_uR(double, double, double, Limiters);

  var_array HLLC_flux(var_array, var_array, var_array, var_array, double, double, int);

  var_array HLLC_flux_Y(var_array, var_array, var_array, var_array, double, double, int);  

  var_array uL_half_update(var_array, var_array, double, double, int);

  var_array uR_half_update(var_array, var_array, double, double, int);

  var_array uL_half_updateY(var_array, var_array, double, double, int);

  var_array uR_half_updateY(var_array, var_array, double, double, int);

private:

  double gamma;
  int nCells;
  EulerEOS eos;

  double deltai_func(double, double, double, double);

  double slope_limiter(double, double, double);

  double minbee(double);

  double superbee(double);

  double van_leer(double);

  double van_albada(double);

  double limiterXi(double, Limiters);

  var_array wavespeed_x(var_array, var_array, var_array, var_array);

  var_array wavespeed_y(var_array, var_array, var_array, var_array);
  
  var_array uHLLC(var_array, var_array, double, double);

  var_array uHLLC_Y(var_array, var_array, double, double);

  var_array fHLLC(var_array, var_array, var_array, var_array, var_array, var_array, double, double, double);
 
};
