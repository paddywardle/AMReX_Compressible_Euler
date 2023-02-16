#include "EulerEOS.h"
#include "Limiters.h"
#include <vector>

class NumericalMethod
{
public:

  NumericalMethod(double, int);

  double reconstruction_uL(double, double, double, Limiters);

  double reconstruction_uR(double, double, double, Limiters);

  var_arr FORCE_flux(var_arr, var_arr, var_arr, var_arr, double, double);

  var_arr HLLC_flux(var_arr, var_arr, var_arr, var_arr, double, double);

  var_arr uL_half_update(var_arr, var_arr, double, double);

  var_arr uR_half_update(var_arr, var_arr, double, double);

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

  double lax_friedrich_flux(double, double, double, double, double, double);

  double richtmyer_flux(double, double, double, double, double, double);

  var_arr wavespeed(var_arr, var_arr, var_arr, var_arr);

  var_arr uHLLC(var_arr, var_arr, double, double);

  var_arr fHLLC(var_arr, var_arr, var_arr, var_arr, var_arr, var_arr, double, double, double);
 
};
