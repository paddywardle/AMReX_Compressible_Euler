#include "common.h"

class EulerEOS
{
public:
  
  EulerEOS(double);

  var_array prim_to_con(var_array);

  var_array con_to_prim(var_array);
  
  double flux_fn_rho(double, double);

  double flux_fn_mom(double, double, double);

  double flux_fn_E(double, double, double);

  var_array Euler_flux_fn(var_array, var_array);

private:

  double gamma;
    
};
