#include "EulerEOS.h"

EulerEOS::EulerEOS(double gamma)
  :gamma(gamma){};

var_array EulerEOS::prim_to_con(var_array u_p)
{
  var_array u_c;

  u_c[0] = u_p[0];
  u_c[1] = u_p[0] * u_p[1];
  u_c[2] = u_p[2] / (gamma - 1.0) + 0.5 * u_p[0] * pow(u_p[1], 2);

  return u_c;
}

var_array EulerEOS::con_to_prim(var_array u_c)
{
  var_array u_p;

  u_p[0] = u_c[0];
  u_p[1] = u_c[1] / u_c[0];
  u_p[2] = (u_c[2] - 0.5 * u_c[0] * pow(u_p[1], 2.0)) * (gamma - 1.0);
  
  return u_p;
}
var_array EulerEOS::Euler_flux_fn(var_array f, var_array f_prim)
{
  
  var_array flux_fn;
  double rho = f[0];
  double energy = f[2];
  double velocity = f_prim[1];
  double pressure = f_prim[2];
      
  flux_fn[0] = flux_fn_rho(rho, velocity);
  flux_fn[1] = flux_fn_mom(rho, velocity, pressure);
  flux_fn[2] = flux_fn_E(energy, velocity, pressure);

  return flux_fn;

}

double EulerEOS::flux_fn_rho(double rho, double v)
{
  // flux function for density
  return rho * v;
}

double EulerEOS::flux_fn_mom(double rho, double v, double p)
{
  // flux function for momentum
  return rho * pow(v, 2.0) + p;
}

double EulerEOS::flux_fn_E(double E, double v, double p)
{
  // flux function for total energy
  return (E + p) * v;
}
