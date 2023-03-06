#include "EulerEOS.h"

EulerEOS::EulerEOS(double gamma)
  :gamma(gamma){};

var_array EulerEOS::prim_to_con(var_array u_p, int dim)
{
  var_array u_c;

  u_c[0] = u_p[0];
  u_c[1] = u_p[0] * u_p[1];
  u_c[3] = u_p[0] * u_c[3];

  if (dim == 2)
    {
      u_c[2] = u_p[2] / (gamma - 1.0) + 0.5 * u_p[0] * (pow(u_p[1], 2) + pow(u_p[3], 2));
    }
  else if (dim == 1)
    {
      u_c[2] = u_p[2] / (gamma - 1.0) + 0.5 * u_p[0] * pow(u_p[1], 2);
    }

  return u_c;
}

var_array EulerEOS::con_to_prim(var_array u_c, int dim)
{
  var_array u_p;

  u_p[0] = u_c[0];
  u_p[1] = u_c[1] / u_c[0];
  u_p[3] = u_c[3] / u_c[0];
  
  if (dim == 2)
    {
      u_p[2] = (u_c[2] - 0.5 * u_c[0] * (pow(u_p[1], 2.0) + pow(u_p[3], 2.0))) * (gamma - 1.0);
    }
  else if (dim == 1)
    {
      u_p[2] = (u_c[2] - 0.5 * u_c[0] * pow(u_p[1], 2.0)) * (gamma - 1.0);
    }   
  
  return u_p;
}
  
var_array EulerEOS::Euler_flux_fn(var_array f, var_array f_prim)
{
  
  var_array flux_fn;
  double rho = f[0];
  double energy = f[2];
  double vx = f_prim[1];
  double vy = f_prim[3];
  double pressure = f_prim[2];
      
  flux_fn[0] = flux_fn_rho(rho, vx);
  flux_fn[1] = flux_fn_mom(rho, vx, pressure);
  flux_fn[2] = flux_fn_E(energy, vx, pressure);
  flux_fn[3] = rho * vx * vy;

  return flux_fn;

}

var_array EulerEOS::Euler_flux_fn_Y(var_array f, var_array f_prim)
{
  
  var_array flux_fn;
  double rho = f[0];
  double energy = f[2];
  double vx = f_prim[1];
  double vy = f_prim[3];
  double pressure = f_prim[2];
      
  flux_fn[0] = flux_fn_rho(rho, vy);
  flux_fn[1] = rho * vx * vy;
  flux_fn[2] = flux_fn_E(energy, vy, pressure);
  flux_fn[3] = flux_fn_mom(rho, vy, pressure);

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
