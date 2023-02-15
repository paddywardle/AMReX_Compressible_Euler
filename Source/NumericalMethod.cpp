#include "NumericalMethod.h"

NumericalMethod::NumericalMethod(double gamma, int nCells)
  :gamma(gamma), nCells(nCells), eos(gamma){};

var_arr NumericalMethod::wavespeed(var_arr uL, var_arr uR, var_arr uL_prim, var_arr uR_prim)
{
  var_arr wavespeeds(1, 3);

  double csL = sqrt((gamma*uL_prim(2))/uL(0));
  double csR = sqrt((gamma*uR_prim(2))/uR(0));

  wavespeeds.col(0) = uL_prim(1) - csL;
  wavespeeds.col(1) = uR_prim(1) + csR;
  wavespeeds.col(2) = (uR_prim(2) - uL_prim(2) + uL(0)*uL_prim(1)*(wavespeeds(0) - uL_prim(1)) - uR(0)*uR_prim(1)*(wavespeeds(1) - uR_prim(1))) / (uL(0)*(wavespeeds(0) - uL_prim(1)) - uR(0)*(wavespeeds(1) - uR_prim(1)));

  return wavespeeds;
}

double NumericalMethod::deltai_func(double u_i, double u_iPlus1, double u_iMinus1, double w=0.0)
{
  // calculates cell delta value for boundary extrapolated reconstruction
  return 0.5 * (1.0 + w) * (u_i - u_iMinus1) + 0.5 * (1.0 - w) * (u_iPlus1 - u_i);
}

double NumericalMethod::slope_limiter(double u_i, double u_iPlus1, double u_iMinus1)
{
  // calculates slope ratio value
  
  if ((u_iPlus1 - u_i) == 0.0){
    return 0.0;
  }
  return (u_i - u_iMinus1) / (u_iPlus1 - u_i);
}

double NumericalMethod::limiterXi(double r, Limiters limiter)
{

  double Xi;
  
  if (limiter == Limiters::Superbee)
    {
      Xi = superbee(r);
    }
  else if (limiter == Limiters::Van_Leer)
    {
      Xi = van_leer(r);
    }
  else if (limiter == Limiters::Van_Albada)
    {
      Xi = van_albada(r);
    }
  else
    {
      Xi = minbee(r);
    }

  return Xi;
}

double NumericalMethod::reconstruction_uL(double u_i, double u_iPlus1, double u_iMinus1, Limiters limiter)
{
  double r = slope_limiter(u_i, u_iPlus1, u_iMinus1);
  double Xi = limiterXi(r, limiter);
  double deltai = deltai_func(u_i, u_iPlus1, u_iMinus1);

  return u_i - 0.5 * Xi * deltai;
}

double NumericalMethod::reconstruction_uR(double u_i, double u_iPlus1, double u_iMinus1, Limiters limiter)
{
  double r = slope_limiter(u_i, u_iPlus1, u_iMinus1);
  double Xi = limiterXi(r, limiter);
  double deltai = deltai_func(u_i, u_iPlus1, u_iMinus1);

  return u_i + 0.5 * Xi * deltai;
}

double NumericalMethod::minbee(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else if (r > 1.0)
    {
      return std::min(1.0, (2.0 / (1.0 + r)));
    }
  else
    {
      return r;
    }
}

double NumericalMethod::superbee(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else if (r > 1.0 && r <= 0.5)
    {
      return 2.0 * r;
    }
  else if (r < 0.5 && r <= 1.0)
    {
      return 1.0;
    }
  else
    {
      return std::min({r, (2.0 / (1.0 + r)), 2.0});
    }
}

double NumericalMethod::van_leer(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else
    {
      return std::min({((2.0 * r) / (1.0 + r)), (2.0 / (1.0 + r))});
    }
}

double NumericalMethod::van_albada(double r)
{
  if (r <= 0.0)
    {
      return 0.0;
    }
  else
    {
      return std::min((r*(1.0 + r) / (1 + pow(r, 2.0))), (2.0 / (1.0 + r)));
    }
}

double NumericalMethod::lax_friedrich_flux(double u_i, double u_iPlus1, double flux_i, double flux_Plus1,  double dt, double dx)
{
  
  double fhalf;

  fhalf = 0.5 * (dx / dt) * (u_i - u_iPlus1) + 0.5 * (flux_Plus1 + flux_i);
  
  return fhalf;
}

double NumericalMethod::richtmyer_flux(double u_i, double u_iPlus1, double flux_i, double flux_Plus1,  double dt, double dx)
{

  double uhalf;
  double fhalf;

  uhalf = 0.5 * (u_i + u_iPlus1) - 0.5 * (dt / dx) * (flux_Plus1 - flux_i);
  
  return uhalf;
}

var_arr NumericalMethod::FORCE_flux(var_arr uLhalf,  var_arr uRhalf, var_arr uLhalf_prim, var_arr uRhalf_prim, double dt, double dx)
{
  
  var_arr fhalf(nCells+1, uLhalf.cols());
  var_arr uRflux = eos.Euler_flux_fn(uRhalf, uRhalf_prim);
  var_arr uLflux = eos.Euler_flux_fn(uLhalf, uLhalf_prim);
  var_arr uhalf(nCells+1, uLhalf.cols());
  var_arr RI_flux(nCells+1, uLhalf.cols());
  var_arr LF_flux(nCells+1, uLhalf.cols());

  for (int i=0; i<nCells+1; i++)
    {
      for (int j=0; j<uLhalf.cols(); j++)
	{
	  uhalf(i, j) = richtmyer_flux(uRhalf(i, j), uLhalf(i+1, j), uRflux(i, j), uLflux(i+1, j), dt, dx);
	  LF_flux(i, j) = lax_friedrich_flux(uRhalf(i, j), uLhalf(i+1, j), uRflux(i, j), uLflux(i+1, j), dt, dx);
	}
    }

  RI_flux = eos.Euler_flux_fn(uhalf, eos.con_to_prim(uhalf));

  for (int i=0; i<nCells+1; i++)
    {
      for (int j=0; j<uLhalf.cols(); j++)
	{
	  fhalf(i, j) = 0.5 * (LF_flux(i, j) + RI_flux(i, j));
	}
    }

  return fhalf;
}

var_arr NumericalMethod::uHLLC(var_arr u, var_arr u_prim, double S, double S_star)
{
  var_arr uHLLC_K(1, 3);

  uHLLC_K(0) = u(0) * ((S - u_prim(1))/(S - S_star));
  uHLLC_K(1) = uHLLC_K(0) * S_star;
  uHLLC_K(2) = uHLLC_K(0) * (u(2)/u(0) + (S_star - u_prim(1)) * (S_star + u_prim(2) / (u(0)*(S - u_prim(1)))));
  /*
  std::cout<<u<<std::endl;
  std::cout<<u_prim<<std::endl;
  std::cout<<S<<" "<<S_star<<std::endl;
  std::cout<<uHLLC_K<<std::endl;
  std::cout<<std::endl;
  */

  return uHLLC_K;
}

var_arr NumericalMethod::fHLLC(var_arr uL, var_arr uR, var_arr uLHLLC, var_arr uRHLLC, var_arr fL, var_arr fR, double SL, double SR, double S_star)
{
  var_arr fHLLC(uLHLLC.rows(), uLHLLC.cols());
  
  if (SL >= 0)
    {
      fHLLC = fL;
    }
  else if (SL < 0 && S_star >= 0)
    {
      fHLLC = fL + SL*(uLHLLC - uL);
    }
  else if (S_star < 0 && SR >= 0)
    {
      fHLLC = fR + SR*(uRHLLC - uR);
    }
  else
    {
      fHLLC = fR;
    }

  return fHLLC;
}

var_arr NumericalMethod::HLLC_flux(var_arr uL, var_arr uR, var_arr uL_prim, var_arr uR_prim, double dt, double dx)
{
  var_arr fhalf(nCells+1, uL.cols());
  var_arr uRflux = eos.Euler_flux_fn(uR, uR_prim);
  var_arr uLflux = eos.Euler_flux_fn(uL, uL_prim);

  var_arr uLHLLC(nCells+1, uL.cols());
  var_arr uRHLLC(nCells+1, uR.cols());

  for (int i=0; i<nCells+1; i++)
    {
      var_arr wavespeeds = wavespeed(uR.row(i), uL.row(i+1), uR_prim.row(i), uL_prim.row(i+1));
      double SL=wavespeeds(0), SR=wavespeeds(1), S_star=wavespeeds(2);

      uLHLLC.row(i) = uHLLC(uR.row(i), uR_prim.row(i), SL, S_star);
      uRHLLC.row(i) = uHLLC(uL.row(i+1), uL_prim.row(i+1), SR, S_star);

      fhalf.row(i) = fHLLC(uR.row(i), uL.row(i+1), uLHLLC.row(i), uRHLLC.row(i), uRflux.row(i), uLflux.row(i+1), SL, SR, S_star);
    }

  return fhalf;
}

var_arr NumericalMethod::uL_half_update(var_arr uL, var_arr uR, double dt, double dx)
{
  
  // data
  var_arr uL_flux = eos.Euler_flux_fn(uL, eos.con_to_prim(uL));
  var_arr uR_flux = eos.Euler_flux_fn(uR, eos.con_to_prim(uR));

  return uL - 0.5 * (dt / dx) * (uR_flux - uL_flux);
  
}

var_arr NumericalMethod::uR_half_update(var_arr uL, var_arr uR, double dt, double dx)
{
  
  // data
  var_arr uL_flux = eos.Euler_flux_fn(uL, eos.con_to_prim(uL));
  var_arr uR_flux = eos.Euler_flux_fn(uR, eos.con_to_prim(uR));

  return uR - 0.5 * (dt / dx) * (uR_flux - uL_flux);
  
}
