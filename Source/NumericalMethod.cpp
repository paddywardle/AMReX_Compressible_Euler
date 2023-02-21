#include "NumericalMethod.h"

NumericalMethod::NumericalMethod(double gamma)
  :gamma(gamma), eos(gamma){};

var_array NumericalMethod::wavespeed(var_array uL, var_array uR, var_array uL_prim, var_array uR_prim)
{
  var_array wavespeeds;

  double csL = sqrt((gamma*uL_prim[2])/uL[0]);
  double csR = sqrt((gamma*uR_prim[2])/uR[0]);

  wavespeeds[0] = uL_prim[1] - csL;
  wavespeeds[1] = uR_prim[1] + csR;
  wavespeeds[2] = (uR_prim[2] - uL_prim[2] + uL[0]*uL_prim[1]*(wavespeeds[0] - uL_prim[1]) - uR[0]*uR_prim[1]*(wavespeeds[1] - uR_prim[1])) / (uL[0]*(wavespeeds[0] - uL_prim[1]) - uR[0]*(wavespeeds[1] - uR_prim[1]));

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

var_array NumericalMethod::uHLLC(var_array u, var_array u_prim, double S, double S_star)
{
  var_array uHLLC_K;

  uHLLC_K[0] = u[0] * ((S - u_prim[1])/(S - S_star));
  uHLLC_K[1] = uHLLC_K[0] * S_star;
  uHLLC_K[2] = uHLLC_K[0] * (u[2]/u[0] + (S_star - u_prim[1]) * (S_star + u_prim[2] / (u[0]*(S - u_prim[1]))));

  return uHLLC_K;
}

var_array NumericalMethod::fHLLC(var_array uL, var_array uR, var_array uLHLLC, var_array uRHLLC, var_array fL, var_array fR, double SL, double SR, double S_star)
{
  var_array fHLLC;

  for (int i=0; i<fHLLC.size(); i++)
    {
  
      if (SL >= 0)
	{
	  fHLLC[i] = fL[i];
	}
      else if (SL < 0 && S_star >= 0)
	{
	  fHLLC[i] = fL[i] + SL*(uLHLLC[i] - uL[i]);
	}
      else if (S_star < 0 && SR >= 0)
	{
	  fHLLC[i] = fR[i] + SR*(uRHLLC[i] - uR[i]);
	}
      else
	{
	  fHLLC[i] = fR[i];
	}
    }

  return fHLLC;
}

var_array NumericalMethod::HLLC_flux(var_array u_i, var_array u_iMinus1, var_array u_iPlus1, var_array u_iPlus2, double dx, double dt)
{

  var_array uL, uR, uLPlus1, uRPlus1;

  for (int i=0; i<uL.size(); i++)
    {
      uL[i] = reconstruction_uL(u_i[i], u_iPlus1[i], u_iMinus1[i], Limiters::Minbee);
      uR[i] = reconstruction_uR(u_i[i], u_iPlus1[i], u_iMinus1[i], Limiters::Minbee);
      uLPlus1[i] = reconstruction_uL(u_iPlus1[i], u_iPlus2[i], u_i[i], Limiters::Minbee);
      uRPlus1[i] = reconstruction_uR(u_iPlus1[i], u_iPlus2[i], u_i[i], Limiters::Minbee);	
    }

  // LOOK at what all of these look like!

  var_array uLhalf = uL_half_update(uL, uR, dt, dx);
  var_array uRhalf = uR_half_update(uL, uR, dt, dx);

  var_array uLhalfPlus1 = uL_half_update(uLPlus1, uRPlus1, dt, dx);
  var_array uRhalfPlus1 = uR_half_update(uLPlus1, uRPlus1, dt, dx);
  
  var_array uLhalf_prim = eos.con_to_prim(uLhalf);
  var_array uRhalf_prim = eos.con_to_prim(uRhalf);

  var_array uLhalfPlus1_prim = eos.con_to_prim(uLhalfPlus1);
  var_array uRhalfPlus1_prim = eos.con_to_prim(uRhalfPlus1);
  
  var_array uRflux = eos.Euler_flux_fn(uRhalf, uRhalf_prim);
  var_array uLflux = eos.Euler_flux_fn(uLhalf, uLhalf_prim);

  var_array uRPlus1flux = eos.Euler_flux_fn(uRhalfPlus1, uRhalfPlus1_prim);
  var_array uLPlus1flux = eos.Euler_flux_fn(uLhalfPlus1, uLhalfPlus1_prim);

  var_array wavespeeds = wavespeed(uRhalf, uLhalfPlus1, uRhalf_prim, uLhalfPlus1_prim); // watch out <- am I doing the correct variables here - check what is being fed from AMReX

  double SL=wavespeeds[0], SR=wavespeeds[1], S_star=wavespeeds[2];
  
  var_array uLHLLC = uHLLC(uRhalf, uRhalf_prim, SL, S_star);
  var_array uRHLLC = uHLLC(uLhalfPlus1, uLhalfPlus1_prim, SR, S_star);

  var_array fhalf = fHLLC(uRhalf, uLhalfPlus1, uLHLLC, uRHLLC, uRflux, uLPlus1flux, SL, SR, S_star);

  return fhalf;
}

var_array NumericalMethod::uL_half_update(var_array uL, var_array uR, double dt, double dx)
{
  // data
  var_array uL_flux = eos.Euler_flux_fn(uL, eos.con_to_prim(uL));
  var_array uR_flux = eos.Euler_flux_fn(uR, eos.con_to_prim(uR));
  var_array uLhalf;

  for (int i=0; i<uLhalf.size(); i++)
    {
      uLhalf[i] = uL[i] - 0.5 * (dt / dx) * (uR_flux[i] - uL_flux[i]);
    }

  return uLhalf;
  
}

var_array NumericalMethod::uR_half_update(var_array uL, var_array uR, double dt, double dx)
{
  
  // data
  var_array uL_flux = eos.Euler_flux_fn(uL, eos.con_to_prim(uL));
  var_array uR_flux = eos.Euler_flux_fn(uR, eos.con_to_prim(uR));
  var_array uRhalf;

  for (int i=0; i<uRhalf.size(); i++)
    {
      uRhalf[i] = uR[i] - 0.5 * (dt / dx) * (uR_flux[i] - uL_flux[i]);
    }

  return uRhalf;
  
}
