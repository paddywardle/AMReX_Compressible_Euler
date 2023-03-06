#include "NumericalMethod.h"

NumericalMethod::NumericalMethod(double gamma, int nCells)
  :gamma(gamma), nCells(nCells), eos(gamma){};

Eigen::ArrayXXf NumericalMethod::PressureApprox(Eigen::ArrayXXf p_guess, Eigen::ArrayXXf u_prim)
{
  // ALTER VARIABLE NAMES
  int Quser=2.0;
  int rows=u_prim.rows();
  std::vector<double, nCells+1> p_g;

  for (int i=0; i<nCells+1; i++)
    {
      Eigen::ArrayXXf uL_prim = u_prim.row(i);
      Eigen::ArrayXXf uR_prim = u_prim.row(i+1);

      double cs_left = sqrt(gamma*uL_prim(2)/uL_prim(0));
      double cs_right = sqrt(gamma*uR_prim(2)/uR_prim(0));

      double rho_int = 0.25 * (uL_prim(0) + uR_prim(0)) * (cs_left * cs_right);
      double ppv = 0.5 * (uL_prim(2) + uR_prim(2)) + 0.5 * (uL_prim(1) - uR_prim(1))*rho_int;
      ppv = std::max(0.0, ppv);
      double minp = std::min(uL_prim(2), uR_prim(2));
      double maxp = std::max(uL_prim(2), uR_prim(2));
      double qmax = maxp/minp;

      if ((qmax <= Quser) && (minp <= ppv) && (maxp >= ppv))
	{
	  p_g[i] = ppv;
	}
      else
	{
	  if (ppv < minp)
	    {
	      // two-rarefaction riemann solver
	      double PQ = pow((uL_prim(2) / uR_prim(2)), (gamma - 1.0) / (2.0 * gamma));
	      double UM = (PQ * uL_prim(1) / cs_left + uR_prim(1)/cs_right + (2.0/(gamma-1.0))*(PQ-1.0))/(PQ/cs_left + 1.0/cs_right);
	      double ptl = 1.0 + ((gamma-1.0)/2.0)*(uL_prim(1)-UM)/cs_left;
	      double ptr = 1.0 + ((gamma-1.0)/2.0)*(UM-uR_prim(1))/cs_right;
	      p_g[i] = 0.5*(pow(uL_prim(2)*ptl, (2.0*gamma)/(gamma-1.0)) + pow(uR_prim(2)*ptr, (2.0*gamma)/(gamma-1.0)));
	    }
	  else
	    {
	      // two-shock riemann solver
	      double G5 = 2.0 / (gamma+1.0);
	      double G6 = (gamma - 1.0) / (gamma + 1.0);
	      double GEL = sqrt((G5/uL_prim(0)) / (G6*uL_prim(2) + ppv));
	      double GER = sqrt((G5/uR_prim(0)) / (G6*uR_prim(2) + ppv));
	      p_g[i] = (GEL * uL_prim(2) + GER * uR_prim(2) - (uR_prim(1) - uL_prim(1))) / (GEL + GER);
	    }
	}
    }
}

Eigen::ArrayXXf NumericalMethod::SoundSpeed(Eigen::ArrayXXf u_prim)
{

  Eigen::ArrayXXf cs(nCells+1, 2);

  for (int i=0; i<nCells+1; i++)
    {
      cs(i, 0) = sqrt(gamma*u_prim(i,2)/u_prim(i,0));
      cs(i, 1) = sqrt(gamma*u_prim(i+1,2)/u_prim(i+1,0));
    }

  return cs;
  
}

Eigen::ArrayXXf NumericalMethod::InitialPressure(Eigen::ArrayXXf u_prim)
{
  Eigen::ArrayXXf p_init(nCells+1);

  for (int i=0; i<p_init.rows(); i++)
    {
      p_init(i,0) = 0.0;
    }

  return p_init;
}

std::array<double,2> NumericalMethod::PressureFlux(double p_prev, Eigen::ArrayXXf uK, double csK)
{
  double AK = 2.0/((gamma+1)*uK(0));
  double BK = (gamma - 1.0)/(gamma + 1.0);
  std::array<double, 2> pK_flux;
  if (p_prev > uK(2))
    {
      pK_flux[0] = (p_prev - uK(2)) * sqrt(AK / (p_prev + BK));
      pK_flux[1] = (1.0 - 0.5*(p_prev - uK(2))/(BK + p_prev)) * sqrt(AK / (p_prev + BK));
    }
  else
    {
      pK_flux[0] = ((2 * csK) / (gamma - 1.0)) * (pow((p_prev/uK(2)), ((gamma-1.0)/(2*gamma))) - 1.0);
      pK_flux[1] = (1.0/(uK(0) * csK)) * pow(p_prev/uK(2),-(gamma+1.0)/(2*gamma)); 
    }

  return pK_flux;
}

Eigen::ArrayXXf NumericalMethod::NewtonRaphson(Eigen::ArrayXXf pressure, Eigen::ArrayXXf u_prim, Eigen::ArrayXXf cs)
{
  for (int i=0; i<nCells+1; i++)
    {
      Eigen::ArrayXXf uL = u_prim.row(i);
      Eigen::ArrayXXf uR = u_prim.row(i+1);
      double p_prev = pressure(i,1);
      double p_new = pressure(i,1);
      double csL = cs(i, 0);
      double csR = cs(i, 1);
      do
	{
	  std::array<double,2> pL_flux = PressureFlux(p_prev, uL, csL);
	  std::array<double,2> pR_flux = PressureFlux(p_prev, uR, csR);
	  p_new = p_prev - ((pL_flux[0] + pR_flux[0] + (uL(1) - uR(1))) / (pL_flux[1] + pR_flux[1]));
	  p_prev = p_new;
	} while ((2.0 * fabs((p_new - p_prev)/(p_new + p_prev))) > 1e-9);

      pressure(i,1) = p_prev;
    }

  return pressure;
  
}
