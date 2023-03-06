#include "common.h"

/*
Eigen::ArrayXXf prim_to_con(Eigen::ArrayXXf u_p)
{
  Eigen::ArrayXXf u_c(u_p.rows(), u_p.cols());

  u_c.col(0) = u_p.col(0);
  u_c.col(1) = u_p.col(0) * u_p.col(1);
  u_c.col(2) = u_p.col(2) / (gamma - 1.0) + 0.5 * u_p.col(0) * pow(u_p.col(1), 2);

  return u_c;
}

Eigen::ArrayXXf con_to_prim(Eigen::ArrayXXf u_c)
{
  Eigen::ArrayXXf u_p(u_c.rows(), u_c.cols());

  u_p.col(0) = u_c.col(0);
  u_p.col(1) = u_c.col(1) / u_c.col(0);
  u_p.col(2) = (u_c.col(2) - 0.5 * u_c.col(0) * pow(u_p.col(1), 2.0)) * (gamma - 1.0);
  
  return u_p;
}
*/
Eigen::ArrayXXf PressureApprox(Eigen::ArrayXXf p_guess, Eigen::ArrayXXf u_prim)
{

  int const nCells = 100;
  double x0=0;
  double x1=1;
  double dx = (x1-x0)/nCells;
  double gamma=1.4;
  // ALTER VARIABLE NAMES
  int Quser=2.0;
  int const rows=u_prim.rows();
  std::vector<double, 101> p_g;

  for (int i=0; i<nCells+1; i++)
    {
      double rhoL = u_prim(i,0);
      double uL = u_prim(i,1);
      double pL = u_prim(i,2);
      double rhoR = u_prim(i+1,0);
      double uR = u_prim(i+1,1);
      double pR = u_prim(i+1,2);

      double cs_left = sqrt(gamma*pL/rhoL);
      double cs_right = sqrt(gamma*pR/rhoR);

      double rho_int = 0.25 * (rhoL + rhoR) * (cs_left * cs_right);
      double ppv = 0.5 * (pL + pR) + 0.5 * (uL - uR)*rho_int;
      ppv = std::max(0.0, ppv);
      double minp = std::min(pL, pR);
      double maxp = std::max(pL, pR);
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
	      double PQ = pow((pL / pR), (gamma - 1.0) / (2.0 * gamma));
	      double UM = (PQ * uL / cs_left + uR/cs_right + (2.0/(gamma-1.0))*(PQ-1.0))/(PQ/cs_left + 1.0/cs_right);
	      double ptl = 1.0 + ((gamma-1.0)/2.0)*(uL-UM)/cs_left;
	      double ptr = 1.0 + ((gamma-1.0)/2.0)*(UM-uR)/cs_right;
	      p_g[i] = 0.5*(pow(pL*ptl, (2.0*gamma)/(gamma-1.0)) + pow(pR*ptr, (2.0*gamma)/(gamma-1.0)));
	    }
	  else
	    {
	      // two-shock riemann solver
	      double G5 = 2.0 / (gamma+1.0);
	      double G6 = (gamma - 1.0) / (gamma + 1.0);
	      double GEL = sqrt((G5/rhoL) / (G6*pL + ppv));
	      double GER = sqrt((G5/rhoR) / (G6*pR + ppv));
	      p_g[i] = (GEL * pL + GER * pR - (uR - uL)) / (GEL + GER);
	    }
	}
    }
}
/*

Eigen::ArrayXXf SoundSpeed(Eigen::ArrayXXf u_prim)
{

  Eigen::ArrayXXf cs(nCells+1, 2);

  for (int i=0; i<nCells+1; i++)
    {
      cs(i, 0) = sqrt(gamma*u_prim(i,2)/u_prim(i,0));
      cs(i, 1) = sqrt(gamma*u_prim(i+1,2)/u_prim(i+1,0));
    }

  return cs;
  
}

Eigen::ArrayXXf InitialPressure(Eigen::ArrayXXf u_prim)
{
  Eigen::ArrayXXf p_init(nCells+1);

  for (int i=0; i<p_init.rows(); i++)
    {
      p_init(i,0) = 0.0;
    }

  return p_init;
}

std::array<double,2> PressureFlux(double p_prev, Eigen::ArrayXXf uK, double csK)
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

Eigen::ArrayXXf NewtonRaphson(Eigen::ArrayXXf pressure, Eigen::ArrayXXf u_prim, Eigen::ArrayXXf cs)
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

Eigen::ArrayXXf InitialConds(double x)
{

  Eigen:: initial_conds;
  
  if (x < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1.0;
    }
  else
    {
      // Density
      initial_conds[0] = 0.125;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.1;
    }

  return initial_conds;
}




void main()
{

  Eigen::ArrayXXf u_prim(nCells+2,3);

  for (int i=0; i<nCells; i++)
    {
      double x= x0 + i * dx;

      if (x < 0.5)
	{
	  // Density
	  u_prim(0) = 1.0;
	  // Velocity
	  u_prim(1) = 0.0;
	  // Pressure
	  u_prim(2) = 1.0;
	}
      else
	{
	  // Density
	  u_prim(0) = 0.125;
	  // Velocity
	  u_prim(1) = 0.0;
	  // Pressure
	  u_prim(2) = 0.1;
	}
    }

  //Eigen::ArrayXXf u_con = prim_to_con(u_prim);

  Eigen::ArrayXXf cs = NumM.SoundSpeed(u_prim);

  Eigen::ArrayXXf p_init = NumM.InitialPressure(u_prim);

  Eigen::ArrayXXf pressure = NumM.NewtonRaphson(p_init, u_prim, cs);

}
