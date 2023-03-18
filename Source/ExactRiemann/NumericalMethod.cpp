#include "NumericalMethod.h"

double gam=1.4;

NumericalMethod::NumericalMethod(double gamma, int nCells)
  :gamma(gamma), nCells(nCells){};

Eigen::ArrayXXf NumericalMethod::leftRar(double p_star, double v_star, double pos, Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim, Eigen::ArrayXXf cs)
{

  // calculate density, velocity and pressure for left rarefaction
  double leftRar = uL_prim(1) - cs(0);
  double rho;
  double velocity;
  double pressure;

  // if in left state
  if (pos <= leftRar)
    {
      rho = uL_prim(0);
      velocity = uL_prim(1);
      pressure = uL_prim(2);
    }
  else
    {
      double c_starL = cs(0) * pow((p_star/uL_prim(2)), (gam - 1.0)/(2.0*gam));
      double st_left = v_star - c_starL;

      // if in left intermediate state
      if (pos > st_left)
	{
	  rho = uL_prim(0) * pow((p_star/uL_prim(2)), (1.0/gam));
	  velocity = v_star;
	  pressure = p_star;
	}
      else
	{
	  double c = (2.0/(gam + 1.0))*(cs(0) + ((gam-1.0)/2.0)*(uL_prim(1) - pos));
	  rho = uL_prim(0) * pow((c/cs(0)), (2.0/(gam-1.0))); 
	  velocity = (2.0/(gam + 1.0))*(cs(0) + ((gam-1.0)/2.0)*uL_prim(1) + pos);
	  pressure = uL_prim(2) * pow((c/cs(0)), (2.0 * gam)/(gam - 1.0));
	}
    }
  
  Eigen::ArrayXXf u(1,3);

  u<<rho,velocity,pressure;

  return u;
  
}

Eigen::ArrayXXf NumericalMethod::leftShock(double p_star, double v_star, double pos, Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim, Eigen::ArrayXXf cs)
{

  // calculate density, velocity and pressure in left shock
  double p_starL = p_star/uL_prim(2);
  double SL = uL_prim(1) - cs(0) * sqrt(((gam+1.0)/(2.0*gam)) * p_starL + ((gam-1.0)/(2.0*gam)));
  double rho;
  double velocity;
  double pressure;

  // if in left state
  if (pos <= SL)
    {
      rho = uL_prim(0);
      velocity = uL_prim(1);
      pressure = uL_prim(2);
    }
  // if in left intermediate state
  else
    {
      rho = uL_prim(0) * (p_starL + ((gam-1.0)/(gam+1.0))) / (p_starL*((gam-1.0)/(gam+1.0)) + 1.0);
      velocity = v_star;
      pressure = p_star;
    }
  
  Eigen::ArrayXXf u(1,3);

  u<<rho,velocity,pressure;

  return u;
  
}

Eigen::ArrayXXf NumericalMethod::rightRar(double p_star, double v_star, double pos, Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim, Eigen::ArrayXXf cs)
{
  // calculate density, velocity and pressure for right rarefaction
  double rightRar = uR_prim(1) + cs(1);
  double rho;
  double velocity;
  double pressure;

  // if in right state
  if (pos >= rightRar)
    {
      rho = uR_prim(0);
      velocity = uR_prim(1);
      pressure = uR_prim(2);
    }
  else
    {
      double c_starR = cs(1) * pow((p_star/uR_prim(2)), (gam - 1.0)/(2.0*gam));
      double st_right = v_star + c_starR;

      // if in right intermediate state
      if (pos <= st_right)
	{
	  rho = uR_prim(0) * pow((p_star/uR_prim(2)), (1.0/gam));
	  velocity = v_star;
	  pressure = p_star;
	}
      else
	{
	  double c = (2.0/(gam + 1.0))*(cs(1) - ((gam-1.0)/2.0)*(uR_prim(1) - pos));
	  rho = uR_prim(0) * pow((c/cs(1)), (2.0/(gam-1.0)));
	  velocity = (2.0/(gam + 1.0))*(-cs(1) + ((gam-1.0)/2.0)*uR_prim(1) + pos);
	  pressure = uR_prim(2) * pow((c/cs(1)), (2.0 * gam)/(gam - 1.0));
	}
    }
  
  Eigen::ArrayXXf u(1,3);

  u<<rho,velocity,pressure;

  return u;
  
}

Eigen::ArrayXXf NumericalMethod::rightShock(double p_star, double v_star, double pos, Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim, Eigen::ArrayXXf cs)
{

  // calculate density, velocity and pressure for right shock
  double p_starR = p_star/uR_prim(2);
  double SR = uR_prim(1) + cs(1) * sqrt(((gam+1.0)/(2.0*gam)) * p_starR + ((gam-1.0)/(2.0*gam)));
  double rho;
  double velocity;
  double pressure;

  // if in right state
  if (pos >= SR)
    {
      rho = uR_prim(0);
      velocity = uR_prim(1);
      pressure = uR_prim(2);
    }
  // if in right intermediate state
  else
    {
       rho = uR_prim(0) * (p_starR + ((gam-1.0)/(gam+1.0))) / (p_starR*((gam-1.0)/(gam+1.0)) + 1.0);
       velocity = v_star;
       pressure = p_star;
    }
  
  Eigen::ArrayXXf u(1,3);

  u<<rho,velocity,pressure;

  return u;
  
}

Eigen::ArrayXXf NumericalMethod::SamplePattern(double p_star, double v_star, double pos, Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim, Eigen::ArrayXXf cs)
{
  Eigen::ArrayXXf u;

  // left of contact discontinuity
  if (pos <= v_star)
    {
      // left rarefaction
      if (p_star <= uL_prim(2))
	{
	  u = leftRar(p_star, v_star, pos, uL_prim, uR_prim, cs);
	}
      //right rarefaction
      else
	{
	  u = leftShock(p_star, v_star, pos, uL_prim, uR_prim, cs);
	}
    }

  // right of contact discontinuity
  else
    {
      // right rarefaction
      if (p_star <= uR_prim(2))
	{
	  u = rightRar(p_star, v_star, pos, uL_prim, uR_prim, cs);
	}
      // right shock
      else
	{
	  u = rightShock(p_star, v_star, pos, uL_prim, uR_prim, cs);
	} 
    }

  return u;
  
}

double NumericalMethod::PressureApprox(Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim)
{
  int Quser=2.0;
  double p_g;

  // sound speed
  double cs_left = sqrt(gam*uL_prim(2)/uL_prim(0));
  double cs_right = sqrt(gam*uR_prim(2)/uR_prim(0));

  double rho_int = 0.25 * (uL_prim(0) + uR_prim(0)) * (cs_left + cs_right);

  // pressure guess
  double ppv = 0.5 * (uL_prim(2) + uR_prim(2)) + 0.5 * (uL_prim(1) - uR_prim(1))*rho_int;
  ppv = std::max(0.0, ppv);
  double minp = std::min(uL_prim(2), uR_prim(2));
  double maxp = std::max(uL_prim(2), uR_prim(2));
  double qmax = maxp/minp;

  if ((qmax <= Quser) && (minp <= ppv &&  ppv <= maxp))
    {
      p_g = ppv;
    }
  else
    {
      if (ppv < minp)
	{
	  // two-rarefaction solution
	  double PQ = pow((uL_prim(2) / uR_prim(2)), (gam - 1.0) / (2.0 * gam));
	  double UM = (PQ * uL_prim(1) / cs_left + uR_prim(1)/cs_right + (2.0/(gam-1.0))*(PQ-1.0))/(PQ/cs_left + 1.0/cs_right);
	  double ptleft = 1.0 + ((gam-1.0)/2.0)*(uL_prim(1)-UM)/cs_left;
	  double ptright = 1.0 + ((gam-1.0)/2.0)*(UM-uR_prim(1))/cs_right;
	  p_g = 0.5*(pow(uL_prim(2)*ptleft, (2.0*gam)/(gam-1.0)) + pow(uR_prim(2)*ptright, (2.0*gam)/(gam-1.0)));
	}
      else
	{
	  // two-shock solution
	  double Gleft = sqrt(((2.0 / (gam+1.0))/uL_prim(0)) / (((gam - 1.0) / (gam + 1.0))*uL_prim(2) + ppv));
	  double Gright = sqrt(((2.0 / (gam+1.0))/uR_prim(0)) / (((gam - 1.0) / (gam + 1.0))*uR_prim(2) + ppv));
	  p_g = (Gleft * uL_prim(2) + Gright * uR_prim(2) - (uR_prim(1) - uL_prim(1))) / (Gleft + Gright);
	}
    }

  return p_g;
}

Eigen::ArrayXXf NumericalMethod::SoundSpeed(Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim)
{

  // calculate left and right sound speeds
  Eigen::ArrayXXf cs(1, 2);

  cs(0, 0) = sqrt(gam*uL_prim(2)/uL_prim(0));
  cs(0, 1) = sqrt(gam*uR_prim(2)/uR_prim(0));

  return cs;
  
}

std::array<double,2> NumericalMethod::PressureFlux(double p_prev, Eigen::ArrayXXf uK, double csK)
{
  double AK = 2.0/((gam+1)*uK(0));
  double BK = ((gam - 1.0)/(gam + 1.0)) * uK(2);
  std::array<double, 2> pK_flux;
  if (p_prev > uK(2))
    {
      pK_flux[0] = (p_prev - uK(2)) * sqrt(AK / (p_prev + BK));
      pK_flux[1] = (1.0 - 0.5*(p_prev - uK(2))/(BK + p_prev)) * sqrt(AK / (p_prev + BK));
    }
  else
    {
      pK_flux[0] = ((2.0 * csK) / (gam - 1.0)) * (pow((p_prev/uK(2)), ((gam-1.0)/(2.0*gam))) - 1.0);
      pK_flux[1] = (1.0/(uK(0) * csK)) * pow((p_prev/uK(2)),-(gam+1.0)/(2.0*gam)); 
    }

  return pK_flux;
}

std::array<double,2> NumericalMethod::NewtonRaphson(double pressure, Eigen::ArrayXXf uL_prim, Eigen::ArrayXXf uR_prim, Eigen::ArrayXXf cs)
{
  
  double p_prev = pressure;
  double p_new = pressure;
  double pLflux;
  double pRflux;
  double csL = cs(0);
  double csR = cs(1);
  double current_tol;

  // loop while not converged
  do
    {
      // calculate left and right pressure flux
      std::array<double,2> pL_flux = PressureFlux(p_prev, uL_prim, csL);
      std::array<double,2> pR_flux = PressureFlux(p_prev, uR_prim, csR);

      // Newton Raphson update
      p_new = p_prev - ((pL_flux[0] + pR_flux[0] + (uR_prim(1) - uL_prim(1))) / (pL_flux[1] + pR_flux[1]));

      // calculate new residual
      current_tol = (2.0 * fabs((p_new - p_prev)/(p_new + p_prev)));
      p_prev = p_new;
      pLflux = pL_flux[0];
      pRflux = pR_flux[0];
    } while (current_tol > 1e-9);

  // calcualte intermediate velocity in star state based on intermediate pressure
  double v_star = 0.5 * (uL_prim(1) + uR_prim(1) + pRflux - pLflux);

  std::array<double,2> star_vals;
  star_vals[0] = p_new;
  star_vals[1] = v_star;
  
  return star_vals;
  
}
