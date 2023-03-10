#include "EulerSystem.h"

EulerSystem::EulerSystem()
  :SettingsData(){};

void EulerSystem::initial_conds()
{
  InitialCondTests Tests;
  for (int i=0; i<u.rows(); i++)
    {
      double x = x0 + i * dx;

      var_array current = Tests.Test1D_1(x);

      u_prim(i,0) = current[0];
      u_prim(i,1) = current[1];
      u_prim(i,2) = current[2];
    }
}

void EulerSystem::resize_matrix()
{
  u.resize(nCells+4, 3);
  uPlus1.resize(nCells+4, 3);
  u_prim.resize(nCells+4, 3);
}

void EulerSystem::outputFile(std::string outputName)
{
  std::ofstream output(outputName);

  for (int i=1; i<u.rows()-2; i++)
    {
      double x = x0 + (i-1) * dx;
      output<<x<<" ";
      for (int j=0; j<u.cols(); j++)
	{
	  output<<u(i, j)<<" ";
	}
      output<<std::endl;
    }
}

void EulerSystem::Solver1D()
{

  // set current time to simulation start time
  double t = tStop;

  // Instance of helper classes
  NumericalMethod NumM(gamma, nCells);

  Eigen::ArrayXXf uL_prim = u_prim.row(0);
  Eigen::ArrayXXf uR_prim = u_prim.row(nCells);
  Eigen::ArrayXXf cs = NumM.SoundSpeed(uL_prim, uR_prim);
  double pressure = NumM.PressureApprox(uL_prim, uR_prim);
  std::array<double,2> star_vals = NumM.NewtonRaphson(pressure, uL_prim, uR_prim, cs);
  double p_star = star_vals[0];
  double v_star = star_vals[1];
  double x_val=0;

  for (int i=0; i<nCells+1; i++)
    {
      x_val += dx;
      double pos = (x_val - 0.5) / t;
      u.row(i) = NumM.SamplePattern(p_star, v_star, pos, u_prim.row(i), u_prim.row(i+1), cs);
    }
}

void EulerSystem::run()
{
  
  resize_matrix();
    
  initial_conds();

  Solver1D();

}
