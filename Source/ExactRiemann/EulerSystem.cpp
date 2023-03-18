#include "EulerSystem.h"

EulerSystem::EulerSystem()
  :SettingsData(){};

void EulerSystem::initial_conds()
{
  // Set initial conditions for simulation
  InitialCondTests Tests;
  for (int i=0; i<u.rows(); i++)
    {
      double x = x0 + (double(i) + 0.5) * dx;

      var_array current = Tests.Test1D_3(x);

      u_prim(i,0) = current[0];
      u_prim(i,1) = current[1];
      u_prim(i,2) = current[2];
    }
}

void EulerSystem::resize_matrix()
{
  // resize matrices to have ghost cells
  u.resize(nCells+4, 4);
  uPlus1.resize(nCells+4, 3);
  u_prim.resize(nCells+4, 3);
}

void EulerSystem::outputFile(std::string outputName)
{
  // output simulation results
  std::ofstream output(outputName);

  for (int i=1; i<u.rows()-3; i++)
    {
      double x = x0 + (double(i-1)+0.5) * dx;
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

  // Calculate sound speed
  Eigen::ArrayXXf uL_prim = u_prim.row(0);
  Eigen::ArrayXXf uR_prim = u_prim.row(nCells+1);
  Eigen::ArrayXXf cs = NumM.SoundSpeed(uL_prim, uR_prim);

  // approximate pressure guess
  double pressure = NumM.PressureApprox(uL_prim, uR_prim);

  // calculate intermediate pressure in star states using Newton-Raphson method
  std::array<double,2> star_vals = NumM.NewtonRaphson(pressure, uL_prim, uR_prim, cs);
  double p_star = star_vals[0];
  double v_star = star_vals[1];
  double x_val=0;

  Eigen::ArrayXXf u_current;

  // loop through cells and sample the value in each point based on wave patterns
  for (int i=0; i<nCells+2; i++)
    {
      x_val += dx;
      double pos = (x_val - 0.5) / t;
      u_current = NumM.SamplePattern(p_star, v_star, pos, uL_prim, uR_prim, cs);
      u(i,0) = u_current(0);
      u(i,1) = u_current(1);
      u(i,2) = u_current(2);
      u(i,3) = u_current(2) / ((gamma-1.0)*u_current(0));
    }
}

void EulerSystem::run()
{

  // call simulation functions
  resize_matrix();
    
  initial_conds();

  Solver1D();

}
