#include "EulerSystem.h"

EulerSystem::EulerSystem()
  :SettingsData(){};

void EulerSystem::initial_conds()
{
  for (int i=0; i<u.rows(); i++)
    {
      double x = x0 + i * dx;
      
      if (x <= 0.5)
	{
	  // density
	  u_prim(i, 0) = 1.0;
	  // velocity
	  u_prim(i, 1) = 0.0;
	  // pressure
	  u_prim(i, 2) = 1.0;
	}
      else
	{
	  // density
	  u_prim(i, 0) = 0.125;
	  // velocity
	  u_prim(i, 1) = 0.0;
	  // pressure
	  u_prim(i, 2) = 0.1;
	}
    }
}

void EulerSystem::resize_matrix()
{
  u.resize(nCells+2, 3);
  u_prim.resize(nCells+2, 3);
}

void EulerSystem::outputFile(std::string outputName)
{
  std::ofstream output(outputName);

  for (int i=1; i<u_prim.rows()-2; i++)
    {
      double x = x0 + (i-1) * dx;
      output<<x<<" ";
      for (int j=0; j<u_prim.cols(); j++)
	{
	  output<<u_prim(i, j)<<" ";
	}
      output<<std::endl;
    }
}

void EulerSystem::Solver1D()
{
  // sets initial conditions into instance variables
  initial_conds();

  // set current time to simulation start time
  double t = tStart;

  // Instance of helper classes
  EulerEOS eos(gamma);
  NumericalMethod NumM(gamma, nCells);

  Eigen::ArrayXXf cs = NumM.SoundSpeed(u_prim);

  Eigen::ArrayXXf p_init = NumM.InitialPressure(u_prim);

  Eigen::ArrayXXf pressure = NumM.NewtonRaphson(p_init, u_prim, cs);

}

void EulerSystem::run()
{

  EulerEOS eos(gamma);
  
  resize_matrix();
    
  initial_conds();
  
  u = eos.prim_to_con(u_prim);
  
  uPlus1 = eos.prim_to_con(u_prim);

  Solver1D();

  u_prim = eos.con_to_prim(u);

}
