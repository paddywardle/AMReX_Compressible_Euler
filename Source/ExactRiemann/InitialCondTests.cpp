#include "InitialCondTests.h"

var_array InitialCondTests::Test1D_1(double x)
{

  var_array initial_conds;
  
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

var_array InitialCondTests::Test1D_2(double x)
{

  var_array initial_conds;
  
  if (x < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = -2.0;
      // Pressure
      initial_conds[2] = 0.4;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 2.0;
      // Pressure
      initial_conds[2] = 0.4;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_3(double x)
{

  var_array initial_conds;
  
  if (x < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1000.0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.01;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_4(double x)
{

  var_array initial_conds;
  
  if (x < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.01;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 100.0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_5(double x)
{

  var_array initial_conds;
  
  if (x < 0.5)
    {
      // Density
      initial_conds[0] = 5.99924;
      // Velocity
      initial_conds[1] = 19.5975;
      // Pressure
      initial_conds[2] = 460.894;
    }
  else
    {
      // Density
      initial_conds[0] = 5.99242;
      // Velocity
      initial_conds[1] = -6.19633;
      // Pressure
      initial_conds[2] = 46.0950;
    }

  return initial_conds;
}


