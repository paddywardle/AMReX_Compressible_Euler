#include "InitialCondTests.h"

var_array InitialCondTests::Test1D_1(double x)
{

  var_array initial_conds;
  
  if (x < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity x
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1.0;
      // Velocity y
      initial_conds[3] = 0.0;
    }
  else
    {
      // Density
      initial_conds[0] = 0.125;
      // Velocity x
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.1;
      // Velocity y
      initial_conds[3] = 0.0;
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
      // Velocity y
      initial_conds[3] = 0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 2.0;
      // Pressure
      initial_conds[2] = 0.4;
      // Velocity y
      initial_conds[3] = 0;
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
      // Velocity y
      initial_conds[3] = 0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.01;
      // Velocity y
      initial_conds[3] = 0;
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
      // Velocity y
      initial_conds[3] = 0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 100.0;
      // Velocity y
      initial_conds[3] = 0;
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
      // Velocity y
      initial_conds[3] = 0;
    }
  else
    {
      // Density
      initial_conds[0] = 5.99242;
      // Velocity
      initial_conds[1] = -6.19633;
      // Pressure
      initial_conds[2] = 46.0950;
      // Velocity y
      initial_conds[3] = 0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_1Y(double y)
{

  var_array initial_conds;
  
  if (y < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity x
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1.0;
      // Velocity y
      initial_conds[3] = 0.0;
    }
  else
    {
      // Density
      initial_conds[0] = 0.125;
      // Velocity x
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.1;
      // Velocity y
      initial_conds[3] = 0.0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_2Y(double y)
{

  var_array initial_conds;
  
  if (y < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.4;
      // Velocity y
      initial_conds[3] = -2.0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.4;
      // Velocity y
      initial_conds[3] = 2.0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_3Y(double y)
{

  var_array initial_conds;
  
  if (y < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1000.0;
      // Velocity y
      initial_conds[3] = 0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.01;
      // Velocity y
      initial_conds[3] = 0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_4Y(double y)
{

  var_array initial_conds;
  
  if (y < 0.5)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.01;
      // Velocity y
      initial_conds[3] = 0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 100.0;
      // Velocity y
      initial_conds[3] = 0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1D_5Y(double y)
{

  var_array initial_conds;
  
  if (y < 0.5)
    {
      // Density
      initial_conds[0] = 5.99924;
      // Velocity
      initial_conds[1] = 0;
      // Pressure
      initial_conds[2] = 460.894;
      // Velocity y
      initial_conds[3] = 19.5975;
    }
  else
    {
      // Density
      initial_conds[0] = 5.99242;
      // Velocity
      initial_conds[1] = 0;
      // Pressure
      initial_conds[2] = 46.0950;
      // Velocity y
      initial_conds[3] = -6.19633;
    }

  return initial_conds;
}

var_array InitialCondTests::Test1_diag(double x, double y)
{

  var_array initial_conds;
  
  if (x > y)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity x
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1.0;
      // Velocity y
      initial_conds[3] = 0.0;
    }
  else
    {
      // Density
      initial_conds[0] = 0.125;
      // Velocity x
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.1;
      // Velocity y
      initial_conds[3] = 0.0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test2_diag(double x, double y)
{

  var_array initial_conds;
  
  if (x > y)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = -2.0;
      // Pressure
      initial_conds[2] = 0.4;
      // Velocity y
      initial_conds[3] = -2.0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 2.0;
      // Pressure
      initial_conds[2] = 0.4;
      // Velocity y
      initial_conds[3] = 2.0;
    }

  return initial_conds;
}


var_array InitialCondTests::Test3_diag(double x, double y)
{

  var_array initial_conds;
  
  if (x > y)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1000.0;
      // Velocity y
      initial_conds[3] = 0.0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.01;
      // Velocity y
      initial_conds[3] = 0.0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test4_diag(double x, double y)
{

  var_array initial_conds;
  
  if (x > y)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.01;
      // Velocity y
      initial_conds[3] = 0.0;
    }
  else
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 100.0;
      // Velocity y
      initial_conds[3] = 0.0;
    }

  return initial_conds;
}

var_array InitialCondTests::Test5_diag(double x, double y)
{

  var_array initial_conds;
  
  if (x > y)
    {
      // Density
      initial_conds[0] = 5.99924;
      // Velocity
      initial_conds[1] = 19.5975;
      // Pressure
      initial_conds[2] = 460.894;
      // Velocity y
      initial_conds[3] = 19.5975;
    }
  else
    {
      // Density
      initial_conds[0] = 5.99242;
      // Velocity
      initial_conds[1] = -6.19633;
      // Pressure
      initial_conds[2] = 46.0950;
      // Velocity y
      initial_conds[3] = -6.19633;
    }

  return initial_conds;
}

var_array InitialCondTests::CylindricalExplosion(double x, double y)
{

  var_array initial_conds;
  
  if (sqrt(pow(x-1,2) + pow(y-1,2))<0.4)
    {
      // Density
      initial_conds[0] = 1.0;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 1.0;
      // Velocity y
      initial_conds[3] = 0.0;
    }
  else
    {
      // Density
      initial_conds[0] = 0.125;
      // Velocity
      initial_conds[1] = 0.0;
      // Pressure
      initial_conds[2] = 0.1;
      // Velocity y
      initial_conds[3] = 0.0;
    }

  return initial_conds;
}

