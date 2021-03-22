#include "CoordinateSystem.h"
#include <unsupported/Eigen/MatrixFunctions>

CoordinateSystem::CoordinateSystem(const vec &origin, const mat &axes)
{
  auto generator = axes.log();
  axis_generator()(0)=generator(1,0);
  axis_generator()(1)=generator(2,0);
  axis_generator()(2)=generator(2,1);
  this->origin()(0)=origin(0);
  this->origin()(1)=origin(1);
  this->origin()(2)=origin(2);
}
const CoordinateSystem::mat CoordinateSystem::axes() const {
  Eigen::Matrix3d generator;
  generator << 0,axis_generator()(0),axis_generator()(1),
  -axis_generator()(0),0,axis_generator()(2),
  -axis_generator()(1),-axis_generator()(2),0;
  return generator.exp();
}
