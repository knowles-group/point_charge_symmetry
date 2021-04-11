#include "Group.h"
#include <molpro/point_charge_symmetry/CoordinateSystem.h>
#include <regex>

namespace molpro::point_charge_symmetry {

static CoordinateSystem s_default_coordinate_system;

Group::Group(CoordinateSystem& coordinate_system) : m_coordinate_system(coordinate_system) {}

// Group::Group(std::string name) : m_coordinate_system(s_group_default_coordinate_system), m_name(std::move(name)) {}

Group::Group() : m_coordinate_system(s_group_default_coordinate_system) {}

Group::Group(CoordinateSystem& coordinate_system, const Group& source)
    : m_coordinate_system(coordinate_system), m_name(source.m_name) {
  for (const auto& m : source.m_members) {
    m_members.emplace_back(m->clone(m_coordinate_system));
  }
}

Group::Group(CoordinateSystem& coordinate_system, std::string name, bool generators_only)
    : m_coordinate_system(coordinate_system) {
  auto& g = *this;
  using vec = CoordinateSystem::vec;
  const vec xaxis{1, 0, 0};
  const vec yaxis{0, 1, 0};
  const vec zaxis{0, 0, 1};
  const auto all = not generators_only;
  g.name() = name;
  if (all)
    g.add(Identity());
  std::smatch m;
  if (name == "Ih")
    g.add(Rotation(zaxis, 7)); // TODO implement
  if (name == "Oh") {
    for (int axis = 0; axis < 3; axis++) {
      for (int count = 0; count < 3; count += 2) {
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, true, count));
        if (all)
          g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, count));
      }
      g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 2, true, 1));
      g.add(Reflection(Eigen::Matrix3d::Identity().col(axis)));
    }
    for (int corner = 0; corner < (all ? 4 : 2); corner++) {
      auto sq2 = std::sqrt(1 / double(2));
      g.add(Rotation(vec{std::cos((2 * corner + 1) * acos(double(-1)) / 4),
                         std::sin((2 * corner + 1) * acos(double(-1)) / 4), sq2},
                     3));
    }
  }
  if (name == "Td") {
    for (int axis = 0; axis < 3; axis++) {
      if (all)
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 2, true, 1));
      if (all)
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, 1));
      if (all)
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, 3));
      if (all)
        g.add(Reflection(Eigen::Matrix3d::Identity().col((axis + 1) % 3) +
            Eigen::Matrix3d::Identity().col((axis + 2) % 3)));
      g.add(Reflection(Eigen::Matrix3d::Identity().col((axis + 1) % 3) -
          Eigen::Matrix3d::Identity().col((axis + 2) % 3)));
    }
    for (int corner = 0; corner < (all ? 4 : 2); corner++) {
      auto sq2 = std::sqrt(1 / double(2));
      for (int count = 1; count < 3; count++)
        g.add(Rotation(vec{std::cos((2 * corner + 1) * acos(double(-1)) / 4),
                           std::sin((2 * corner + 1) * acos(double(-1)) / 4), sq2},
                       3, true, count));
    }
  }
  if (name == "Cinfv" or name == "Dinfh") { // representative only
    m_name = name;
    name.replace(1, 3, "11");
  }

  if (std::regex_match(name, m, std::regex{"[CD][0-9]*[02468]h"}) or
      std::regex_match(name, m, std::regex{"D[0-9]*[13579]d"}) or
      std::regex_match(name, m, std::regex{"Ih|Th|Oh|Dinfh|Ci"}))
    g.add(Inversion());

  if (std::regex_match(name, m, std::regex{"Dinfh|Ih|Cs|[CD][1-9][0-9]*h"}))
    g.add(Reflection(zaxis));

  if (std::regex_match(name, m, std::regex{"([C])([1-9][0-9]*)([v])"}) or
      std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)([hd])"})) {
    auto order = std::stoi(m.str(2));
    auto angle = std::acos(double(-1)) / order;
    for (int count = 0; count < order; count++)
      g.add(Reflection({std::cos(count * angle), std::sin(count * angle), 0}));
  }

  if (std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)([^v])?"})) {
    auto order = std::stoi(m.str(2));
    //    for (double angle = 0; angle < std::acos(double(-1)) - 1e-10; angle += std::acos(double(-1)) / order)
    for (int count = 0; count < (all ? order : 2); count++) {
      double angle = (std::regex_match(name, m, std::regex{"D[0-9]*[02468]d"}) ? (count + 0.5) : count) *
          std::acos(double(-1)) / order;
      g.add(Rotation({std::cos(angle), std::sin(angle), 0}, 2));
    }
  }

  if (std::regex_match(name, m, std::regex{"([CD])([1-9][0-9]*)([hvd]*)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 1; count < order; count++)
      g.add(Rotation(zaxis, order, true, count));
  }
  if (std::regex_match(name, m, std::regex{"([CD])([0-9]*[02468])(h)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order / 2; count++)
      if (count * 2 + 1 != order / 2)
        g.add(Rotation(zaxis, order, false, count * 2 + 1));
    if (order > 2)
      for (int count = 0; count < ((order % 4) ? order / 2 : order / 4); count++)
        if (count * 2 + 1 != ((order % 4) ? order / 2 : order / 4))
          g.add(Rotation(zaxis, order / 2, false, count * 2 + 1));
  }
  if (std::regex_match(name, m, std::regex{"([CD])([0-9]*[13579])(h)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order; count++)
      if (count != (order - 1) / 2)
        g.add(Rotation(zaxis, order, false, count * 2 + 1));
  }
  if (std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)(d)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order; count++)
      if (2 * count != (order - 1))
        g.add(Rotation(zaxis, 2 * order, false, count * 2 + 1));
  }
  if (std::regex_match(
      name, m,
      std::regex{
          "([S])([0-9]*[02468])"})) { // could be done more prettily with explicit proper rotations and inversion
    auto order = std::stoi(m.str(2));
    for (int count = 1; count < order; count++)
      g.add(Rotation(zaxis, order, false, count));
  }
}

Group::Group(const std::string& name, bool generators_only) : Group(s_default_coordinate_system, name, generators_only) {}

} // namespace molpro::point_charge_symmetry
