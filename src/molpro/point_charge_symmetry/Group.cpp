#include "Group.h"
#include "CoordinateSystem.h"
#include <iostream>
#include <regex>

namespace molpro::point_charge_symmetry {

static CoordinateSystem s_default_coordinate_system;

Group::Group(CoordinateSystem& coordinate_system) : m_coordinate_system(coordinate_system) {}

// Group::Group(std::string name) : m_coordinate_system(s_group_default_coordinate_system), m_name(std::move(name)) {}

Group::Group() : m_coordinate_system(s_group_default_coordinate_system) {}

Group::Group(CoordinateSystem& coordinate_system, const Group& source)
    : m_coordinate_system(coordinate_system), m_name(source.m_name) {
  for (const auto& m : source.m_members) {
    m_members.emplace(m->clone(m_coordinate_system));
  }
}

Group::Group(CoordinateSystem& coordinate_system, std::string name, bool generators_only)
    : m_coordinate_system(coordinate_system), m_name(name) {
  //  std::cout << "Group::Group " << name << " " << generators_only << std::endl;
  using vec = CoordinateSystem::vec;
  const vec zaxis{0, 0, 1};
  const auto all = not generators_only;
  if (all)
    add(Identity());

  std::smatch m;

  if (all and (std::regex_match(name, m, std::regex{"[CD][0-9]*[02468]h"}) or
               std::regex_match(name, m, std::regex{"D[0-9]*[13579]d"}) or
               std::regex_match(name, m, std::regex{"Ih|Th|Oh|Dinfh"})))
    add(Inversion());
  if (name == "Ci")
    add(Inversion());

  if (name[0] == 'I') {
    auto gold = (1 + std::sqrt(double(5))) / 2;
    std::vector<vec> pentagons;
    for (int axis = 0; axis < 3; axis++) {
      auto yaxis = Eigen::Matrix3d::Identity().col((axis + 1) % 3);
      auto zaxis = Eigen::Matrix3d::Identity().col((axis + 2) % 3);
      for (int flip = -1; flip < 2; flip += 2) {
        pentagons.push_back((yaxis * flip + gold * zaxis) / std::sqrt(2 + gold));
        if (pentagons.back().dot(pentagons.front()) < 0)
          pentagons.back() = -pentagons.back();
        for (int count = 1; count < 6; count++) {
          add(Rotation(pentagons.back(), 5, true, count));
          if (name == "Ih")
            add(Rotation(pentagons.back(), 10, false, 2 * count - 1));
        }
        //                std::cout << "pentagon "<<pentagons.back().transpose()<<std::endl;
      }
    }
    std::swap(pentagons[1], pentagons[2]); // so that the 5 satellite pentagons are in order
                                           //    for (const auto& p1:pentagons) {
                                           //      std::cout <<"pentagon overlaps";
                                           //      for (const auto& p2:pentagons)
                                           // std::cout << " "<<p2.dot(p1);
    //      std::cout << std::endl;
    //    }
    if (all) {
      for (int i = 0; i < 5; i++) {
        vec bisector, joiner;
        bisector = (pentagons[0] + pentagons[i + 1]) / std::sqrt(double(2));
        joiner = (pentagons[0] - pentagons[i + 1]);
        add(Rotation(bisector, 2));
        if (name == "Ih")
          add(Reflection(joiner.cross(bisector)));
        bisector = (pentagons[(i + 1) % 5 + 1] + pentagons[i + 1]) / std::sqrt(double(2));
        joiner = (pentagons[(i + 1) % 5 + 1] - pentagons[i + 1]);
        add(Rotation(bisector, 2));
        if (name == "Ih")
          add(Reflection(joiner.cross(bisector)));
        bisector = (-pentagons[(i + 2) % 5 + 1] + pentagons[i + 1]) / std::sqrt(double(2));
        joiner = (-pentagons[(i + 2) % 5 + 1] - pentagons[i + 1]);
        add(Rotation(bisector, 2));
        if (name == "Ih")
          add(Reflection(joiner.cross(bisector)));
      }
      for (int i = 0; i < 5; i++)
        for (int count = 1; count < 3; count++) {
          add(Rotation((pentagons[0] + pentagons[i + 1] + pentagons[(i + 1) % 5 + 1]) / std::sqrt(double(3)), 3, true,
                       count));
          if (name == "Ih")
            add(Rotation((pentagons[0] + pentagons[i + 1] + pentagons[(i + 1) % 5 + 1]) / std::sqrt(double(3)), 6,
                         false, 4 * count - 3));
          add(Rotation((-pentagons[(i + 3) % 5 + 1] + pentagons[i + 1] + pentagons[(i + 1) % 5 + 1]) /
                           std::sqrt(double(3)),
                       3, true, count));
          if (name == "Ih")
            add(Rotation((-pentagons[(i + 3) % 5 + 1] + pentagons[i + 1] + pentagons[(i + 1) % 5 + 1]) /
                             std::sqrt(double(3)),
                         6, false, 4 * count - 3));
        }
    }
    //    std::cout <<*this<<std::endl;
  }
  if (name[0] == 'O') {
    for (int axis = 0; axis < 3; axis++) {
      for (int count = 1; count < 4; count += 2) {
        if (all or (name == "O" and axis == 0 and count == 1))
          add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, true, count));
        if (name == "Oh" and (all or (axis == 0 and count == 1)))
          add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, count));
      }
      if (all)
        add(Rotation(Eigen::Matrix3d::Identity().col(axis), 2, true, 1));
      if (name == "Oh" and all)
        add(Reflection(Eigen::Matrix3d::Identity().col(axis)));
      Eigen::Vector3d vec{1, 1, 1};
      vec(axis) = 0;
      if (all)
        add(Rotation(vec, 2));
      if (name == "Oh" and all)
        add(Reflection(vec));
      vec((axis + 1) % 3) = -1;
      if (all)
        add(Rotation(vec, 2));
      if (name == "Oh" and all)
        add(Reflection(vec));
    }
    for (int corner = 0; corner < (all ? 4 : 2); corner++) {
      auto sq2 = std::sqrt(1 / double(2));
      for (int count = 1; count < 3; count++) {
        const vec& axis = vec{std::cos((2 * corner + 1) * acos(double(-1)) / 4),
                              std::sin((2 * corner + 1) * acos(double(-1)) / 4), sq2};
        if (all or (name == "O" and corner == 0 and count == 1))
          add(Rotation(axis, 3, true, count));
        if (name == "Oh" and (all or (corner == 0 and count == 1)))
          add(Rotation(axis, 6, false, count * 4 - 3));
      }
    }
  }
  if (name[0] == 'T') {
    for (int axis = 0; axis < 3; axis++) {
      if (name == "Th" and (all or axis == 0))
        add(Reflection(Eigen::Matrix3d::Identity().col(axis)));
      if (all or (name == "T" and axis == 0))
        add(Rotation(Eigen::Matrix3d::Identity().col(axis), 2, true, 1));
      if (name == "Td" and (all or axis == 0))
        add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, 1));
      if (name == "Td" and all)
        add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, 3));
      if (name == "Td" and all)
        add(Reflection(Eigen::Matrix3d::Identity().col((axis + 1) % 3) +
                       Eigen::Matrix3d::Identity().col((axis + 2) % 3)));
      if (name == "Td" and all)
        add(Reflection(Eigen::Matrix3d::Identity().col((axis + 1) % 3) -
                       Eigen::Matrix3d::Identity().col((axis + 2) % 3)));
    }
    for (int corner = 0; corner < (all ? 4 : 2); corner++) {
      auto sq2 = std::sqrt(1 / double(2));
      for (int count = 1; count < 3; count++) {
        const vec axis = vec{std::cos((2 * corner + 1) * acos(double(-1)) / 4),
                             std::sin((2 * corner + 1) * acos(double(-1)) / 4), sq2};
        if (all or (corner == 0 and count == 1))
          add(Rotation(axis, 3, true, count));
        if (name == "Th" and all)
          add(Rotation(axis, 6, false, count * 4 - 3));
      }
    }
  }
  if (name == "Cinfv" or name == "Dinfh") { // representative only
    m_name = name;
    name.replace(1, 3, "10");
  }

  if (std::regex_match(name, m, std::regex{"Dinfh|Cs|[CD][1-9][0-9]*h"}))
    add(Reflection(zaxis));

  if (std::regex_match(name, m, std::regex{"([C])([1-9][0-9]*)([v])"}) or
      std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)([d])"}) or
      std::regex_match(name, m, std::regex{"([D])([1-9]*[02468])([h])"})) {
    auto order = std::stoi(m.str(2));
    auto angle = std::acos(double(-1)) / order;
    for (int count = 0; count < (all ? order : 1); count++)
      add(Reflection({std::cos(count * angle), std::sin(count * angle), 0}));
  }
  if (std::regex_match(name, m, std::regex{"([D])([1-9]*[13579])([h])"})) {
    auto order = std::stoi(m.str(2));
    auto angle = std::acos(double(-1)) / order;
    for (int count = 0; count < (all ? order : 1); count++)
      add(Reflection({std::cos((count + double(0.5)) * angle), std::sin((count + double(0.5)) * angle), 0}));
  }

  if (std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)([^v])?"})) {
    auto order = std::stoi(m.str(2));
    //    for (double angle = 0; angle < std::acos(double(-1)) - 1e-10; angle += std::acos(double(-1)) / order)
    for (int count = 0; count < (all ? order : (std::regex_match(name, m, std::regex{"D[0-9]*h"}) ? 0 : 1)); count++) {
      //      std::cout << "count=" << count << std::endl;
      double angle = (std::regex_match(name, m, std::regex{"D[0-9]*[02468]d"}) ? (count + 0.5) : count) *
                     std::acos(double(-1)) / order;
      //      std::cout << "adding C2 for " << name << ", angle=" << angle * 180 / std::acos(double(-1)) << std::endl;
      add(Rotation({std::cos(angle), std::sin(angle), 0}, 2));
    }
  }

  if (std::regex_match(name, m, std::regex{"([CD])([1-9][0-9]*)([hvd]*)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 1; count < (all ? order : 2); count++)
      add(Rotation(zaxis, order, true, count));
  }
  if (all and std::regex_match(name, m, std::regex{"([CD])([0-9]*[02468])(h)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < (all ? order / 2 : 1); count++)
      if (count * 2 + 1 != order / 2)
        add(Rotation(zaxis, order, false, count * 2 + 1));
    if (order > 2)
      for (int count = 0; count < ((order % 4) ? order / 2 : order / 4); count++)
        if (count * 2 + 1 != ((order % 4) ? order / 2 : order / 4))
          add(Rotation(zaxis, order / 2, false, count * 2 + 1));
  }
  if (all and std::regex_match(name, m, std::regex{"([CD])([0-9]*[13579])(h)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order; count++)
      if (count != (order - 1) / 2)
        add(Rotation(zaxis, order, false, count * 2 + 1));
  }
  if (std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)(d)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order; count++)
      if (2 * count != (order - 1))
        add(Rotation(zaxis, 2 * order, false, count * 2 + 1));
  }
  if (std::regex_match(
          name, m,
          std::regex{
              "([S])([0-9]*[02468])"})) { // could be done more prettily with explicit proper rotations and inversion
    auto order = std::stoi(m.str(2));
    for (int count = 1; count < order; count++)
      if (count % 2)
        add(Rotation(zaxis, order, false, count));
      else if (count == order / 2)
        add(Inversion());
      else
        add(Rotation(zaxis, order / 2, true, count / 2));
  }
}

Group::Group(const std::string& name, bool generators_only)
    : Group(s_default_coordinate_system, name, generators_only) {}

Rotation Group::highest_rotation(bool proper, size_t index) const {
  size_t count = 0;
  int order = 0;
  for (const auto& member : m_members)
    if (member->order() > order && (member->proper() or not proper))
      order = member->order();
  if (order > 0)
    for (const auto& member : m_members)
      if (member->order() == order && member->count() == 1 && (member->proper() or not proper) && index == count++)
        return dynamic_cast<Rotation&>(*member);
  return Rotation({0, 0, 1}, 1);
}

Group generate(const Group& generator) {
  //  auto g = generator;
  Group g(generator.name(), true);
  g.add(Identity());
  for (size_t last_size = 0; last_size < g.size();) {
    last_size = g.size();
    //    std::cout << "generate iterate last_size="<<last_size<<std::endl;
    const auto oldg = Group(const_cast<CoordinateSystem&>(g.coordinate_system()), g);
    //    std::cout << "oldg:\n"<<oldg<<std::endl;
    for (const auto& o1 : oldg)
      for (const auto& o2 : oldg)
        g.add((*o1) * (*o2));
    //    std::cout << "generate iterate new size="<<g.size()<<std::endl;
    //    std::cout << "group at end of iteration\n"<<g<<std::endl;
  }
  return g;
}

const Operator& Group::operator[](size_t index) const {
  size_t i = 0;
  for (const auto& member : m_members) {
    if (++i > index)
      return *member;
  }
  throw std::runtime_error("index out of range");
}
} // namespace molpro::point_charge_symmetry
