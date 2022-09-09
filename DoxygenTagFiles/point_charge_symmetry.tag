<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.1">
  <compound kind="struct">
    <name>molpro::point_charge_symmetry::Atom</name>
    <filename>structmolpro_1_1point__charge__symmetry_1_1Atom.html</filename>
    <member kind="function">
      <type></type>
      <name>Atom</name>
      <anchorfile>structmolpro_1_1point__charge__symmetry_1_1Atom.html</anchorfile>
      <anchor>a6df76defae80f8ae60bb8a79a6fad011</anchor>
      <arglist>(Eigen::Vector3d r, double q, std::string name=&quot;&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>Eigen::Vector3d</type>
      <name>position</name>
      <anchorfile>structmolpro_1_1point__charge__symmetry_1_1Atom.html</anchorfile>
      <anchor>a09bdd998db8a2817ee4b29c72c29b097</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>charge</name>
      <anchorfile>structmolpro_1_1point__charge__symmetry_1_1Atom.html</anchorfile>
      <anchor>aef07128f455f9b14c2f86447ae140e42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>structmolpro_1_1point__charge__symmetry_1_1Atom.html</anchorfile>
      <anchor>ac98125bb9701e93fbcf479a9d6053216</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::CoordinateSystem</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</filename>
    <member kind="typedef">
      <type>std::array&lt; double, 6 &gt;</type>
      <name>parameters_t</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a0766b74cbd6d7967fa16a8a2fe2e75cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Eigen::Vector3d</type>
      <name>vec</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a74aafae6b66ae6f6d409439efc2513a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Eigen::Matrix3d</type>
      <name>mat</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a9df61c083097e7ea89d5ad62b9d2bed5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>CoordinateSystem &amp;</type>
      <name>operator=</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a6bc15950ac3015351f84df56b0cd1a0f</anchor>
      <arglist>(const CoordinateSystem &amp;source)</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Map&lt; vec &gt;</type>
      <name>origin</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a9b039a24e349034c374b9cb7f4d1d05b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Map&lt; const vec &gt;</type>
      <name>origin</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a21c11f56c490c4a40bd39b633efafd7e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const mat</type>
      <name>axes</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>af09e7dbfdf94bbaad095b78584f95b3d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const std::array&lt; mat, 3 &gt;</type>
      <name>axes_gradient</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>aabe61582198af046f7e52d6526796cf5</anchor>
      <arglist>(int displacements=0, double step=1e-4) const</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Map&lt; const vec &gt;</type>
      <name>axis_generator</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a484319a969de8edfdb832149715f5b0f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Map&lt; vec &gt;</type>
      <name>axis_generator</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>ac74eb7ad54bed486b7f244d34154395a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CoordinateSystem</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a5342335452fc787affb699d87ff3b4f9</anchor>
      <arglist>(const RotationParameterType rotation_parameter_type=RotationParameterType::Euler, const vec &amp;origin=vec::Zero(), const mat &amp;axes=mat::Identity())</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>from_axes</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>afcfbe64cf72c1d7035fff172f6de2761</anchor>
      <arglist>(const mat &amp;axes=mat::Identity()) const</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>data</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a4ccd1e6888ba516e390de1cf9bdd480a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const double *</type>
      <name>data</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>aae13ca05b633b3e1b3f181de62ffdf22</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a8113cdab9cf926e6b9b42982a5feb4d6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>vec</type>
      <name>to_local</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a06e7cc9e3ccbd0ebbfebcbbe53172f7b</anchor>
      <arglist>(const vec &amp;source) const</arglist>
    </member>
    <member kind="function">
      <type>vec</type>
      <name>to_global</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a319038a09cc3f37b65d4bc1a5c9d590f</anchor>
      <arglist>(const vec &amp;source) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cycle_axes</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a96581e272ac7324db6d3c667b9555da8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::array&lt; std::array&lt; double, 2 &gt;, 3 &gt;</type>
      <name>rotation_generator_ranges</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a4359f2e7d6f87ef5a7f00a44df7d03ad</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable">
      <type>parameters_t</type>
      <name>m_parameters</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a84ddd78d489ed7b9cee5299f01b13c06</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const RotationParameterType</type>
      <name>m_rotation_parameter_type</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a6f25e0750bf7b0885e868de574e514ac</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>m_rotation_parameter_scale</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>a880dd02ce8bf89ec2efa20e9ad5c0a9b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>m_axis_permutation_rot90_next</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1CoordinateSystem.html</anchorfile>
      <anchor>ad19380dc340caec7d84a2cc9d1f1f837</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Group</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Group.html</filename>
    <member kind="typedef">
      <type>std::vector&lt; std::unique_ptr&lt; Operator &gt; &gt;::iterator</type>
      <name>iterator</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a9d152eb6590e2c5b211c11c18ff2f78e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::vector&lt; std::unique_ptr&lt; Operator &gt; &gt;::const_iterator</type>
      <name>const_iterator</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a1f3dc95b58a5978ec3caa68d31815ae4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>aefb9e01d5cf09d0fdbcaea9fe2b4b126</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a86f39b7dd54d4ce75f732ff6e3386395</anchor>
      <arglist>(CoordinateSystem &amp;coordinate_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a7834b7b2403ca79706e755442cef1c79</anchor>
      <arglist>(CoordinateSystem &amp;coordinate_system, std::string name, bool generators_only=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a948236b5980d157955cd311089686f9d</anchor>
      <arglist>(const std::string &amp;name, bool generators_only=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a7506a49596064a77f008cd9edeb3dcf5</anchor>
      <arglist>(CoordinateSystem &amp;coordinate_system, const Group &amp;source)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>aebca0e169a19eaa8f71c1e6a9171660e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::string &amp;</type>
      <name>name</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a7e43c4589f0329fab245740920def9c1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a10fae350429d0a75617918427000c390</anchor>
      <arglist>(Identity op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>aba563349d133957cdd618617c8efe1a2</anchor>
      <arglist>(Inversion op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a5c10250502abddcd11c40f6649d3f391</anchor>
      <arglist>(Reflection op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a68149c054d363b983995c0e24d3b0c24</anchor>
      <arglist>(const Rotation &amp;op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>add7e62eb589fc689f8442409be4cac9f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const CoordinateSystem &amp;</type>
      <name>coordinate_system</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a9b60594c0f1f9583cf19bdff8e210b7d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>CoordinateSystem::parameters_t &amp;</type>
      <name>coordinate_system_parameters</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a9191a4d4de96864f043a17624cf3314c</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>iterator</type>
      <name>begin</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a8ae4d5befe0e5d5bf06667d644a237bf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const_iterator</type>
      <name>begin</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a7100e1fd9ed9b375735ad2647cfaca7b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>iterator</type>
      <name>end</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a826ef92bafca26beb776822ac16bd8d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const_iterator</type>
      <name>end</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a3b5d41ef30fa818d0a6b5c94ab94df7c</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>size</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>aeb652ce8fdcbab1590cad8c2cbb4271a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Rotation</type>
      <name>highest_rotation</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a619e53a4e820b38cb96593c519b569ed</anchor>
      <arglist>(bool proper=true, size_t index=0) const</arglist>
    </member>
    <member kind="function">
      <type>Operator &amp;</type>
      <name>operator[]</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a497684b31e896cd58968d599523c9197</anchor>
      <arglist>(size_t index)</arglist>
    </member>
    <member kind="function">
      <type>const Operator &amp;</type>
      <name>operator[]</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a182234782fe1755a0c9d70186a445769</anchor>
      <arglist>(size_t index) const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>CoordinateSystem &amp;</type>
      <name>m_coordinate_system</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a5981a8cec50271e54f2dbf0db83391a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; std::unique_ptr&lt; Operator &gt; &gt;</type>
      <name>m_members</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>aaac4a21c988a5e1e2c1c9b1b46c21757</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>m_name</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Group.html</anchorfile>
      <anchor>a862e78d8c77beb38f4e9d4d0689db638</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Identity</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Identity.html</filename>
    <base>molpro::point_charge_symmetry::Operator</base>
    <member kind="function">
      <type></type>
      <name>Identity</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Identity.html</anchorfile>
      <anchor>a0bcd073fc741ca424b9692cef2a07287</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Identity</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Identity.html</anchorfile>
      <anchor>a8f2f9c160c9d1bf879c4ebbaafce636d</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system)</arglist>
    </member>
    <member kind="function">
      <type>vec</type>
      <name>operator_local</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Identity.html</anchorfile>
      <anchor>af2fc7f3c83b0bda5d1918755f34ad843</anchor>
      <arglist>(vec v) const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Identity.html</anchorfile>
      <anchor>a1e9014f64bc42a8ad0ed63af7426d36b</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Identity.html</anchorfile>
      <anchor>ae658e6423038e7195fcc3e4102636880</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system) const override</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Identity.html</anchorfile>
      <anchor>ab38079fe54be3bccff340d7db758d530</anchor>
      <arglist>(const std::string &amp;title, bool coordinate_frame=false) const override</arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Identity.html</anchorfile>
      <anchor>a2697825715974a353728f0d4d5658112</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Inversion</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</filename>
    <base>molpro::point_charge_symmetry::Operator</base>
    <member kind="function">
      <type></type>
      <name>Inversion</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>a39d3bf7828e39af47d18230bce5627d1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Inversion</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>a8353b36e227997d1ccf43d26279313da</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system)</arglist>
    </member>
    <member kind="function">
      <type>vec</type>
      <name>operator_local</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>ab3704a78409fcb14a3743e5f1b8b3e57</anchor>
      <arglist>(vec v) const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>a06010c3f3e6fb780766f2e5c52bd4cd3</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>aff1965ad7b0063c6ef565909d7666228</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system) const override</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>a3b90a2ab1031de84c97e9979358c6e33</anchor>
      <arglist>(const std::string &amp;title, bool coordinate_frame=false) const override</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>proper</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>aaeacadfbfe31488d166285746cdb5f3f</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Inversion.html</anchorfile>
      <anchor>a2697825715974a353728f0d4d5658112</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Molecule</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</filename>
    <member kind="function">
      <type></type>
      <name>Molecule</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>aa11f66935c6058ab9b6e5a467b58f3b3</anchor>
      <arglist>(const std::string &amp;filename)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Molecule</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a2fbbb47c6bd3dd56a461d24a95c73c0f</anchor>
      <arglist>(const Eigen::MatrixXd &amp;coordinates, const Eigen::VectorXd &amp;charges)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>aecaeec41b7c88521b5ebe791184ef2b8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Vector3d</type>
      <name>centre_of_charge</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a847def33d738272a7303eeeaed44254f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Matrix3d</type>
      <name>inertia_tensor</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a9ea4c6c4df5e88a2e5b54659b7598a81</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Matrix3d</type>
      <name>inertial_axes</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a8e7b4c1ea465cc5c42e0de28c4319a05</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a2d9cf5b90d93d8afe1fcb1c229cff003</anchor>
      <arglist>(const std::string &amp;filename, const std::string &amp;title=&quot;&quot;, const std::string &amp;format=&quot;xyz&quot;)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>size</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a3725f821fbd7f57382a41914706d3986</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Vector3d</type>
      <name>findaxis</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a9b45d61102dc23275ed13f34656f984c</anchor>
      <arglist>(int order) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>randomise</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a0e82308d19831fa7e0eb680a4f06feba</anchor>
      <arglist>(double amplitude)</arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Atom &gt;</type>
      <name>m_atoms</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a203f22c659931fd192607af4c18b46e9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>m_title</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Molecule.html</anchorfile>
      <anchor>a4cc69ad89625abb8645a4d3b5d2344b8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Operator</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Operator.html</filename>
    <member kind="typedef">
      <type>CoordinateSystem::vec</type>
      <name>vec</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>ae3272314c978fb9df9d4b88ab4c67580</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>CoordinateSystem::mat</type>
      <name>mat</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a1ce1b799e1e3c52486733f12880a645d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Operator</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>ae1a861c35385f2198a9422139e4ad180</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Operator</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>ae7e5408788aacb588a806f52ec73a301</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Operator</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a48ca02743a6b1847aac84a1686128864</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system)</arglist>
    </member>
    <member kind="function">
      <type>vec</type>
      <name>operator()</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a4365a6ab9e14c4fcec66ddf16ebdbb45</anchor>
      <arglist>(vec v) const</arglist>
    </member>
    <member kind="function">
      <type>std::array&lt; Operator::vec, 6 &gt;</type>
      <name>operator_gradient</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a3bfef8b36616aed81cd6558aa99a5e32</anchor>
      <arglist>(vec v, int numerical=0, double step=2e-3) const</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual vec</type>
      <name>operator_local</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>aaac505b23c26e8ddb095bd142020e479</anchor>
      <arglist>(vec v) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a87cc4d2be30529bae58e3902cb0d590b</anchor>
      <arglist>(const std::string &amp;title, bool coordinate_frame=false) const</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>name</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a60fa2a5a9884f41d2766f643a272f1ad</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>ae2a492ae9030d34065024229dda232d4</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a0e703885947bb8888e32a6b370ec769d</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system) const =0</arglist>
    </member>
    <member kind="function">
      <type>const CoordinateSystem &amp;</type>
      <name>coordinate_system</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>ade097fa711cf69417ce54d85db305c79</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>order</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a1e9e5fcdd66cfac70ed71da13e60851b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>proper</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a29cd6ae83b0de6082586b7022afcb858</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const CoordinateSystem &amp;</type>
      <name>m_coordinate_system</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a99f00f9ab327c17009812c7253b24c71</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>m_name</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a7af58a306d03c37434b273430bcdb605</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="protected">
      <type>friend class</type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Operator.html</anchorfile>
      <anchor>a2697825715974a353728f0d4d5658112</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Problem_optimise_frame</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Problem__optimise__frame.html</filename>
    <member kind="function">
      <type></type>
      <name>Problem_optimise_frame</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Problem__optimise__frame.html</anchorfile>
      <anchor>a1f9fbe2839aac28a059050710387cd91</anchor>
      <arglist>(SymmetryMeasure &amp;sm)</arglist>
    </member>
    <member kind="function">
      <type>value_t</type>
      <name>residual</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Problem__optimise__frame.html</anchorfile>
      <anchor>aeae31fd6467c816d2bc02580bbd08a83</anchor>
      <arglist>(const container_t &amp;parameters, container_t &amp;residual) const override</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>diagonals</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Problem__optimise__frame.html</anchorfile>
      <anchor>ac7f0aef0e5a1fe6423619f283594a0da</anchor>
      <arglist>(container_t &amp;d) const override</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Projector</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Projector.html</filename>
    <member kind="function">
      <type></type>
      <name>Projector</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Projector.html</anchorfile>
      <anchor>a053a1e3b0e4f1caa137d9e33bb5c8e2d</anchor>
      <arglist>(const Group &amp;group, const Molecule &amp;molecule)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Projector</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Projector.html</anchorfile>
      <anchor>a8d836c34d3701325232089c61741671c</anchor>
      <arglist>()=delete</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; double &gt;</type>
      <name>symmetric</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Projector.html</anchorfile>
      <anchor>a2b2c8d8c64043cf2afd070f8c3605018</anchor>
      <arglist>(std::vector&lt; double &gt; vector) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove_symmetric</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Projector.html</anchorfile>
      <anchor>aec7ce2396aa654c345b3b153d790f551</anchor>
      <arglist>(std::vector&lt; double &gt; vector) const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const size_t</type>
      <name>m_n3</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Projector.html</anchorfile>
      <anchor>a16015db0cf1df216f57a110447dd7085</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; double &gt;</type>
      <name>m_V</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Projector.html</anchorfile>
      <anchor>a2b7e1cf7e3df03798c0a3f148be59325</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Reflection</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</filename>
    <base>molpro::point_charge_symmetry::Operator</base>
    <member kind="function">
      <type></type>
      <name>Reflection</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a871aa930da07e30d2caef567a70e62b1</anchor>
      <arglist>(vec normal)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Reflection</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a7c31f984e13160b4c5f94f110d74ad3a</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system, vec normal)</arglist>
    </member>
    <member kind="function">
      <type>vec</type>
      <name>operator_local</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a8c50095d94aa1d2b6eee1857dfebeb55</anchor>
      <arglist>(vec v) const override</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a9c0834a79af715b56c38bbea055cb56b</anchor>
      <arglist>(const std::string &amp;title, bool coordinate_frame=false) const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a3063061df0f5e524c2c1eab628f57556</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a8ae983c77f944f2b289db0921688fa90</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system) const override</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>proper</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a4f8d056747846cae1548cd3684107d8b</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vec</type>
      <name>m_normal</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>ad903d80f4f45d470eec7abab9277a4af</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Reflection.html</anchorfile>
      <anchor>a2697825715974a353728f0d4d5658112</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::Rotation</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</filename>
    <base>molpro::point_charge_symmetry::Operator</base>
    <member kind="function">
      <type></type>
      <name>Rotation</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a67fbd2680c9b5d1fa67ffb0d6d1b6981</anchor>
      <arglist>(vec axis, int order=2, bool proper=true, int count=1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Rotation</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a6e103d5cfb4aeeb56808c72b5c5107c4</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system, vec axis, int order=2, bool proper=true, int count=1)</arglist>
    </member>
    <member kind="function">
      <type>vec</type>
      <name>operator_local</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a8513dd1f3bed673a91ed7a24cb6a6994</anchor>
      <arglist>(vec v) const override</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>acf7118318038b3dba385bac53b833fe9</anchor>
      <arglist>(const std::string &amp;title, bool coordinate_frame=false) const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>acdb6f3d5ef624f5ba496e12ba768d32c</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>Operator *</type>
      <name>clone</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>afe898a83df2f11094b98f04d953810dd</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system) const override</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>order</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>ac090990e3f79ae6119dec12f37a740ed</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>proper</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a48dcae7888fa5c8020895ef36d16634c</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>const vec &amp;</type>
      <name>axis</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a96cb927ff508b7a8a3a47ec8a8e3741b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vec</type>
      <name>m_axis</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>abb6ee136f2a0bbb5f657ebe8c56bc334</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>m_order</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>ab8dc095c03eb0f72db997dc95b9ce161</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>m_proper</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a07f01f5be2969b39a601b2969715fad5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>m_count</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a4ef3446238b6e82bfbf1a555acdbe319</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>Group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1Rotation.html</anchorfile>
      <anchor>a2697825715974a353728f0d4d5658112</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::point_charge_symmetry::SymmetryMeasure</name>
    <filename>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</filename>
    <member kind="function">
      <type></type>
      <name>SymmetryMeasure</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a82c2fdf8da032b3cf7a14b1456312316</anchor>
      <arglist>(const Molecule &amp;molecule, const Group &amp;group)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset_neighbours</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a8023ad7c62357750b4eb7c5f6faa3ed6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a3a7f99ee3bd2e4ec336805282220d047</anchor>
      <arglist>(int operator_index=-1, int functional_form=0, int verbosity=-1) const</arglist>
    </member>
    <member kind="function">
      <type>CoordinateSystem::parameters_t</type>
      <name>coordinate_system_gradient</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a51b1b5744abb7f2b75b404b0456b55a3</anchor>
      <arglist>(int operator_index=-1, int functional_form=0) const</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; double &gt;</type>
      <name>atom_gradient</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a379514c18f85459fa850ebec2c6fc851</anchor>
      <arglist>(int operator_index=-1, int functional_form=0) const</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>ad5c5038a81d74ce75c4c8d22c0237ca4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>adopt_inertial_axes</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>aee71c982c7898266708bb6c540898045</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>refine_frame</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>aa1db904a538e1746b2f3c8f71d573ac6</anchor>
      <arglist>(int verbosity=-1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>spherical_top</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>ad15003f9f3540e24f0439b9ee9bf733c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>symmetric_top</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a3b89cbe811f16c26c874ac9643799a6c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Molecule</type>
      <name>refine</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a0f1a3705f25b2fd01ad7f6ffba8da86b</anchor>
      <arglist>(double distance_penalty=0.001, bool project=false, int repeat=1) const</arglist>
    </member>
    <member kind="function">
      <type>CoordinateSystem::vec</type>
      <name>inertia_principal_values</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>ae40dda5ea89463f45f735360c348eef0</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Atom</type>
      <name>image</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a16ca1db99f14511a90d8dd7c03f76e6b</anchor>
      <arglist>(const Atom &amp;source, const Operator &amp;op) const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>image_neighbour</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>af51b2b4f0830cbc94a3a261df4c9c996</anchor>
      <arglist>(size_t atom_index, const Operator &amp;op)</arglist>
    </member>
    <member kind="function">
      <type>const Group &amp;</type>
      <name>group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>ac611e784f9154efe785c3e2ac001229a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable">
      <type>const Molecule &amp;</type>
      <name>m_molecule</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>afa6bc5ce0ed7f6f62a06b48cc70130da</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const Group &amp;</type>
      <name>m_group</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a68f942e5e77ab84b4f0ca8d8227c50d6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; std::vector&lt; size_t &gt; &gt;</type>
      <name>m_neighbours</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>ae04d437a6b951c7adce7834054abdedb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>CoordinateSystem::vec</type>
      <name>m_inertia_principal_values</name>
      <anchorfile>classmolpro_1_1point__charge__symmetry_1_1SymmetryMeasure.html</anchorfile>
      <anchor>a88a69f971c31e92cc71a63b8379a679c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>molpro</name>
    <filename>namespacemolpro.html</filename>
    <namespace>molpro::point_charge_symmetry</namespace>
  </compound>
  <compound kind="namespace">
    <name>molpro::point_charge_symmetry</name>
    <filename>namespacemolpro_1_1point__charge__symmetry.html</filename>
    <class kind="class">molpro::point_charge_symmetry::Operator</class>
    <class kind="class">molpro::point_charge_symmetry::Reflection</class>
    <class kind="class">molpro::point_charge_symmetry::Rotation</class>
    <class kind="class">molpro::point_charge_symmetry::Inversion</class>
    <class kind="class">molpro::point_charge_symmetry::Identity</class>
    <class kind="class">molpro::point_charge_symmetry::CoordinateSystem</class>
    <class kind="class">molpro::point_charge_symmetry::Group</class>
    <class kind="struct">molpro::point_charge_symmetry::Atom</class>
    <class kind="class">molpro::point_charge_symmetry::Molecule</class>
    <class kind="class">molpro::point_charge_symmetry::SymmetryMeasure</class>
    <class kind="class">molpro::point_charge_symmetry::Problem_optimise_frame</class>
    <class kind="class">molpro::point_charge_symmetry::Projector</class>
    <member kind="enumeration">
      <type></type>
      <name>RotationParameterType</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a98f75420a563d8c7ffdd3bc95ae992ea</anchor>
      <arglist></arglist>
      <enumvalue file="namespacemolpro_1_1point__charge__symmetry.html" anchor="a98f75420a563d8c7ffdd3bc95ae992eaace0be71e33226e4c1db2bcea5959f16b">Log</enumvalue>
      <enumvalue file="namespacemolpro_1_1point__charge__symmetry.html" anchor="a98f75420a563d8c7ffdd3bc95ae992eaa0a7532036415f2491bf5f952220827b8">Euler</enumvalue>
      <enumvalue file="namespacemolpro_1_1point__charge__symmetry.html" anchor="a98f75420a563d8c7ffdd3bc95ae992eaa3743af167c53361d795405561faac2b2">Quaternion</enumvalue>
      <enumvalue file="namespacemolpro_1_1point__charge__symmetry.html" anchor="a98f75420a563d8c7ffdd3bc95ae992eaa9b48a1f4556690f29fc77555c1a2bc24">TanEuler</enumvalue>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a0b01abd3eb56c3fda7a913a12f091271</anchor>
      <arglist>(std::ostream &amp;os, const Operator &amp;op)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>ae6a9ac05d70478a0214767358756539e</anchor>
      <arglist>(std::ostream &amp;os, const CoordinateSystem &amp;op)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>ab36b4b2d6bb3e75f623864de5be7fabd</anchor>
      <arglist>(std::ostream &amp;os, const CoordinateSystem::parameters_t &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>acc179a1598d43f50ae3c437cf4597869</anchor>
      <arglist>(std::ostream &amp;os, const Group &amp;g)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a9d8107f78e634186916a64e3751d82f7</anchor>
      <arglist>(std::ostream &amp;os, const Molecule &amp;op)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>cartesian_distance</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>aae65bd04ad416d22bc7c308cdcf2944f</anchor>
      <arglist>(const Molecule &amp;molecule1, const Molecule &amp;molecule2)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; double, std::vector&lt; double &gt; &gt;</type>
      <name>distance</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a8f28299831c41cc4e92243b4749c24d4</anchor>
      <arglist>(const Molecule &amp;molecule, const Molecule &amp;molecule0)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a98f18ad41b8d5d464a615c3910265489</anchor>
      <arglist>(std::ostream &amp;os, const SymmetryMeasure &amp;sm)</arglist>
    </member>
    <member kind="function">
      <type>Group</type>
      <name>discover_group</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a47866910d4d5ffb316585f43d0671b5c</anchor>
      <arglist>(const Molecule &amp;molecule, CoordinateSystem &amp;coordinate_system, double threshold=1e-10, int verbosity=-1)</arglist>
    </member>
    <member kind="function">
      <type>Group</type>
      <name>discover_group</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>ab153e85b11740778aef2c1133f80e2ea</anchor>
      <arglist>(const Molecule &amp;molecule, double threshold=1e-10, int verbosity=-1)</arglist>
    </member>
    <member kind="function">
      <type>Molecule</type>
      <name>molecule_localised</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a3e17936a15cca8f213f033c322053f09</anchor>
      <arglist>(const CoordinateSystem &amp;coordinate_system, const Molecule &amp;source)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>test_group</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>af050f224f2c9241d4cad07bce28757dc</anchor>
      <arglist>(const Molecule &amp;molecule, const Group &amp;group, double threshold=1e-6, int verbosity=-1)</arglist>
    </member>
    <member kind="function">
      <type>Eigen::Matrix3d</type>
      <name>find_axis_frame</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a79b6188e50aabbff3ff2c4dc996f52f4</anchor>
      <arglist>(const Molecule &amp;molecule, const Group &amp;group)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1point__charge__symmetry.html</anchorfile>
      <anchor>a59495758f9ca07cab0a7bdcce9cca7d7</anchor>
      <arglist>(std::ostream &amp;s, const std::vector&lt; T &gt; &amp;v)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Point Group Discovery and Repair</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md____w_point_charge_symmetry_point_charge_symmetry_README</docanchor>
  </compound>
</tagfile>
