/**
 * Modified based on the example cohesive_extrinsic.cc
 *
 * Modified by:
 *
 * Information for the original file:
 *
 * @file   cohesive_extrinsic.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Mon Jan 18 2016
 *
 * @brief  Test for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef std::map< std::pair<ID, ID>, ID > MaterialRules;
/* -------------------------------------------------------------------------- */
class MaterialRulesMaterialSelector : public MeshDataMaterialSelector<std::string>{
 public:
  MaterialRulesMaterialSelector(
    SolidMechanicsModelCohesive & model, const MaterialRules & rules,
    const ID & mesh_data_id) // what we have here is2 the name of model and also
                             // the name of different materials
    : MeshDataMaterialSelector<std::string>(mesh_data_id, model),
      rules(rules),
      mesh_data_id(mesh_data_id),
      model(model),
      mesh(model.getMesh()),
      mesh_facets(model.getMeshFacets()),
      spatial_dimension(model.getSpatialDimension()) {}

/* -------------------------------------------------------------------------- */
UInt operator()(const Element & element) {
  if (mesh_facets.getSpatialDimension(element.type) ==
      (spatial_dimension - 1)) {
    const std::vector<Element> & element_to_subelement =
        mesh_facets.getElementToSubelement(element.type,
                                           element.ghost_type)(element.element);
    const Element & el1 = element_to_subelement[0];
    const Element & el2 = element_to_subelement[1];

    ID id1 = mesh.getData<std::string>(mesh_data_id, el1.type,
                                       el1.ghost_type)(el1.element);

    ID id2 = id1;
    if (el2 != ElementNull) {
      id2 = mesh.getData<std::string>(mesh_data_id, el2.type,
                                      el2.ghost_type)(el2.element);
    } else {
      std::cout << "Trying to set material\"("
		<< id1 << ", " << id2
		<< ")\" to a boundary facet" << std::endl;
      Array<bool> & facets_check =
	model.getElementInserter().getCheckFacets(element.type, element.ghost_type);
      facets_check(element.element) = false;
    }

    MaterialRules::const_iterator rit = rules.find(std::make_pair(id1, id2));
    if (rit == rules.end()) {
      rit = rules.find(std::make_pair(id2, id1));
    }

    if (rit != rules.end()) {
      return model.getMaterialIndex(rit->second);
    } else {
      DefaultMaterialCohesiveSelector cohesive_selector(model);
      return cohesive_selector(element);
    }
  } else if (element.kind() == _ek_cohesive) {
    DefaultMaterialCohesiveSelector cohesive_selector(model);
    return cohesive_selector(element);
  } else {
    return MeshDataMaterialSelector::operator()(element);
  }
}
  private:
  SolidMechanicsModelCohesive & model;
  ID mesh_data_id;
  const Mesh & mesh;
  const Mesh & mesh_facets;
  UInt spatial_dimension;
  MaterialRules rules;
};
/* -------------------------------------------------------------------------- */

class InternalFieldSetter {
public:
  InternalFieldSetter(SolidMechanicsModel &model, const ID &material)
      : model(model), material(model.getMaterial(material)),
        quadrature_points("quadrature_points:" + material,
                          "internal_field_setter") {
    quadrature_points.initialize(model.getFEEngine(), _all_ghost_types = true,
                                 _element_kind = _ek_not_defined,
                                 _nb_component = model.getSpatialDimension());
    model.getFEEngine().interpolateOnIntegrationPoints(
        model.getMesh().getNodes(), quadrature_points,
        &this->material.getElementFilter());
  }

  template <class Func, class... Ns>
  void setField(const ID &field_id, Func &&func, Ns... ns) {
    auto &&field = this->material.getInternal<Real>(field_id);
    auto &&fe_engine = model.getFEEngine();

    IntegrationPoint qp;
    for (auto ghost_type : ghost_types) {
      qp.ghost_type = ghost_type;
      for (auto type :
           field.elementTypes(_all_dimensions, ghost_type, _ek_not_defined)) {
        qp.type = type;
        auto nb_quadrature_points = fe_engine.getNbIntegrationPoints(type);
        auto &&field_array = field(type, ghost_type);

        for (auto &&data : zip(arange(field_array.size()),
                               make_view(quadrature_points(type, ghost_type),
                                         model.getSpatialDimension()),
                               make_view(field_array, ns...))) {
          qp.element = std::get<0>(data) / nb_quadrature_points;
          qp.num_point = std::get<0>(data) % nb_quadrature_points;

          func(qp, std::get<1>(data), std::get<2>(data));
        }
      }
    }
  }

private:
  SolidMechanicsModel &model;
  Material &material;

  ElementTypeMapArray<Real> quadrature_points;
};

int main(int argc, char *argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 1000;

  Mesh mesh(spatial_dimension);
  mesh.read("triangle_2.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  MaterialRules rules;
  rules[std::make_pair("cathode", "cathode")] = "cathode_cathode";
  rules[std::make_pair("ceramic", "ceramic")] = "ceramic_ceramic";
  rules[std::make_pair("anode", "anode")] = "anode_anode";
  rules[std::make_pair("cathode", "ceramic")] = "cathode_ceramic";
  rules[std::make_pair("anode", "ceramic")] = "anode_ceramic";
  auto material_selector = std::make_shared<MaterialRulesMaterialSelector>(
	model, rules, "physical_names");
   model.setMaterialSelector(material_selector);
   
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = true);

  InternalFieldSetter cathode_setter(model, "cathode");
  InternalFieldSetter ceramic_setter(model, "ceramic");
  InternalFieldSetter anode_setter(model, "anode");

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  CohesiveElementInserter &inserter = model.getElementInserter();
  inserter.setLimit(_y, 0.0, 0.50);
  model.updateAutomaticInsertion();

  Array<Real> &position = mesh.getNodes();
  Array<Real> &velocity = model.getVelocity();
  Array<bool> &boundary = model.getBlockedDOFs();
  Array<Real> &displacement = model.getDisplacement();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    // if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
      if (position(n, 1) < 0.01)
      boundary(n, 1) = true;

      //if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      if (position(n, 0) < 0.01)
      boundary(n, 0) = true;
  }

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("delta_T");
  model.dump();

  /// initial conditions
  //Real loading_rate = 0.5;
  //  Real disp_update = loading_rate * time_step;
  // for (UInt n = 0; n < nb_nodes; ++n) {
  //  velocity(n, 1) = loading_rate * position(n, 1);
  //}

  Vector<Real> a{-1., -1.};
  Vector<Real> b{1., 1.};
  auto direction = (b - a).normalize();
  auto L = (b - a).norm();

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    // /// update displacement on extreme nodes
    // for (UInt n = 0; n < nb_nodes; ++n) {
    //   if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
    //     displacement(n, 1) += disp_update * position(n, 1);
    // }

    // Here I would like to create a temperature field which linearly changes
    // from point (-1,-1) to point (1,1), the change is then gradually
    // increasing to cause some thermal cracks.

    Real scale = 0.1 * s;
    cathode_setter.setField(
        "delta_T", [&](IntegrationPoint qp, Vector<Real> &pos, Real &delta_T) {
          delta_T = pos.dot(direction) / L * scale;
        });
    ceramic_setter.setField(
        "delta_T", [&](IntegrationPoint qp, Vector<Real> &pos, Real &delta_T) {
          delta_T = pos.dot(direction) / L * scale;
        });
    anode_setter.setField(
        "delta_T", [&](IntegrationPoint qp, Vector<Real> &pos, Real &delta_T) {
          delta_T = pos.dot(direction) / L * scale;
        });   
    model.checkCohesiveStress();
    model.solveStep();

    if (s % 10 == 0) {
      model.dump(time_step * s, s);

      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }
  }

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 200 * std::sqrt(2);
  std::cout << Ed << " " << Edt << std::endl;
  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  return EXIT_SUCCESS;
}
