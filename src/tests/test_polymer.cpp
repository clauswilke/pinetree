/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <string>
#include <vector>

#include "choices.hpp"
#include "feature.hpp"
#include "polymer.hpp"

namespace Helper {
bool termination_fired = false;
void ShiftMaskN(Polymer *polymer, int n) {
  for (int i = 0; i < n; i++) {
    polymer->ShiftMask();
  }
}
void EmitTermination(int index, const std::string &pol_name,
                     const std::string &last_gene) {
  termination_fired = true;
}
}

TEST_CASE("Polymer methods", "[Polymer]") {
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Promoter::Ptr prom;
  Terminator::Ptr term;
  prom = std::make_shared<Promoter>(Promoter("p1", 5, 15, interactions));
  std::map<std::string, double> efficiency;
  efficiency["ecolipol"] = 0.5;
  efficiency["rnapol"] = 0.8;
  term = std::make_shared<Terminator>(
      Terminator("t1", 50, 55, interactions, efficiency));

  std::vector<Element::Ptr> elements;
  elements.push_back(prom);
  elements.push_back(term);
  Mask mask = Mask("test_mask", 10, 100, interactions);

  Polymer polymer = Polymer("test_polymer", 1, 100, elements, mask);
  polymer.termination_signal_.Connect(Helper::EmitTermination);

  auto pol = std::make_shared<Polymerase>(Polymerase("ecolipol", 10, 30));

  SECTION("Construction") {
    REQUIRE(polymer.index() == 0);
    REQUIRE(term->IsCovered());
    REQUIRE(!term->WasCovered());
    REQUIRE(!term->WasUncovered());
    REQUIRE(prom->IsCovered());
    REQUIRE(!prom->WasCovered());
    REQUIRE(!prom->WasUncovered());
    REQUIRE(polymer.uncovered("p1") == 0);
    REQUIRE(polymer.uncovered("t1") == 0);
  }

  SECTION("Polymerase binding") {
    Random::seed(22);
    // Promoter should be covered an inaccessible
    REQUIRE_THROWS(polymer.Bind(pol, "p1"));
    Helper::ShiftMaskN(&polymer, 10);
    polymer.Bind(pol, "p1");
    REQUIRE(polymer.uncovered("p1") == 0);
    REQUIRE(polymer.prop_sum() == 30);
  }

  SECTION("Mask shifting") {
    polymer.ShiftMask();
    REQUIRE(prom->IsCovered());
    REQUIRE(term->IsCovered());
    // Shift mask 10 spaces
    Helper::ShiftMaskN(&polymer, 10);
    REQUIRE(!prom->IsCovered());
    REQUIRE(term->IsCovered());
    // Make sure mask can't shift past end of polymer
    Helper::ShiftMaskN(&polymer, 1000);
    REQUIRE(!term->IsCovered());
  }

  SECTION("Test termination") {
    Helper::ShiftMaskN(&polymer, 10);
    polymer.Bind(pol, "p1");
    // Record old propensity
    double old_prop_sum = polymer.prop_sum();
    polymer.Terminate(pol, "gene1");
    // Make sure propensity has changed and termination signal fired
    REQUIRE(polymer.prop_sum() != old_prop_sum);
    REQUIRE(Helper::termination_fired);
  }

  SECTION("Test coverings") {
    polymer.UncoverElement("p1");
    REQUIRE(polymer.uncovered("p1") == 1);
    polymer.CoverElement("p1");
    REQUIRE(polymer.uncovered("p1") == 0);
  }
}