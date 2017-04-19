/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <string>
#include <vector>

#include "polymer.hpp"
#include "feature.hpp"
#include "choices.hpp"

namespace Helper
{
void ShiftMaskN(Polymer *polymer, int n)
{
  for (int i = 0; i < n; i++)
  {
    polymer->ShiftMask();
  }
}
}

TEST_CASE("Polymer methods", "[Polymer]")
{
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Promoter::Ptr prom;
  Terminator::Ptr term;
  prom = std::make_shared<Promoter>(Promoter("p1", 5, 15, interactions));
  term = std::make_shared<Terminator>(Terminator("t1", 50, 55, interactions));

  std::vector<Element::Ptr> elements;
  elements.push_back(prom);
  elements.push_back(term);
  Mask mask = Mask("test_mask", 10, 100, interactions);

  Polymer polymer = Polymer("test_polymer", 1, 100, elements, mask);

  SECTION("Construction")
  {
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

  SECTION("Polymerase binding")
  {
    Random::seed(22);
    Helper::ShiftMaskN(&polymer, 10);
    auto pol = std::make_shared<Polymerase>(Polymerase("ecolipol", 10, 30));
    polymer.Bind(pol, "p1");
    REQUIRE(polymer.uncovered("p1") == 0);
    REQUIRE(polymer.prop_sum() == 30);
  }

  SECTION("Mask shifting")
  {
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
}