/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <string>
#include <vector>

#include "polymer.hpp"
#include "feature.hpp"

TEST_CASE("Polymer construction", "[Polymer]")
{
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Promoter prom = Promoter("p1", 1, 10, interactions);
  Terminator term = Terminator("t1", 50, 55, interactions);

  std::vector<Element> elements;
  elements.push_back(prom);
  elements.push_back(term);
  Mask mask = Mask("test_mask", 20, 100, interactions);

  Polymer polymer = Polymer("test_polymer", 1, 100, elements, mask);

  REQUIRE(polymer.index() == 0);
  REQUIRE(polymer.uncovered("p1") == 1);
  REQUIRE(polymer.uncovered("t1") == 0);
}