/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <string>
#include <vector>

#include "polymer.hpp"
#include "feature.hpp"

TEST_CASE("Polymer construction", "[Polymer]")
{
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Element feat = Element(name, start, stop, interactions);

  std::vector<Element> elements;
  elements.push_back(feat);
  Mask mask = Mask("test_mask", 50, 100, interactions);

  Polymer polymer = Polymer("test_polymer", 1, 100, elements, mask);

  REQUIRE(polymer.get_index() == 0);
}