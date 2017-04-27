#include <catch.hpp>

#include "simulation.hpp"

TEST_CASE("SpeciesReaction methods", "[Reaction]") {
  auto reaction = std::make_shared<SpeciesReaction>(
      1000, std::vector<std::string>{"reactant1", "reactant2"},
      std::vector<std::string>{"product1", "product2"});
  SpeciesTracker::Register(reaction);

  SECTION("Initialization") {
    REQUIRE_THROWS(SpeciesReaction(
        1000, std::vector<std::string>{"reactant1", "reactant2", "reactant3"},
        std::vector<std::string>{"product1", "product2"}));
  }
}