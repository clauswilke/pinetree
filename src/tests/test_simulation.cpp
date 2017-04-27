#include <catch.hpp>

#include "simulation.hpp"

TEST_CASE("SpeciesReaction methods", "[Reaction]") {
  auto reaction = std::make_shared<SpeciesReaction>(
      1000, std::vector<std::string>{"reactant1", "reactant2"},
      std::vector<std::string>{"product1", "product2"});
  auto &tracker = SpeciesTracker::Instance();
  tracker.Clear();
  tracker.Register(reaction);
  tracker.Increment("reactant1", 2);
  tracker.Increment("reactant2", 3);

  SECTION("Initialization and registration") {
    REQUIRE_THROWS(SpeciesReaction(
        1000, std::vector<std::string>{"reactant1", "reactant2", "reactant3"},
        std::vector<std::string>{"product1", "product2"}));
    // Make sure that maps have been set up appropriately
    auto reactions = tracker.FindReactions("reactant1");
    REQUIRE(std::find(reactions.begin(), reactions.end(), reaction) !=
            reactions.end());
    reactions = tracker.FindReactions("reactant2");
    REQUIRE(std::find(reactions.begin(), reactions.end(), reaction) !=
            reactions.end());
    reactions = tracker.FindReactions("product1");
    REQUIRE(std::find(reactions.begin(), reactions.end(), reaction) !=
            reactions.end());
    reactions = tracker.FindReactions("product2");
    REQUIRE(std::find(reactions.begin(), reactions.end(), reaction) !=
            reactions.end());
  }

  SECTION("Propensity calculation") {
    REQUIRE((1000 * 2 * 3) / (double(6.0221409e+23) * double(8e-15)) ==
            reaction->CalculatePropensity());
  }

  SECTION("Execution") {
    reaction->Execute();
    REQUIRE(tracker.species("reactant1") == 1);
    REQUIRE(tracker.species("reactant2") == 2);
    REQUIRE(tracker.species("product1") == 1);
    REQUIRE(tracker.species("product1") == 1);
  }
}

TEST_CASE("SpeciesTracker methods", "[SpeciesTracker]") {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Clear();

  SECTION("Increment species") {
    tracker.Increment("reactant1", 1);
    REQUIRE(tracker.species("reactant1") == 1);
    tracker.Increment("reactant2", 1);
    REQUIRE(tracker.species("reactant2") == 1);
    tracker.Increment("reactant1", 2);
    REQUIRE(tracker.species("reactant1") == 3);
  }

  SECTION("Add polymer") {
    Polymer::Ptr polymer;
    tracker.Add("promoter1", polymer);
    auto polymers = tracker.FindPolymers("promoter1");
    REQUIRE(std::find(polymers.begin(), polymers.end(), polymer) !=
            polymers.end());
    tracker.Add("promoter2", polymer);
    polymers = tracker.FindPolymers("promoter2");
    REQUIRE(std::find(polymers.begin(), polymers.end(), polymer) !=
            polymers.end());
  }
}