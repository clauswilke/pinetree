/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <string>
#include <vector>

#include "feature.hpp"

TEST_CASE("Feature construction", "[Feature]") {
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Feature feat = Feature(name, start, stop, interactions);

  // Check that interactions have been set up correctly
  REQUIRE(feat.check_interaction("ecolipol") == true);
  REQUIRE(feat.check_interaction("some other pol") == false);
}

TEST_CASE("Polymerase construction and movement", "[Polymerase]") {
  std::string name = "testing!";
  int footprint = 40;
  int speed = 60;
  Polymerase pol = Polymerase(name, footprint, speed);

  REQUIRE(pol.get_start() == 0);
  REQUIRE(pol.get_stop() == 39);

  pol.move();
  REQUIRE(pol.get_start() == 1);
  REQUIRE(pol.get_stop() == 40);

  pol.move_back();
  REQUIRE(pol.get_start() == 0);
  REQUIRE(pol.get_stop() == 39);
}

TEST_CASE("Mask construction and movement", "[Mask]") {
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Mask mask = Mask(name, start, stop, interactions);

  // Check that interactions have been set up correctly
  REQUIRE(mask.check_interaction("ecolipol") == true);
  REQUIRE(mask.check_interaction("some other pol") == false);

  REQUIRE(mask.get_start() == 1);
  REQUIRE(mask.get_stop() == 10);

  mask.recede();
  REQUIRE(mask.get_start() == 2);
  REQUIRE(mask.get_stop() == 10);
}

TEST_CASE("Element construction and state changes", "[Element]") {
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Element elem = Element(name, start, stop, interactions);

  SECTION("Saving the state") {
    elem.cover();
    elem.cover();
    elem.cover();
    REQUIRE(elem.is_covered());
    REQUIRE(elem.was_covered());

    elem.save_state();
    REQUIRE(!elem.was_covered());

    elem.uncover();
    elem.uncover();
    elem.uncover();
    REQUIRE(!elem.is_covered());
    REQUIRE(elem.was_uncovered());

    elem.save_state();
    REQUIRE(!elem.was_uncovered());
  }

  SECTION("Uncovering and covering an element") {
    // Make sure covering counts don't go below zero
    elem.uncover();
    elem.uncover();
    elem.uncover();
    REQUIRE(!elem.is_covered());
    REQUIRE(!elem.is_covered());

    elem.cover();
    REQUIRE(elem.is_covered());
    REQUIRE(elem.was_covered());

    elem.uncover();
    REQUIRE(!elem.is_covered());
    REQUIRE(!elem.was_covered());
    REQUIRE(!elem.was_uncovered());
  }
}
