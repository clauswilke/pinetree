/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <string>
#include <vector>

#include "feature.hpp"

TEST_CASE("Feature construction", "[Feature]")
{
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Feature feat = Feature(name, start, stop, interactions);

  // Check that interactions have been set up correctly
  REQUIRE(feat.CheckInteraction("ecolipol") == true);
  REQUIRE(feat.CheckInteraction("some other pol") == false);
}

TEST_CASE("Polymerase construction and movement", "[Polymerase]")
{
  std::string name = "testing!";
  int footprint = 40;
  int speed = 60;
  Polymerase pol = Polymerase(name, footprint, speed);

  REQUIRE(pol.start() == 0);
  REQUIRE(pol.stop() == 39);

  pol.Move();
  REQUIRE(pol.start() == 1);
  REQUIRE(pol.stop() == 40);

  pol.MoveBack();
  REQUIRE(pol.start() == 0);
  REQUIRE(pol.stop() == 39);
}

TEST_CASE("Mask construction and movement", "[Mask]")
{
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Mask mask = Mask(name, start, stop, interactions);

  // Check that interactions have been set up correctly
  REQUIRE(mask.CheckInteraction("ecolipol") == true);
  REQUIRE(mask.CheckInteraction("some other pol") == false);

  REQUIRE(mask.start() == 1);
  REQUIRE(mask.stop() == 10);

  mask.Recede();
  REQUIRE(mask.start() == 2);
  REQUIRE(mask.stop() == 10);
}

TEST_CASE("Element construction and state changes", "[Element]")
{
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Element elem = Element(name, start, stop, interactions);

  SECTION("Saving the state")
  {
    elem.Cover();
    elem.Cover();
    elem.Cover();
    REQUIRE(elem.IsCovered());
    REQUIRE(elem.WasCovered());

    elem.SaveState();
    REQUIRE(!elem.WasCovered());

    elem.Uncover();
    elem.Uncover();
    elem.Uncover();
    REQUIRE(!elem.IsCovered());
    REQUIRE(elem.WasUncovered());

    elem.SaveState();
    REQUIRE(!elem.WasUncovered());
  }

  SECTION("Uncovering and covering an element")
  {
    // Make sure covering counts don't go below zero
    elem.Uncover();
    elem.Uncover();
    elem.Uncover();
    REQUIRE(!elem.IsCovered());
    REQUIRE(!elem.IsCovered());

    elem.Cover();
    REQUIRE(elem.IsCovered());
    REQUIRE(elem.WasCovered());

    elem.Uncover();
    REQUIRE(!elem.IsCovered());
    REQUIRE(!elem.WasCovered());
    REQUIRE(!elem.WasUncovered());
  }
}

TEST_CASE("Terminator construction", "[Terminator]")
{
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Terminator term = Terminator(name, start, stop, interactions);
  REQUIRE(!term.readthrough());
  REQUIRE(term.reading_frame() == -1);
  // Test for reading frame checking
  term.set_reading_frame(1);
  REQUIRE(term.CheckInteraction("ecolipol", 1));
  REQUIRE(!term.CheckInteraction("ecolipol", 0));
  REQUIRE(!term.CheckInteraction("nopol", 1));
}
