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
void MovePolymeraseN(Polymer *polymer, Polymerase::Ptr pol, int n) {
  for (int i = 0; i < n; i++) {
    polymer->Move(pol);
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
  efficiency["ecolipol"] = 0.6;
  efficiency["rnapol"] = 1.0;
  term = std::make_shared<Terminator>(
      Terminator("t1", 50, 55, interactions, efficiency));

  std::vector<Element::Ptr> elements;
  elements.push_back(prom);
  elements.push_back(term);
  std::vector<std::string> mask_interactions = {"ecolipol"};
  Mask mask = Mask("test_mask", 10, 100, mask_interactions);

  Polymer polymer = Polymer("test_polymer", 1, 100, elements, mask);
  polymer.termination_signal_.Connect(Helper::EmitTermination);

  auto pol = std::make_shared<Polymerase>(Polymerase("ecolipol", 10, 30));
  auto pol2 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 30));
  auto pol3 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 30));

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

  SECTION("Moving polymerase") {
    // Shift mask back to expose promoter
    Helper::ShiftMaskN(&polymer, 10);
    polymer.Bind(pol, "p1");
    // Make sure promoter is recorded as covered
    REQUIRE(polymer.uncovered("p1") == 0);
    REQUIRE(prom->IsCovered());
    // Move polymerase 10 spaces and re-expose promoter
    Helper::MovePolymeraseN(&polymer, pol, 11);
    REQUIRE(polymer.uncovered("p1") == 1);
    REQUIRE(!prom->IsCovered());
    // Check for collisions between polymerases
    polymer.Bind(pol2, "p1");
    Helper::MovePolymeraseN(&polymer, pol2, 15);
    REQUIRE(pol2->stop() + 1 == pol->start());
    // Remove pol and make sure pol2 can't shift mask'
    polymer.Terminate(pol, "gene1");
    int old_stop = pol2->stop();
    Helper::MovePolymeraseN(&polymer, pol2, 20);
    // Pol should only be able to move 10 spaces
    REQUIRE(pol2->stop() == old_stop + 10);
    // Expose terminator
    REQUIRE(polymer.uncovered("t1") == 0);
    Helper::ShiftMaskN(&polymer, 60);
    REQUIRE(polymer.uncovered("t1") == 1);
    // Test termination, pol should detach as soon as it hits terminator
    Helper::MovePolymeraseN(&polymer, pol2, 25);
    REQUIRE_THROWS(polymer.Move(pol2));

    // Now check for readthrough
    polymer.Bind(pol, "p1");
    Random::seed(55);
    Helper::MovePolymeraseN(&polymer, pol, 42);
    REQUIRE(term->readthrough());
    Helper::MovePolymeraseN(&polymer, pol, 10);
    // Now run polymerase off of the end of the transcript
    Helper::MovePolymeraseN(&polymer, pol, 35);
    REQUIRE_THROWS(polymer.Move(pol));
    // Run polymerase through terminator again, this time it should terminate
    Random::seed(19);
    polymer.Bind(pol, "p1");
    Helper::MovePolymeraseN(&polymer, pol, 36);
    REQUIRE(!term->readthrough());
    REQUIRE_THROWS(polymer.Move(pol));
  }
}

TEST_CASE("Polymer methods with multipromoter", "[Polymer]") {
  std::vector<std::string> interactions = {"ecolipol"};
  Promoter::Ptr prom1, prom2, prom3;
  Terminator::Ptr term;
  prom1 = std::make_shared<Promoter>(Promoter("p1", 5, 15, interactions));
  prom2 = std::make_shared<Promoter>(Promoter("p2", 16, 20, interactions));
  prom3 = std::make_shared<Promoter>(Promoter("p3", 21, 30, interactions));
  std::map<std::string, double> efficiency;
  efficiency["ecolipol"] = 1.0;
  term = std::make_shared<Terminator>(
      Terminator("t1", 31, 33, interactions, efficiency));

  std::vector<Element::Ptr> elements;
  elements.push_back(prom1);
  elements.push_back(prom2);
  elements.push_back(prom3);
  elements.push_back(term);
  std::vector<std::string> mask_interactions = {"ecolipol"};
  Mask mask = Mask("test_mask", 100, 100, mask_interactions);

  Polymer polymer = Polymer("test_polymer", 1, 100, elements, mask);
  polymer.termination_signal_.Connect(Helper::EmitTermination);

  auto pol = std::make_shared<Polymerase>(Polymerase("ecolipol", 10, 30));
  auto pol2 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 30));
  auto pol3 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 30));

  SECTION("Moving polymerase") {
    polymer.Bind(pol, "p1");
    Helper::MovePolymeraseN(&polymer, pol, 7);
    REQUIRE(prom1->IsCovered());
    REQUIRE(prom2->IsCovered());
    REQUIRE(prom3->IsCovered());

    Helper::MovePolymeraseN(&polymer, pol, 4);
    REQUIRE(!prom1->IsCovered());
    REQUIRE(prom2->IsCovered());
    REQUIRE(prom3->IsCovered());
  }
}

TEST_CASE("Genome methods", "[Polymer]") {
  // Set up transcript template
  std::vector<Element::Ptr> transcript_template;
  std::vector<std::string> interactions = {"ribosome"};
  std::map<std::string, double> efficiency;
  efficiency["ribosome"] = 1.0;
  Promoter::Ptr rbs1, rbs2, rbs3;
  Terminator::Ptr stop1, stop2, stop3;
  rbs1 = std::make_shared<Promoter>(Promoter("rbs", 0, 10, interactions));
  transcript_template.push_back(rbs1);
  stop1 = std::make_shared<Terminator>(
      Terminator("stop", 209, 210, interactions, efficiency));
  transcript_template.push_back(stop1);
  rbs2 = std::make_shared<Promoter>(Promoter("rbs", 215, 230, interactions));
  transcript_template.push_back(rbs2);
  stop2 = std::make_shared<Terminator>(
      Terminator("stop", 269, 270, interactions, efficiency));
  transcript_template.push_back(stop2);
  rbs3 = std::make_shared<Promoter>(Promoter("rbs", 285, 300, interactions));
  transcript_template.push_back(rbs3);
  stop3 = std::make_shared<Terminator>(
      Terminator("stop", 599, 600, interactions, efficiency));
  transcript_template.push_back(stop3);
  // Set up genome
  interactions = {"ecolipol", "rnapol"};
  Promoter::Ptr prom;
  Terminator::Ptr term;
  prom = std::make_shared<Promoter>(Promoter("p1", 5, 15, interactions));
  efficiency["ecolipol"] = 1.0;
  term = std::make_shared<Terminator>(
      Terminator("t1", 601, 605, interactions, efficiency));
  std::vector<Element::Ptr> elements;
  elements.push_back(prom);
  elements.push_back(term);
  std::vector<std::string> mask_interactions = {"ecolipol"};
  Mask mask = Mask("test_mask", 10, 100, mask_interactions);
  Genome genome =
      Genome("test_genome", 700, elements, transcript_template, mask);
  genome.termination_signal_.Connect(Helper::EmitTermination);

  auto pol = std::make_shared<Polymerase>(Polymerase("ecolipol", 10, 30));

  SECTION("Bind polymerase to genome") {
    genome.Bind(pol, "p1");
    // tests here
  }
}