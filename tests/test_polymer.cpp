/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <iostream>
#include <string>
#include <vector>

#include "choices.hpp"
#include "feature.hpp"
#include "polymer.hpp"

namespace Helper {
bool termination_fired = false;
bool transcript_fired = false;
Transcript::Ptr my_transcript;
void ShiftMaskN(Polymer::Ptr polymer, int n) {
  for (int i = 0; i < n; i++) {
    polymer->ShiftMask();
  }
}
void MovePolymeraseN(Polymer::Ptr polymer, Polymerase::Ptr pol, int n) {
  int index = polymer->PolymeraseIndex(pol);
  for (int i = 0; i < n; i++) {
    polymer->Move(index);
  }
}
void EmitTermination(int index, const std::string &pol_name,
                     const std::string &last_gene) {
  termination_fired = true;
}
void EmitTranscript(const Transcript::Ptr &transcript) {
  transcript_fired = true;
  transcript->InitElements();
  my_transcript = transcript;
}
}

TEST_CASE("Polymer methods", "[Polymer]") {
  std::map<std::string, double> interactions {
    {"ecolipol", 1.0}, 
    {"rnapol", 1.0},
  };
  Promoter::Ptr prom;
  Terminator::Ptr term;
  prom = std::make_shared<Promoter>("p1", 5, 15, interactions);
  std::map<std::string, double> efficiency;
  efficiency["ecolipol"] = 0.6;
  efficiency["rnapol"] = 1.0;
  term = std::make_shared<Terminator>("t1", 50, 55, efficiency);

  std::vector<Element::Ptr> elements;
  elements.push_back(prom);
  elements.push_back(term);
  std::map<std::string, double> mask_interactions {
    {"ecolipol", 1.0}, 
  };
  Mask mask = Mask("test_mask", 10, 100, mask_interactions);

  auto polymer =
      std::make_shared<Polymer>("test_polymer", 1, 100, elements, mask);
  polymer->termination_signal_.Connect(Helper::EmitTermination);
  polymer->InitElements();

  auto pol = std::make_shared<Polymerase>("ecolipol", 10, 30);
  auto pol2 = std::make_shared<Polymerase>("rnapol", 10, 30);
  auto pol3 = std::make_shared<Polymerase>("rnapol", 10, 30);

  SECTION("Construction") {
    REQUIRE(polymer->index() == 0);
    REQUIRE(term->IsCovered());
    REQUIRE(!term->WasCovered());
    REQUIRE(!term->WasUncovered());
    REQUIRE(prom->IsCovered());
    REQUIRE(!prom->WasCovered());
    REQUIRE(!prom->WasUncovered());
    REQUIRE(polymer->uncovered("p1") == 0);
    REQUIRE(polymer->uncovered("t1") == 0);
  }

  SECTION("Polymerase binding") {
    Random::seed(22);
    // Promoter should be covered an inaccessible
    REQUIRE_THROWS(polymer->Bind(pol, "p1"));
    Helper::ShiftMaskN(polymer, 10);
    polymer->Bind(pol, "p1");
    REQUIRE(polymer->uncovered("p1") == 0);
    REQUIRE(polymer->prop_sum() == 30);
  }

  SECTION("Mask shifting") {
    polymer->ShiftMask();
    REQUIRE(prom->IsCovered());
    REQUIRE(term->IsCovered());
    // Shift mask 10 spaces
    Helper::ShiftMaskN(polymer, 10);
    REQUIRE(!prom->IsCovered());
    REQUIRE(term->IsCovered());
    // Make sure mask can't shift past end of polymer
    Helper::ShiftMaskN(polymer, 1000);
    REQUIRE(!term->IsCovered());
  }

  SECTION("Test termination") {
    Helper::ShiftMaskN(polymer, 10);
    polymer->Bind(pol, "p1");
    // Record old propensity
    double old_prop_sum = polymer->prop_sum();
    polymer->Terminate(pol, "gene1");
    // Make sure propensity has changed and termination signal fired
    REQUIRE(polymer->prop_sum() != old_prop_sum);
    REQUIRE(Helper::termination_fired);
  }

  SECTION("Test coverings") {
    polymer->UncoverElement("p1");
    REQUIRE(polymer->uncovered("p1") == 1);
    polymer->CoverElement("p1");
    REQUIRE(polymer->uncovered("p1") == 0);
  }

  SECTION("Moving polymerase") {
    // Shift mask back to expose promoter
    Helper::ShiftMaskN(polymer, 11);
    polymer->Bind(pol, "p1");
    // Make sure promoter is recorded as covered
    REQUIRE(polymer->uncovered("p1") == 0);
    REQUIRE(prom->IsCovered());
    // Move polymerase 10 spaces and re-expose promoter
    Helper::MovePolymeraseN(polymer, pol, 11);
    REQUIRE(polymer->uncovered("p1") == 1);
    REQUIRE(!prom->IsCovered());
    // Check for collisions between polymerases
    polymer->Bind(pol2, "p1");
    Helper::MovePolymeraseN(polymer, pol2, 15);
    REQUIRE(pol2->stop() + 1 == pol->start());
    // Remove pol and make sure pol2 can't shift mask'
    polymer->Terminate(pol, "gene1");
    int old_stop = pol2->stop();
    Helper::MovePolymeraseN(polymer, pol2, 20);
    // Pol should only be able to move 10 spaces
    REQUIRE(pol2->stop() == old_stop + 10);
    // Expose terminator
    REQUIRE(polymer->uncovered("t1") == 0);
    Helper::ShiftMaskN(polymer, 60);
    REQUIRE(polymer->uncovered("t1") == 1);
    // Test termination, pol should detach as soon as it hits terminator
    Helper::MovePolymeraseN(polymer, pol2, 25);
    REQUIRE_THROWS(polymer->Move(polymer->PolymeraseIndex(pol2)));

    // Now check for readthrough
    polymer->Bind(pol, "p1");
    Random::seed(55);
    Helper::MovePolymeraseN(polymer, pol, 42);
    REQUIRE(term->readthrough());
    Helper::MovePolymeraseN(polymer, pol, 10);
    // Now run polymerase off of the end of the transcript
    Helper::MovePolymeraseN(polymer, pol, 35);
    REQUIRE_THROWS(polymer->Move(polymer->PolymeraseIndex(pol)));
    // Run polymerase through terminator again, this time it should terminate
    Random::seed(19);
    polymer->Bind(pol, "p1");
    Helper::MovePolymeraseN(polymer, pol, 36);
    REQUIRE(!term->readthrough());
    REQUIRE_THROWS(polymer->Move(polymer->PolymeraseIndex(pol)));
  }
}

TEST_CASE("Polymer methods with multipromoter", "[Polymer]") {
  std::map<std::string, double> interactions {
    {"ecolipol", 1.0}
  };
  Promoter::Ptr prom1, prom2, prom3;
  Terminator::Ptr term;
  prom1 = std::make_shared<Promoter>("p1", 5, 15, interactions);
  prom2 = std::make_shared<Promoter>("p2", 16, 20, interactions);
  prom3 = std::make_shared<Promoter>("p3", 21, 30, interactions);
  std::map<std::string, double> efficiency;
  efficiency["ecolipol"] = 1.0;
  term = std::make_shared<Terminator>("t1", 31, 33, efficiency);

  std::vector<Element::Ptr> elements;
  elements.push_back(prom1);
  elements.push_back(prom2);
  elements.push_back(prom3);
  elements.push_back(term);
  std::map<std::string, double> mask_interactions {
    {"ecolipol", 1.0}, 
  };
  Mask mask = Mask("test_mask", 100, 100, mask_interactions);

  auto polymer =
      std::make_shared<Polymer>("test_polymer", 1, 100, elements, mask);
  polymer->termination_signal_.Connect(Helper::EmitTermination);
  polymer->InitElements();

  auto pol = std::make_shared<Polymerase>("ecolipol", 10, 30);
  auto pol2 = std::make_shared<Polymerase>("rnapol", 10, 30);
  auto pol3 = std::make_shared<Polymerase>("rnapol", 10, 30);

  SECTION("Moving polymerase") {
    polymer->Bind(pol, "p1");
    Helper::MovePolymeraseN(polymer, pol, 7);
    REQUIRE(prom1->IsCovered());
    REQUIRE(prom2->IsCovered());
    REQUIRE(prom3->IsCovered());

    Helper::MovePolymeraseN(polymer, pol, 4);
    REQUIRE(!prom1->IsCovered());
    REQUIRE(prom2->IsCovered());
    REQUIRE(prom3->IsCovered());
  }
}

TEST_CASE("Genome methods", "[Polymer]") {
  // Set up transcript template
  Element::VecPtr transcript_template;
  std::map<std::string, double> interactions {
    {"ribosome", 1.0}, 
  };
  std::map<std::string, double> efficiency;
  efficiency["ribosome"] = 1.0;
  transcript_template.push_back(
      std::make_shared<Promoter>("rbs", 0, 10, interactions));
  transcript_template.push_back(
      std::make_shared<Terminator>("stop", 209, 210, efficiency));
  transcript_template.push_back(
      std::make_shared<Promoter>("rbs", 215, 230, interactions));
  transcript_template.push_back(
      std::make_shared<Terminator>("stop", 269, 270, efficiency));
  transcript_template.push_back(
      std::make_shared<Promoter>("rbs", 285, 300, interactions));
  transcript_template.push_back(
      std::make_shared<Terminator>("stop", 599, 600, efficiency));
  // Set up genome
  interactions = {
    {"ecolipol", 1.0}, 
    {"rnapol", 1.0},
  };
  Promoter::Ptr prom;
  Terminator::Ptr term;
  prom = std::make_shared<Promoter>("p1", 5, 15, interactions);
  efficiency["ecolipol"] = 1.0;
  term = std::make_shared<Terminator>("t1", 601, 605, efficiency);
  Element::VecPtr elements;
  elements.push_back(prom);
  elements.push_back(term);
 std::map<std::string, double> mask_interactions {
    {"ecolipol", 1.0}, 
  };
  Mask mask = Mask("test_mask", 10, 700, mask_interactions);
  auto genome = std::make_shared<Genome>("test_genome", 700, elements,
                                         transcript_template, mask);
  genome->termination_signal_.Connect(Helper::EmitTermination);
  genome->transcript_signal_.Connect(Helper::EmitTranscript);

  genome->InitElements();

  auto pol = std::make_shared<Polymerase>("ecolipol", 10, 30);

  SECTION("Bind polymerase to genome") {
    REQUIRE(!Helper::transcript_fired);
    Helper::ShiftMaskN(genome, 20);
    genome->Bind(pol, "p1");
    REQUIRE(Helper::transcript_fired);
    REQUIRE(Helper::my_transcript->stop() == genome->stop());
  }

  SECTION("Polymerase movement and transcript unmasking") {
    Helper::ShiftMaskN(genome, 20);
    genome->Bind(pol, "p1");
    Helper::MovePolymeraseN(genome, pol, 500);
    REQUIRE(Helper::my_transcript->uncovered("stop") == 2);
    REQUIRE(Helper::my_transcript->uncovered("rbs") == 2);
    Helper::MovePolymeraseN(genome, pol, 87);
    REQUIRE(Helper::my_transcript->uncovered("stop") == 3);
  }
}

TEST_CASE("Variable translation rates", "[Polymer]") {
  std::map<std::string, double> interactions {
    {"ecolipol", 1.0}, 
    {"rnapol", 1.0},
  };
  Promoter::Ptr prom;
  Terminator::Ptr term;
  prom = std::make_shared<Promoter>("p1", 5, 15, interactions);
  std::map<std::string, double> efficiency;
  efficiency["ecolipol"] = 0.6;
  efficiency["rnapol"] = 1.0;
  term = std::make_shared<Terminator>("t1", 50, 55, efficiency);

  std::vector<Element::Ptr> elements;
  elements.push_back(prom);
  elements.push_back(term);
  std::map<std::string, double> mask_interactions {
    {"ecolipol", 1.0}, 
  };
  Mask mask = Mask("test_mask", 10, 100, mask_interactions);

  std::vector<double> weights(100, 0.1);
  for (int i = 0; i < 50; i++) {
    weights[i] = 2.0;
  }
  auto polymer = std::make_shared<Polymer>("test_polymer", 1, 100, elements,
                                           mask, weights);
  polymer->termination_signal_.Connect(Helper::EmitTermination);
  polymer->InitElements();

  auto pol = std::make_shared<Polymerase>("ecolipol", 10, 30);
  auto pol2 = std::make_shared<Polymerase>("rnapol", 10, 30);
  auto pol3 = std::make_shared<Polymerase>("rnapol", 10, 30);

  SECTION("Execution") {
    Random::seed(22);
    // Shift mask back to expose promoter
    Helper::ShiftMaskN(polymer, 10);
    polymer->Bind(pol, "p1");
    // Make sure promoter is recorded as covered
    REQUIRE(polymer->uncovered("p1") == 0);
    REQUIRE(prom->IsCovered());
    // Move polymerase 10 spaces and re-expose promoter
    Helper::MovePolymeraseN(polymer, pol, 11);
    REQUIRE(polymer->uncovered("p1") == 1);
    REQUIRE(!prom->IsCovered());
    // Check for collisions between polymerases
    polymer->Bind(pol2, "p1");
    for (int i = 0; i < 100; i++) {
      polymer->Execute();
    }
  }
}