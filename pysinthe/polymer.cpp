#include "polymer.hpp"
#include "choices.hpp"
#include "tracker.hpp"

#include <iostream>

Polymer::Polymer(const std::string &name, int start, int stop,
                 const Element::VecPtr &elements, const Mask &mask)
    : index_(0), name_(name), start_(start), stop_(stop), elements_(elements),
      mask_(mask), prop_sum_(0) {
  weights_ = std::vector<double>(stop - start + 2, 1.0);
}

Polymer::Polymer(const std::string &name, int start, int stop,
                 const Element::VecPtr &elements, const Mask &mask,
                 const std::vector<double> &weights)
    : index_(0), name_(name), start_(start), stop_(stop), elements_(elements),
      mask_(mask), weights_(weights), prop_sum_(0) {
  if (weights_.size() != (stop_ - start_ + 2)) {
    throw std::length_error("Weights vector is not the correct size. " +
                            std::to_string(weights.size()) + " " +
                            std::to_string(stop_ - start_ + 2));
  }
}

void Polymer::InitElements() {
  for (auto &element : elements_) {
    element->cover_signal_.ConnectMember(shared_from_this(),
                                         &Polymer::CoverElement);
    element->uncover_signal_.ConnectMember(shared_from_this(),
                                           &Polymer::UncoverElement);
    // Cover masked elements and set up signals for covered/uncovered
    // elements
    if (Intersect(*element, mask_)) {
      element->Cover();
      element->SaveState();
      if (uncovered_.count(element->name()) == 0) {
        uncovered_[element->name()] = 0;
      }
    } else {
      if (uncovered_.count(element->name()) == 0) {
        uncovered_[element->name()] = 1;
      } else {
        uncovered_[element->name()]++;
      }
    }
  }
}

void Polymer::Bind(Polymerase::Ptr pol, const std::string &promoter_name) {
  // Make a list of free promoters that pol can bind
  bool found = false;
  Element::VecPtr element_choices;
  for (auto &elem : elements_) {
    if (elem->name() == promoter_name && !elem->IsCovered()) {
      element_choices.push_back(elem);
      found = true;
    }
  }
  // Error checking
  if (!found) {
    std::string err = "Polymerase " + pol->name() +
                      " could not find free promoter " + promoter_name +
                      " to bind in the polyemr " + name_;
    throw std::runtime_error(err);
  }
  // Randomly select promoter.
  Element::Ptr elem = Random::WeightedChoice(element_choices);
  // More error checking.
  if (!elem->CheckInteraction(pol->name())) {
    std::string err = "Polymerase " + pol->name() +
                      " does not interact with promoter " + promoter_name;
    throw std::runtime_error(err);
  }
  // Update polymerase coordinates
  // (TODO: refactor; pol doesn't need to expose footprint/stop position)
  pol->set_start(elem->start());
  pol->set_stop(elem->start() + pol->footprint() - 1);
  // Find index of element
  auto it = std::find(elements_.begin(), elements_.end(), elem);
  pol->set_left_most_element(it - elements_.begin());
  // More error checking.
  if (pol->stop() > mask_.start()) {
    std::string err = "Polymerase " + pol->name() +
                      " will overlap with mask upon promoter binding. This may "
                      "cause the polymerase to stall and produce unexpected "
                      "behavior.";
    throw std::runtime_error(err);
  }
  elem->Cover();
  elem->SaveState();
  // Cover promoter in cache
  CoverElement(elem->name());
  // Add polymerase to this polymer
  Insert(pol);

  // Report some data to tracker
  if (elem->name() == "rbs") {
    auto &tracker = SpeciesTracker::Instance();
    tracker.IncrementRibo(elem->gene(), 1);
  }
}

void Polymer::Execute() {
  if (prop_sum_ == 0) {
    throw std::runtime_error(
        "Attempting to execute polymer with reaction propensity of 0.");
  }
  int pol_index = Choose();
  Move(pol_index);
}

void Polymer::ShiftMask() {
  if (mask_.start() >= mask_.stop()) {
    return;
  }
  int index = -1;
  for (size_t i = 0; i < elements_.size(); i++) {
    if (Intersect(mask_, *elements_[i])) {
      elements_[i]->SaveState();
      elements_[i]->Uncover();
      index = i;
      break;
    }
  }
  mask_.Recede();
  if (index == -1) {
    return;
  }
  if (Intersect(mask_, *elements_[index])) {
    elements_[index]->Cover();
  }
  elements_[index]->CheckState();
  elements_[index]->SaveState();
}

int Polymer::PolymeraseIndex(Polymerase::Ptr pol) {
  auto it = std::find(polymerases_.begin(), polymerases_.end(), pol);
  int index = it - polymerases_.begin();
  return index;
}

void Polymer::Terminate(Polymerase::Ptr pol, const std::string &last_gene) {
  auto it = std::find(polymerases_.begin(), polymerases_.end(), pol);
  int index = it - polymerases_.begin();
  prop_sum_ -= prop_list_[index];
  termination_signal_.Emit(index_, pol->name(), last_gene);
  polymerases_.erase(it);
  prop_list_.erase(prop_list_.begin() + index);
  if (prop_list_.size() != polymerases_.size()) {
    throw std::runtime_error("Prop list not correct size.");
  }
}

void Polymer::CoverElement(const std::string &species_name) {
  uncovered_[species_name]--;
  if (uncovered_[species_name] < 0) {
    std::string err = "Cached count of uncovered element " + species_name +
                      " cannot be a negative value";
    throw std::runtime_error(err);
  }
}

void Polymer::UncoverElement(const std::string &species_name) {
  uncovered_[species_name]++;
}

void Polymer::Insert(Polymerase::Ptr pol) {
  // Find where in vector polymerase should insert; use a lambda function
  // to make comparison between polymerase pointers
  auto it = std::upper_bound(polymerases_.begin(), polymerases_.end(), pol,
                             [](Polymerase::Ptr a, Polymerase::Ptr b) {
                               return a->start() < b->start();
                             });
  // Record position for prop_list_
  // NOTE: iterators become invalid as soon as a vector is changed!!
  // Attempting to use an iterator twice will lead to a segfault.
  auto prop_it = (it - polymerases_.begin()) + prop_list_.begin();
  // Add polymerase to this polymer
  polymerases_.insert(it, pol);
  // Cache polymerase speed
  double weight = weights_[pol->stop() - start_];
  // Update total move propensity of this polymer
  prop_sum_ += weight * pol->speed();
  prop_list_.insert(prop_it, weight * pol->speed());
  if (prop_list_.size() != polymerases_.size()) {
    throw std::runtime_error("Prop list not correct size.");
  }
}

int Polymer::Choose() {
  double prop_list_sum;
  for (auto &n : prop_list_) {
    prop_list_sum += n;
  }
  if (prop_list_sum != prop_sum_) {
    throw std::runtime_error("Prop sum and list not equal " +
                             std::to_string(prop_list_sum) + " " +
                             std::to_string(prop_sum_));
  }
  if (prop_list_.size() == 0) {
    std::string err = "There are no active polymerases on polymer " + name_ +
                      std::to_string(prop_sum_);
    throw std::runtime_error(err);
  }
  return Random::WeightedChoiceIndex(polymerases_, prop_list_);
}

void Polymer::Move(int pol_index) {
  // Error checking to make sure that pol is in vector
  if (pol_index >= polymerases_.size()) {
    std::string err = "Attempting to move unbound polymerase with index " +
                      std::to_string(pol_index) + " on polymer " + name_;
    throw std::runtime_error(err);
  }
  if (pol_index >= prop_list_.size()) {
    throw std::runtime_error(
        "Prop list vector index is invalid (before move).");
  }
  // Find which elements this polymerase is covering and temporarily uncover
  // them
  auto pol = polymerases_[pol_index];
  UncoverElements(pol);
  // Move polymerase
  pol->Move();
  // Resolve any collisions between polymerases or with mask
  bool pol_collision = ResolveCollisions(pol);
  bool mask_collision = ResolveMaskCollisions(pol);

  if (pol_index >= prop_list_.size()) {
    throw std::runtime_error(
        "Prop list vector index is invalid (before move 2).");
  }

  // Check for uncoverings
  bool terminated = RecoverElements(pol);

  // Check to see if it's safe to broadcast that this polymerase has moved
  if (!pol_collision && !mask_collision && !terminated) {
    if (pol_index >= prop_list_.size()) {
      throw std::runtime_error(
          "Prop list vector index is invalid (before move 3).");
    }
    // Terminate polymerase if it's run off the end of the polymer
    if (pol->stop() > stop_) {
      // TODO: Why are we sending the stop position of the polymer ?
      Terminate(pol, std::to_string(stop_));
    } else {
      // Update propensity for new codon
      if ((pol->stop() - start_) >= weights_.size()) {
        throw std::runtime_error("Weight is missing for this position.");
      }

      double prop_list_sum = 0;
      for (auto &n : prop_list_) {
        prop_list_sum += n;
      }

      // std::cout << "Just moved! Prop sum and list not equal " +
      //                  std::to_string(prop_list_sum) + " " +
      //                  std::to_string(prop_sum_);
      if (prop_list_.size() != polymerases_.size()) {
        throw std::runtime_error("Prop list not correct size.");
      }
      if ((pol->stop() - start_) >= weights_.size()) {
        throw std::runtime_error("Weights vector index is invalid.");
      }
      if (pol_index >= prop_list_.size()) {
        throw std::runtime_error("Prop list vector index is invalid.");
      }

      double weight = weights_[pol->stop() - start_];
      double new_speed = weight * pol->speed();
      // std::cout << std::to_string(pol->speed());
      double diff = new_speed - prop_list_[pol_index];
      prop_sum_ += diff;
      prop_list_[pol_index] = new_speed;
      pol->move_signal_.Emit();
      prop_list_sum = 0;
      for (auto &n : prop_list_) {
        prop_list_sum += n;
      }
      if (prop_list_sum != prop_sum_) {
        throw std::runtime_error(
            "Just moved! Prop sum and list not equal " +
            std::to_string(prop_list_sum) + " " + std::to_string(prop_sum_) +
            " " + std::to_string(diff) + " " + std::to_string(new_speed) + " " +
            std::to_string(weight));
      }
    }
  }
}

void Polymer::UncoverElements(Polymerase::Ptr pol) {
  int save_index = pol->left_most_element();
  if (save_index < 0) {
    std::string err = "`left_most_element` for polymerase " + pol->name() +
                      " is not defined.";
    throw std::runtime_error(err);
  }
  while (elements_[save_index]->start() <= pol->stop()) {
    if (Intersect(*pol, *elements_[save_index])) {
      elements_[save_index]->SaveState();
      elements_[save_index]->Uncover();
    }
    save_index++;
    if (save_index >= elements_.size()) {
      break;
    }
  }
}

bool Polymer::RecoverElements(Polymerase::Ptr pol) {
  int old_index = pol->left_most_element();
  if (old_index < 0) {
    std::string err = "`left_most_element` for polymerase " + pol->name() +
                      " is not defined.";
    throw std::runtime_error(err);
  }
  bool reset_index = true;
  bool terminating = false;
  while (elements_[old_index]->start() <= pol->stop()) {
    if (Intersect(*pol, *elements_[old_index])) {
      if (reset_index) {
        pol->set_left_most_element(old_index);
        reset_index = false;
      }
      if (terminating) {
        elements_[old_index]->Uncover();
      } else {
        elements_[old_index]->Cover();
        if (ResolveTermination(pol, elements_[old_index])) {
          old_index = pol->left_most_element() - 1;
          terminating = true;
        }
      }
    }
    elements_[old_index]->CheckState();
    elements_[old_index]->SaveState();
    old_index++;
    if (old_index >= elements_.size()) {
      break;
    }
  }
  return terminating;
}

bool Polymer::ResolveTermination(Polymerase::Ptr pol, Element::Ptr element) {
  // TODO: find less janky way to do this
  Terminator::Ptr term = std::dynamic_pointer_cast<Terminator>(element);
  if (element->type() != "terminator") {
    return false;
  } else {
    if (!term->CheckInteraction(pol->name(), pol->reading_frame())) {
      return false;
    }
    if (term->readthrough()) {
      return false;
    }
    if (pol->stop() != term->start()) {
      return false;
    }
    double random_num = Random::random();
    if (random_num <= term->efficiency(pol->name())) {
      if (PolymeraseIndex(pol) >= prop_list_.size()) {
        throw std::runtime_error(
            "Prop list vector index is invalid (before termination).");
      }
      // Fire Emit signal until entire terminator is uncovered
      int dist = term->stop() - pol->stop() + 1;
      for (int i = 0; i < dist; i++) {
        pol->move_signal_.Emit();
      }
      Terminate(pol, term->gene());
      return true;
    } else {
      term->set_readthrough(true);
      return false;
    }
  }
}

bool Polymer::ResolveMaskCollisions(Polymerase::Ptr pol) {
  if (mask_.start() > stop_) {
    // Is there still a mask?
    return false;
  }
  if (Intersect(*pol, mask_)) {
    if (pol->stop() - mask_.start() > 1) {
      std::string err =
          "Polymerase " + pol->name() +
          " is overlapping mask by more than one position on polymer";
      throw std::runtime_error(err);
    }
    if (mask_.CheckInteraction(pol->name())) {
      ShiftMask();
    } else {
      pol->MoveBack();
      return true;
    }
  }
  return false;
}

bool Polymer::ResolveCollisions(Polymerase::Ptr pol) {
  bool collision = false;
  auto it = std::find(polymerases_.begin(), polymerases_.end(), pol);
  int index = it - polymerases_.begin();
  if (index + 1 > polymerases_.size() - 1) {
    return collision;
  }
  // We only need to check the polymerase one position ahead of this polymerase
  if (Intersect(*pol, *polymerases_[index + 1])) {
    if (pol->stop() > polymerases_[index + 1]->start()) {
      std::string err =
          "Polymerase " + pol->name() +
          " (start: " + std::to_string(pol->start()) +
          ", stop: " + std::to_string(pol->stop()) +
          ", index: " + std::to_string(index) + ") is overlapping polymerase " +
          polymerases_[index + 1]->name() +
          " (start: " + std::to_string(polymerases_[index + 1]->start()) +
          ", stop: " + std::to_string(polymerases_[index + 1]->stop()) +
          ", index: " + std::to_string(index + 1) +
          ") by more than one position on polymer " + name_;
      throw std::runtime_error(err);
    }
    pol->MoveBack();
    collision = true;
  }
  return collision;
}

bool Polymer::Intersect(const Feature &elem1, const Feature &elem2) {
  return (elem1.stop() >= elem2.start()) && (elem2.stop() >= elem1.start());
}

Transcript::Transcript(const std::string &name, int start, int stop,
                       const Element::VecPtr &elements, const Mask &mask,
                       const std::vector<double> &weights)
    : Polymer(name, start, stop, elements, mask, weights) {}

void Transcript::Bind(Polymerase::Ptr pol, const std::string &promoter_name) {
  // Bind polymerase just like in parent Polymer
  Polymer::Bind(pol, promoter_name);
  // Set the reading frame of the polymerase
  // TODO: should the reading frame be set by the polymerase start position?
  pol->set_reading_frame(pol->start() % 3);
}

Genome::Genome(const std::string &name, int length,
               const Element::VecPtr &elements,
               const Element::VecPtr &transcript_template, const Mask &mask)
    : Polymer(name, 1, length, elements, mask),
      transcript_template_(transcript_template) {
  transcript_weights_ = std::vector<double>(length + 1, 1.0);
  // Sort transcript template
  std::sort(transcript_template_.begin(), transcript_template_.end(),
            [](Element::Ptr elem1, Element::Ptr elem2) {
              return elem1->start() < elem2->start();
            });
  // Sort genomic elements
  std::sort(elements_.begin(), elements_.end(),
            [](Element::Ptr elem1, Element::Ptr elem2) {
              return elem1->start() < elem2->start();
            });
}

Genome::Genome(const std::string &name, int length,
               const Element::VecPtr &elements,
               const Element::VecPtr &transcript_template, const Mask &mask,
               const std::vector<double> &transcript_weights)
    : Polymer(name, 1, length, elements, mask),
      transcript_template_(transcript_template),
      transcript_weights_(transcript_weights) {
  // Sort transcript template
  std::sort(transcript_template_.begin(), transcript_template_.end(),
            [](Element::Ptr elem1, Element::Ptr elem2) {
              return elem1->start() < elem2->start();
            });
  // Sort genomic elements
  std::sort(elements_.begin(), elements_.end(),
            [](Element::Ptr elem1, Element::Ptr elem2) {
              return elem1->start() < elem2->start();
            });
}

void Genome::Bind(Polymerase::Ptr pol, const std::string &promoter_name) {
  // Bind polymerase
  Polymer::Bind(pol, promoter_name);
  // Construct a transcript starting from *end* of polymerase
  Transcript::Ptr transcript = BuildTranscript(pol->stop(), stop_);
  // Connect polymerase movement signal to transcript, so that transcript knows
  // when to expose new elements
  // TODO: figure out if this could cause a memory leak
  // What happens if transcript gets deleted?
  pol->move_signal_.ConnectMember(transcript, &Transcript::ShiftMask);
  // Fire new transcript signal (adds transcript to Simulation)
  transcript_signal_.Emit(transcript);
}

Transcript::Ptr Genome::BuildTranscript(int start, int stop) {
  Element::VecPtr elements;
  for (const auto &elem : transcript_template_) {
    if (elem->start() >= start && elem->stop() <= stop) {
      // Construct a *copy* and insert into elements
      elements.emplace_back(elem->Clone());
    }
  }
  Transcript::Ptr transcript;
  Mask mask = Mask("mask", start, stop, std::vector<std::string>());
  // We need to used the standard shared_ptr constructor here because the
  // constructor of Transcript needs to know its address in memory to wire
  // signals appropriately.
  transcript = std::make_shared<Transcript>("rna", start_, stop_, elements,
                                            mask, transcript_weights_);
  return transcript;
}