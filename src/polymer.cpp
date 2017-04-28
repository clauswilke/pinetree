#include "polymer.hpp"
#include "choices.hpp"

Polymer::Polymer(const std::string &name, int start, int stop,
                 const Element::VecPtr &elements, const Mask &mask)
    : index_(0), name_(name), start_(start), stop_(stop), elements_(elements),
      mask_(mask), prop_sum_(0) {
  for (auto &element : elements_) {
    element->cover_signal_.ConnectMember(this, &Polymer::CoverElement);
    element->uncover_signal_.ConnectMember(this, &Polymer::UncoverElement);
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
  // // Update total move propensity of this polymer
  prop_sum_ += pol->speed();
}

void Polymer::Execute() {
  if (prop_sum_ == 0) {
    throw std::runtime_error(
        "Attempting to execute polymer with reaction propensity of 0.");
  }
  auto pol = Choose();
  Move(pol);
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

void Polymer::Terminate(Polymerase::Ptr pol, const std::string &last_gene) {
  prop_sum_ -= pol->speed();
  auto it = std::find(polymerases_.begin(), polymerases_.end(), pol);
  int index = it - polymerases_.begin();
  termination_signal_.Emit(index, pol->name(), last_gene);
  polymerases_.erase(it);
  prop_list_.erase(prop_list_.begin() + index);
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
  // Find where in vector polymerase should insert; use a lambada function
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
  prop_list_.insert(prop_it, pol->speed());
}

Polymerase::Ptr Polymer::Choose() {
  if (prop_list_.size() == 0) {
    std::string err = "There are no active polymerases on polymer " + name_;
    throw std::runtime_error(err);
  }
  return Random::WeightedChoice(polymerases_, prop_list_);
}

void Polymer::Move(Polymerase::Ptr pol) {
  // Error checking to make sure that pol is in vector
  auto it = std::find(polymerases_.begin(), polymerases_.end(), pol);
  if (it == polymerases_.end()) {
    std::string err = "Attempting to move unbound polymerase " + pol->name() +
                      " on polymer " + name_;
    throw std::runtime_error(err);
  }
  // Find which elements this polymerase is covering and temporarily uncover
  // them
  UncoverElements(pol);
  // Move polymerase
  pol->Move();
  // Resolve any collisions between polymerases or with mask
  bool pol_collision = ResolveCollisions(pol);
  bool mask_collision = ResolveMaskCollisions(pol);
  // Check to see if it's safe to broadcast that this polymerase has moved
  if (!pol_collision && !mask_collision) {
    pol->move_signal_.Emit();
  }
  // Check for uncoverings
  RecoverElements(pol);
  // Terminate polymerase if it's run off the end of the polymer
  if (pol->stop() > stop_) {
    // TODO: Why are we sending the stop position of the polymer ?
    Terminate(pol, std::to_string(stop_));
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

void Polymer::RecoverElements(Polymerase::Ptr pol) {
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
          "Polymerase " + pol->name() + " (start: " +
          std::to_string(pol->start()) + ", stop: " +
          std::to_string(pol->stop()) + ", index: " + std::to_string(index) +
          ") is overlapping polymerase " + polymerases_[index + 1]->name() +
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
                       const Element::VecPtr &elements, const Mask &mask)
    : Polymer(name, start, stop, elements, mask) {}

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
      transcript_template_(transcript_template) {}

void Genome::Bind(Polymerase::Ptr pol, const std::string &promoter_name) {
  // Bind polymerase
  Polymer::Bind(pol, promoter_name);
  // Construct a transcript starting from *end* of polymerase
  Transcript::Ptr transcript = BuildTranscript(pol->stop(), stop_);
  // Connect polymerase movement signal to transcript, so that transcript knows
  // when to expose new elements
  // TODO: figure out if this could cause a memory leak
  // What happens if transcript gets deleted?
  pol->move_signal_.ConnectMember(&(*transcript), &Transcript::ShiftMask);
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
  transcript =
      std::make_shared<Transcript>("rna", start_, stop_, elements, mask);
  return transcript;
}