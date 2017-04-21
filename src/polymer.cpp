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
  if (!CheckInteraction(pol->name())) {
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
}

void Polymer::UncoverElement(const std::string &species_name) {
  uncovered_[species_name]++;
}

void Polymer::Insert(Polymerase::Ptr pol) {
  auto it = std::upper_bound(polymerases_.begin(), polymerases_.end(), pol);
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
  return Random::WeightedChoice(polymerases_, prop_list_);
}

void Polymer::UncoverElements(Polymerase::Ptr pol) {
  int save_index = pol->left_most_element();
  // TODO: Error checking
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

bool Polymer::ResolveTermination(Polymerase::Ptr pol, Terminator::Ptr element) {
  if (!element->CheckInteraction(pol->name(), pol->reading_frame())) {
    return false;
  }
  if (element->readthrough()) {
    return false;
  }
  if (pol->stop() != element->start()) {
    return false;
  }
  double random_num = Random::random();
  if (random_num <= element->efficiency(pol->name())) {
    Terminate(pol, element->gene());
    return true;
  } else {
    element->set_readthrough(true);
    return false;
  }
}

bool Polymer::ResolveMaskCollisions(Polymerase::Ptr pol) {}

bool Polymer::Intersect(const Feature &elem1, const Feature &elem2) {
  return (elem1.stop() >= elem2.start()) && (elem2.stop() >= elem1.start());
}