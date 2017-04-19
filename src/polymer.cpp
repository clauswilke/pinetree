#include "polymer.hpp"
#include "choices.hpp"

Polymer::Polymer(const std::string &name, int start, int stop,
                 const Element::VecPtr &elements,
                 const Mask &mask)
    : index_(0),
      name_(name),
      start_(start),
      stop_(stop),
      elements_(elements),
      mask_(mask),
      prop_sum_(0)
{
    for (auto &element : elements_)
    {
        element->cover_signal_.ConnectMember(this, &Polymer::CoverElement);
        element->uncover_signal_.ConnectMember(this, &Polymer::UncoverElement);
        // Cover masked elements and set up signals for covered/uncovered
        // elements
        if (Intersect(*element, mask_))
        {
            element->Cover();
            element->SaveState();
            if (uncovered_.count(element->name()) == 0)
            {
                uncovered_[element->name()] = 0;
            }
        }
        else
        {
            if (uncovered_.count(element->name()) == 0)
            {
                uncovered_[element->name()] = 1;
            }
            else
            {
                uncovered_[element->name()]++;
            }
        }
    }
}

void Polymer::Bind(Polymerase::Ptr pol,
                   const std::string &promoter_name)
{
    // Make a list of free promoters that pol can bind
    bool found = false;
    Element::VecPtr element_choices;
    for (auto &elem : elements_)
    {
        if (elem->name() == promoter_name && !elem->IsCovered())
        {
            element_choices.push_back(elem);
            found = true;
        }
    }
    // Error checking.
    // Randomly select promoter.
    Element::Ptr elem = Random::WeightedChoice(element_choices);
    // More error checking.
    // Update polymerase coordinates
    // (TODO: refactor; pol doesn't need to expose footprint/stop position)
    pol->set_start(elem->start());
    pol->set_stop(elem->start() + pol->footprint() - 1);
    // Find index of element
    auto it = std::find(elements_.begin(), elements_.end(), elem);
    pol->set_left_most_element(it - elements_.begin());
    // More error checking.
    elem->Cover();
    elem->SaveState();
    // Cover promoter in cache
    CoverElement(elem->name());
    // Add polymerase to this polymer
    Insert(pol);
    // Update total move propensity of this polymer
    prop_sum_ += pol->speed();
}

void Polymer::ShiftMask()
{
    if (mask_.start() >= mask_.stop())
    {
        return;
    }

    int index = -1;
    for (size_t i = 0; i < elements_.size(); i++)
    {
        if (Intersect(mask_, *elements_[i]))
        {
            elements_[i]->SaveState();
            elements_[i]->Uncover();
            index = i;
            break;
        }
    }
    mask_.Recede();
    if (index == -1)
    {
        return;
    }
    if (Intersect(mask_, *elements_[index]))
    {
        elements_[index]->Cover();
    }
    elements_[index]->CheckState();
    elements_[index]->SaveState();
}

void Polymer::Insert(Polymerase::Ptr pol)
{
    auto it = std::upper_bound(polymerases_.begin(), polymerases_.end(), pol);
    polymerases_.insert(it, pol);
    auto prop_it = prop_list_.begin() + (it - polymerases_.begin());
    prop_list_.insert(prop_it, pol->speed());
}

bool Polymer::Intersect(const Feature &elem1, const Feature &elem2)
{
    return (elem1.stop() >= elem2.start()) && (elem2.stop() >= elem1.start());
}