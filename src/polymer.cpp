#include "polymer.hpp"

Polymer::Polymer(const std::string &name, int start, int stop,
                 const std::vector<Element> &elements, const Mask &mask)
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
        element.cover_signal_.ConnectMember(this, &Polymer::CoverElement);
        element.uncover_signal_.ConnectMember(this, &Polymer::UncoverElement);
        // Cover masked elements and set up signals for covered/uncovered
        // elements
        if (Intersect(element, mask_))
        {
            element.Cover();
            element.SaveState();
            if (uncovered_.count(element.name()) == 0)
            {
                uncovered_[element.name()] = 0;
            }
        }
        else
        {
            if (uncovered_.count(element.name()) == 0)
            {
                uncovered_[element.name()] = 1;
            }
            else
            {
                uncovered_[element.name()]++;
            }
        }
    }
}

void Polymer::Insert(const Polymerase &pol)
{
    auto it = std::upper_bound(polymerases_.begin(), polymerases_.end(), pol);
    polymerases_.insert(it, pol);
    auto prop_it = prop_list_.begin() + (it - polymerases_.begin());
    prop_list_.insert(prop_it, pol.speed());
}

bool Polymer::Intersect(const Feature &elem1, const Feature &elem2)
{
    return (elem1.stop() >= elem2.start()) && (elem2.stop() >= elem1.start());
}