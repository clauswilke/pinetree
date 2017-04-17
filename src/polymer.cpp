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
        element.cover_signal.connect_member(this, &Polymer::cover_element);
    }
}