/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <algorithm>
#include <string>
#include <vector>

#include "feature.hpp"

Feature::Feature(const std::string &name, int start, int stop,
                 const std::vector<std::string> &interactions)
    : name_(name),
      start_(start),
      stop_(stop),
      interactions_(interactions)
{
}

bool Feature::check_interaction(const std::string &name)
{
  return std::find(interactions_.begin(), interactions_.end(), name) !=
         interactions_.end();
}

Polymerase::Polymerase(const std::string &name, int footprint, int speed)
    : name_(name),
      footprint_(footprint),
      speed_(speed)
{
  start_ = 0;
  stop_ = footprint - 1;
  left_most_element_ = 0;
  bound_ = 0;
  type_ = "polymerase";
  reading_frame_ = 0;
}

void Polymerase::move()
{
  start_++;
  stop_++;
}

void Polymerase::move_back()
{
  start_--;
  stop_--;
}

Mask::Mask(const std::string &name, int start, int stop,
           const std::vector<std::string> &interactions)
    : Feature(name, start, stop, interactions)
{
}

Element::Element(const std::string &name, int start, int stop,
                 const std::vector<std::string> &interactions)
    : Feature(name, start, stop, interactions),
      covered_(0),
      old_covered_(0)
{
}

Promoter::Promoter(const std::string &name, int start, int stop,
                   const std::vector<std::string> &interactions)
    : Element(name, start, stop, interactions)
{
  type_ = "promoter";
}

void Promoter::check_state()
{
  if (was_covered())
  {
    cover_signal.emit(name_);
  }
  else if (was_uncovered())
  {
    uncover_signal.emit(name_);
  }
}

Terminator::Terminator(const std::string &name, int start, int stop,
                       const std::vector<std::string> &interactions)
    : Element(name, start, stop, interactions)
{
  readthrough_ = false;
  reading_frame_ = -1;
}

void Terminator::check_state()
{
  if (was_uncovered())
  {
    readthrough_ = false;
    uncover_signal.emit(name_);
  }
  else if (was_covered())
  {
    cover_signal.emit(name_);
  }
}

bool Terminator::check_interaction(const std::string &name, int reading_frame)
{
  return (reading_frame == reading_frame_ && Element::check_interaction(name));
}
