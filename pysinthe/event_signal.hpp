/**
 * Derived from this blog post:
 * http://simmesimme.github.io/tutorials/2015/09/20/signal-slot
 */

#ifndef SIGNAL_HPP
#define SIGNAL_HPP

#include <functional>
#include <map>
#include <memory>

// A signal object may call multiple slots with the
// same signature. You can connect functions to the signal
// which will be called when the emit() method on the
// signal object is invoked. Any argument passed to emit()
// will be passed to the given functions.

template <typename... Args> class Signal {

public:
  Signal() : current_id_(0) {}

  // copy creates new signal
  Signal(Signal const &other) : current_id_(0) {}

  // connects a member function to this Signal
  template <typename T> int ConnectMember(T *inst, void (T::*func)(Args...)) {
    return Connect([=](Args... args) { (inst->*func)(args...); });
  }

  // connects a const member function to this Signal
  template <typename T>
  int ConnectMember(T *inst, void (T::*func)(Args...) const) {
    return Connect([=](Args... args) { (inst->*func)(args...); });
  }

  // smart pointer support
  template <typename T>
  int ConnectMember(std::shared_ptr<T> inst, void (T::*func)(Args...)) {
    return Connect([=](Args... args) { std::mem_fn(func)(inst, args...); });
  }

  // smart pointer support
  template <typename T>
  int ConnectMember(std::shared_ptr<T> inst, void (T::*func)(Args...) const) {
    return Connect([=](Args... args) { std::mem_fn(func)(inst, args...); });
  }

  // connects a std::function to the signal. The returned
  // value can be used to disconnect the function again
  int Connect(std::function<void(Args...)> const &slot) const {
    slots_.insert(std::make_pair(++current_id_, slot));
    return current_id_;
  }

  // disconnects a previously connected function
  void Disconnect(int id) const { slots_.erase(id); }

  // disconnects all previously connected functions
  void DisconnectAll() const { slots_.clear(); }

  // calls all connected functions
  void Emit(Args... p) {
    for (auto it : slots_) {
      it.second(p...);
    }
  }

  // assignment creates new Signal
  Signal &operator=(Signal const &other) { DisconnectAll(); }

private:
  mutable std::map<int, std::function<void(Args...)>> slots_;
  mutable int current_id_;
};

#endif /* SIGNAL_HPP */