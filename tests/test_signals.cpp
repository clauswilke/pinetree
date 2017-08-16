#include <catch.hpp>

#include <iostream>
#include <string>
#include <vector>

#include "event_signal.hpp"

TEST_CASE("Signal connection and firing", "[Signals]") {
  class Button {
  public:
    Signal<> on_click;
  };

  class Message {
  public:
    bool fired = false;
    void display() {
      fired = true;
      // std::cout << "Hello World!" << std::endl;
    }
  };

  auto button = std::make_shared<Button>();
  auto message = std::make_shared<Message>();

  button->on_click.ConnectMember(message, &Message::display);
  auto button2 = button;
  button2->on_click.Emit();
  REQUIRE(message->fired == true);
}