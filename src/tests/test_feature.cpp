#include <catch.h>
#include <string>
#include <vector>
#include "../feature.h"

TEST_CASE("Feature construction.", "[feature]") {
  std::string name = "testing!";
  int start = 1;
  int stop = 10;
  std::vector<std::string> interactions = {"ecolipol", "rnapol"};
  Feature feat = Feature(name, start, stop, interactions);
  REQUIRE(feat.check_interaction(std::string("ecolipol")));
}
