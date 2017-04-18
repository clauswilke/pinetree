/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <catch.hpp>

#include <string>
#include <vector>

#include "choices.hpp"

TEST_CASE("Random number generation and weighted choice", "[Random]")
{
    Random rng = Random();
    rng.seed(5);
    REQUIRE(rng.random() == Approx(0.0551801));

    std::vector<int> pop;
    for (int i = 0; i < 51; i++)
    {
        pop.push_back(i);
    }
    std::vector<double> weights(pop.begin(), pop.end());
    int out = rng.WeightedChoice(pop, weights);
    REQUIRE(out == 46);
    out = rng.WeightedChoice(pop, weights);
    REQUIRE(out == 30);
}