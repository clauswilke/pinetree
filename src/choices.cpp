#include "choices.hpp"

std::mt19937 Random::gen_;
std::uniform_real_distribution<> Random::dis_(0, 1);

Random::Random()
{
    std::random_device rd;
    gen_.seed(rd());
}

void Random::seed(int seed)
{
    gen_.seed(seed);
}

double Random::random()
{
    return dis_(gen_);
}