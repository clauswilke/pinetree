#include <catch.hpp>
#include "choices.hpp"
#include "feature.hpp"
#include "model.hpp"
#include "polymer.hpp"
#include "reaction.hpp"
#include "tracker.hpp"

//Test genome constructor
TEST_CASE("Genome Construction")
{
    Genome plasmid = Genome("T7", 305, 0.01, 0.001, 9.0, 20.0);
    Genome plasmid2 = Genome("T7", 305);
    
    REQUIRE(plasmid.transcript_degradation_rate() == 0.01);
    REQUIRE(plasmid.transcript_degradation_rate_ext() == 0.001);
    REQUIRE(plasmid.rnase_speed() == 9);
    REQUIRE(plasmid.rnase_footprint() == 20);

    REQUIRE(plasmid2.transcript_degradation_rate() == 0.0);
    REQUIRE(plasmid2.transcript_degradation_rate_ext() == 0.0);
    REQUIRE(plasmid2.rnase_speed() == 0.0);
    REQUIRE(plasmid2.rnase_footprint() == 0);
}