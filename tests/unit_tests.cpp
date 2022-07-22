#include "./lib/catch.hpp"
#include "choices.hpp"
#include "feature.hpp"
#include "model.hpp"
#include "polymer.hpp"
#include "reaction.hpp"
#include "tracker.hpp"

TEST_CASE("Genome construction")
{
    auto plasmid = std::shared_ptr<Genome>(new Genome("T7", 305, 0.001, 9.0, 20, 0.01));
    auto plasmid2 = std::shared_ptr<Genome>(new Genome("T7", 305));
    
    REQUIRE(plasmid->transcript_degradation_rate() == 0.01);
    REQUIRE(plasmid->transcript_degradation_rate_ext() == 0.001);
    REQUIRE(plasmid->rnase_speed() == 9);
    REQUIRE(plasmid->rnase_footprint() == 20);

    //Test the genome default constructor
    REQUIRE(plasmid2->transcript_degradation_rate() == 0.0);
    REQUIRE(plasmid2->transcript_degradation_rate_ext() == 0.0);
    REQUIRE(plasmid2->rnase_speed() == 0.0);
    REQUIRE(plasmid2->rnase_footprint() == 0);
}

TEST_CASE("Add mask to genome")
{
    auto plasmid = std::shared_ptr<Genome>(new Genome("T7", 305));
    std::vector<std::string> interactions = {"rnapol", "ecolipol"};
    plasmid->AddMask(15, interactions);

    REQUIRE(plasmid->GetMask().start() == 15);
    REQUIRE(plasmid->GetMask().stop() == 305);
    REQUIRE(plasmid->GetMask().CheckInteraction("rnapol"));
    REQUIRE(plasmid->GetMask().CheckInteraction("ecolipol"));
}

TEST_CASE("Add fixed elements to a genome")
{
    //Create a genome and add a promoter
    auto plasmid = std::shared_ptr<Genome>(new Genome("T7", 305));
    std::map<std::string, double> interactions = {{"rnapol", 2e8}};
    plasmid->AddPromoter("phi1", 1, 10, interactions);
    CHECK(plasmid->bindings().size() == 1);
    CHECK(plasmid->GetBindingIntervals().size() == 1);
    
    //Test the AddTerminator method
    std::map<std::string, double> efficiency = {{"rnapol", 1.0}};
    plasmid->AddTerminator("t1", 304, 305, efficiency);
    CHECK(plasmid->GetReleaseIntervals().size() == 1);

    //Add two genes to the genome
    //There should now be three binding sites ('bindings') -
    //one polymerase binding site and two ribosome binding sites
    plasmid->AddGene("rnapol", 26, 225, 11, 26, 1e7);
    plasmid->AddGene("proteinX", 241, 280, 226, 241, 1e7);
    REQUIRE(plasmid->bindings().size() == 3);
}

TEST_CASE("MobileElementManager Insert") 
{
    //The MobileElementManager constructor expects a vector of weights
    //but the values have no bearing on this particular test
    std::vector<double> weights = {30, 1.0};
    MobileElementManager manager = MobileElementManager(weights);

    //Make two polymerases and give them unique start positions

    // TANVI'S EDITED SECTION
    // @TODO: Change Polymerase constructor function calls for these rnapols
    // Old code commented out 
    //auto polymerase1 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));
    //auto polymerase2 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));

    auto polymerase1 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));
    auto polymerase2 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));

    polymerase1->start(1); polymerase2->start(51);
    
    manager.Insert(polymerase1, std::shared_ptr<Polymer>());
    manager.Insert(polymerase2, std::shared_ptr<Polymer>());

    //The MobileElementManager should have the polymerases ordered by their start 
    //positions
    REQUIRE(manager.pol_count() == 2);
    REQUIRE(manager.GetPol(0)->start() == 1);
    REQUIRE(manager.GetPol(1)->start() == 51);
}

TEST_CASE("MobileElementManager Delete")
{
    //The MobileElementManager constructor expects a vector of weights
    //but the values have no bearing on this particular test
    std::vector<double> weights = {30, 1.0};
    MobileElementManager manager = MobileElementManager(weights);
    
    // TANVI'S EDITED SECTION
    // @TODO: Change Polymerase constructor function calls for these rnapols
    // Old code removes 0 and false from Polymerase constructor call. 
    auto polymerase1 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));
    auto polymerase2 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));
    auto polymerase3 = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));
    polymerase1->start(1); polymerase2->start(51); polymerase3->start(75);
    
    manager.Insert(polymerase1, std::shared_ptr<Polymer>());
    manager.Insert(polymerase2, std::shared_ptr<Polymer>());
    manager.Insert(polymerase3, std::shared_ptr<Polymer>());
    CHECK(manager.pol_count() == 3);
    manager.Delete(1);

    REQUIRE(manager.pol_count() == 2);
    REQUIRE(manager.GetPol(0)->start() == 1);
    REQUIRE(manager.GetPol(1)->start() == 75);

    manager.Delete(0);
    REQUIRE(manager.pol_count() == 1);
    REQUIRE(manager.GetPol(0)->start() == 75);
}

//Tests the Polymer::Bind method, which takes a mobile element and
//the name of a target binding site and tries to attach the mobile element to 
//the target
TEST_CASE("Attach a polymerase to a registered genome")
{
    //Setup - create a genome with a single promoter that interacts with rnapol
    //The genome must be registered with a model object before Bind can be called -
    //failing to do so results in undefined behavior
    int promoter_start = 2;
    auto sim = std::shared_ptr<Model>(new Model(1.1e-15));
    auto plasmid = std::shared_ptr<Genome>(new Genome("T7", 305));
    std::map<std::string, double> interactions = {{"rnapol", 2e8}};
    std::map<std::string, double> efficiency = {{"rnapol", 1.0}};
    plasmid->AddPromoter("phi1", promoter_start, 10, interactions);
    sim->RegisterGenome(plasmid);

    //Create a polymerase called 'rnapol' and bind it to the genome
    // TANVI'S EDITED SECTION - old Polymerase constructor code does not have 0 and false as arguments
    auto polymerase = std::make_shared<Polymerase>(Polymerase("rnapol", 10, 40));
    plasmid->Bind(polymerase, "phi1");
    
    //The genome should now have a single mobile element attached to it that
    //has the same start position as the promoter
    CHECK(plasmid->num_attached() == 1);
    REQUIRE(plasmid->attached_pol_start(0) == promoter_start);
}
