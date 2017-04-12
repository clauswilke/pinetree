import unittest

from unittest.mock import MagicMock
import random
from pysinthe import feature, polymer, _USE_CPP_


class TestPolymerMethods(unittest.TestCase):
    # # #
    # SET UP AND HELPER FUNCITONS
    # # #
    def setUp(self):
        # Set up promoter
        self.promoter = feature.Promoter("p1", 5, 15, ["ecolipol", "rnapol"])
        self.promoter.cover_signal.connect(self.fire_block)
        self.promoter.uncover_signal.connect(self.fire_promoter)
        # Set up terminator
        self.terminator = feature.Terminator("t1", 50, 55,
                                             {"rnapol": {"efficiency": 1.0},
                                              "ecolipol": {"efficiency": 0.6}
                                              })
        # Set up mask
        self.mask = feature.Mask("m1", 10, 100, ["ecolipol"])
        # Set up polymerases
        self.pol1 = feature.Polymerase("ecolipol", 10, 30)
        self.pol1.start = 20
        self.pol1.stop = 30
        self.pol2 = feature.Polymerase("rnapol", 10, 30)
        self.pol2.start = 60
        self.pol2.stop = 70
        self.pol3 = feature.Polymerase("rnapol", 10, 30)
        self.pol3.start = 40
        self.pol3.stop = 50

        # Set up attributes for detecting when signals have been fired
        self.termination_fired = False
        self.promoter_fired = 0
        self.block_fired = 0

    def construct_polymer(self):
        # Construct a polymer from elements set up in setUp()
        poly = polymer.Polymer("g1", 1, 100, [self.promoter, self.terminator],
                               self.mask)
        poly.termination_signal.connect(self.fire_termination)
        return poly

    def construct_polymer_multi_promoter(self):
        promoter1 = feature.Promoter("p1", 5, 15, ["ecolipol"])
        promoter2 = feature.Promoter("p2", 16, 20, [])
        promoter3 = feature.Promoter("p3", 21, 30, [])
        term = feature.Terminator("t1", 31, 33,
                                  {"ecolipol": {"efficiency": 1.0}}
                                  )
        elements = [promoter1, promoter2, promoter3, term]
        poly = polymer.Polymer("test", 1, 100, elements)
        return promoter1, promoter2, promoter3, term, poly

    def shift_mask_n(self, polymer, n):
        # Shift mask n times
        for i in range(n):
            polymer.shift_mask()

    def move_polymerase_n(self, polymer, pol, n):
        # Move polymerase n times
        for i in range(n):
            polymer._move_polymerase(pol)

    def fire_termination(self, index, pol_name, gene_name):
        self.termination_fired = True

    def fire_promoter(self, name):
        self.promoter_fired += 1

    def fire_block(self, name):
        self.block_fired += 1
    # # #
    # END OF SET UP AND HELPER FUNCTIONS
    # # #

    def test_init(self):
        polymer = self.construct_polymer()
        # Check that all appropriate elements have been covered
        self.assertTrue(self.terminator.is_covered())
        self.assertFalse(self.terminator.was_covered())
        self.assertFalse(self.terminator.was_uncovered())
        self.assertTrue(self.promoter.is_covered())
        self.assertFalse(self.promoter.was_covered())
        self.assertFalse(self.promoter.was_uncovered())
        self.assertEqual(polymer.uncovered["p1"], 0)
        self.assertEqual(polymer.uncovered["t1"], 0)

    def test_bind_polymerase(self):
        self.setUp()
        polymer = self.construct_polymer()
        random.seed(22)
        # Promoter should be covered and inaccessible
        self.assertRaises(RuntimeError,
                          polymer.bind_polymerase, self.pol1, "p1")
        # Shift mask back 10 positions
        self.shift_mask_n(polymer, 10)
        # Bind polymerase
        polymer.bind_polymerase(self.pol1, "p1")
        # Check changes in coverings and positions
        self.assertEqual(self.pol1.start, 5)
        self.assertEqual(self.pol1.stop, 14)
        self.assertEqual(polymer.uncovered["p1"], 0)
        self.assertEqual(polymer.prop_sum, 30)

    def test_execute(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Shift mask out of the way
        self.shift_mask_n(polymer, 100)
        # Arrange polymerases along polymer
        polymer.bind_polymerase(self.pol1, "p1")
        self.move_polymerase_n(polymer, self.pol1, 35)
        polymer.bind_polymerase(self.pol2, "p1")
        self.move_polymerase_n(polymer, self.pol2, 25)
        polymer.bind_polymerase(self.pol3, "p1")
        # Execute reactions (movement on polymerase)
        for i in range(30):
            polymer.execute()
        # Once all of the polymerases have terminated, execute should throw
        # a RuntimeError
        with self.assertRaises(RuntimeError):
            for i in range(30):
                polymer.execute()

    def test_shift_mask(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Record old mask position
        old_mask_start = polymer.mask.start
        # Shift mask and check that coordinates have changed
        polymer.shift_mask()
        self.assertEqual(old_mask_start + 1, polymer.mask.start)
        # Check that coverings have not changed
        self.assertTrue(polymer.elements[0].is_covered())
        self.assertTrue(polymer.elements[1].is_covered())
        # Shift mask 10 more spaces
        self.shift_mask_n(polymer, 10)
        # Check that promoter has been uncovered
        self.assertFalse(polymer.elements[0].is_covered())
        self.assertTrue(polymer.elements[1].is_covered())
        # Make sure mask can't get shifted beyond its end position
        self.shift_mask_n(polymer, 1000)
        self.assertEqual(polymer.mask.start, polymer.mask.stop)
        self.assertFalse(polymer.elements[1].is_covered())

    def test_terminate(self):
        self.setUp()
        polymer = self.construct_polymer()
        self.shift_mask_n(polymer, 10)
        polymer.bind_polymerase(self.pol1, "p1")
        # Record old propensity
        old_prop_sum = polymer.prop_sum
        polymer.terminate(self.pol1, "gene1")
        # Make sure propensity has changed and termination signal fired
        self.assertNotEqual(polymer.prop_sum, old_prop_sum)
        self.assertTrue(self.termination_fired)
        # Make sure pol object has actually been removed from polymer
        with self.assertRaises(ValueError):
            polymer.polymerases.index(self.pol1)

    def test_count_uncovered(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Check that cached counts match true counts
        count = 0
        for element in polymer.elements:
            if element.name == "p1" and not element.is_covered():
                count += 1
        self.assertEqual(polymer.count_uncovered("p1"), count)
        count = 0
        for element in polymer.elements:
            if element.name == "t1" and not element.is_covered():
                count += 1
        self.assertEqual(polymer.count_uncovered("t1"), count)

    def test_cover_element(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Initialize count
        polymer.uncovered["p1"] = 0
        # Shouldn't be allow uncovered counts to drop below zero
        with self.assertRaises(RuntimeError):
            polymer.cover_element("p1")
        # Covering element should only decrease count by 1
        polymer.uncovered["p1"] = 5
        polymer.cover_element("p1")
        self.assertEqual(polymer.uncovered["p1"], 4)

    def test_uncover_element(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Inverse of test_cover_element()
        polymer.uncovered["p1"] = 0
        polymer.uncover_element("p1")
        self.assertEqual(polymer.uncovered["p1"], 1)

    def test_calculate_propensity(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Make sure cached propensity matches true propensity
        prop_sum = 0
        for pol in polymer.polymerases:
            prop_sum += pol.speed
        self.assertEqual(prop_sum, polymer.calculate_propensity())

    def test_insert_polymerase(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Add polymerases to polymer
        polymer._insert_polymerase(self.pol2)
        self.assertEqual(polymer.polymerases.index(self.pol2), 0)
        # List of polymerases should be maintained in order that they are
        # positioned on polymer, from left to right
        polymer._insert_polymerase(self.pol1)
        self.assertEqual(polymer.polymerases.index(self.pol1), 0)
        self.assertEqual(polymer.polymerases.index(self.pol2), 1)
        polymer._insert_polymerase(self.pol3)
        self.assertEqual(polymer.polymerases.index(self.pol1), 0)
        self.assertEqual(polymer.polymerases.index(self.pol2), 2)
        self.assertEqual(polymer.polymerases.index(self.pol3), 1)
        # Trying to insert the same pol object twice should throw an error
        with self.assertRaises(RuntimeError):
            polymer._insert_polymerase(self.pol2)

    def test_choose_polymerase(self):
        self.setUp()
        polymer = self.construct_polymer()
        # There are no polymerases on the transcript so trying to randomly
        # choose a polymerase should throw an error
        with self.assertRaises(RuntimeError):
            polymer._choose_polymerase()
        # Move mask back
        self.shift_mask_n(polymer, 100)
        # Arrange polymerases along polymer
        polymer.bind_polymerase(self.pol1, "p1")
        self.move_polymerase_n(polymer, self.pol1, 35)
        polymer.bind_polymerase(self.pol2, "p1")
        self.move_polymerase_n(polymer, self.pol2, 25)
        polymer.bind_polymerase(self.pol3, "p1")
        # Set seed and randomly select polymerases by propensity
        random.seed(13)
        self.assertEqual(polymer._choose_polymerase(), self.pol3)
        self.assertEqual(polymer._choose_polymerase(), self.pol1)
        self.assertEqual(polymer._choose_polymerase(), self.pol1)

    def test_move_polymerase(self):
        # Start with clean polymer
        self.setUp()
        polymer = self.construct_polymer()
        # Shift mask back to expose promoter
        self.shift_mask_n(polymer, 10)
        polymer.bind_polymerase(self.pol1, "p1")
        # Make sure promoter is recorded as covered
        self.assertEqual(polymer.uncovered["p1"], 0)
        self.assertTrue(polymer.elements[0].is_covered())
        # Move polymerase 10 spaces and re-expose promoter
        self.move_polymerase_n(polymer, self.pol1, 11)
        self.assertEqual(polymer.uncovered["p1"], 1)
        self.assertFalse(polymer.elements[0].is_covered())
        # Make sure that the mask has also moved appropriately
        self.assertEqual(polymer.mask.start, self.pol1.stop + 1)

        # Check collisions between polymerases
        polymer.bind_polymerase(self.pol2, "p1")
        self.move_polymerase_n(polymer, self.pol2, 15)
        # One polymerase should not be able to pass another
        self.assertEqual(self.pol2.stop + 1, self.pol1.start)
        self.assertFalse(polymer.elements_intersect(self.pol1, self.pol2))

        # Remove pol1
        polymer.terminate(self.pol1, "gene1")
        # Make sure that pol2 is unable to shift mask back
        old_mask_start = polymer.mask.start
        self.move_polymerase_n(polymer, self.pol2, 20)
        self.assertEqual(old_mask_start, polymer.mask.start)
        self.assertEqual(polymer.mask.start, self.pol2.stop + 1)

        # Move back mask to expose terminator
        self.shift_mask_n(polymer, 60)
        # Terminator should be uncovered
        self.assertFalse(polymer.elements[1].is_covered())
        # pol2 should still be on polymer
        self.assertTrue(self.pol2 in polymer.polymerases)

        # Test termination, pol should detach as soon as it hits terminator
        self.move_polymerase_n(polymer, self.pol2, 25)
        # Although removed pol should have intersected with terminator just
        # before termination
        self.assertTrue(
            polymer.elements_intersect(self.pol2, polymer.elements[1])
        )
        # Make sure polymerase has been removed from internal list
        with self.assertRaises(ValueError):
            polymer.polymerases.index(self.pol2)

        # Now check for readthrough
        polymer.bind_polymerase(self.pol1, "p1")
        random.seed(19)
        self.move_polymerase_n(polymer, self.pol1, 42)
        self.assertTrue(polymer.elements[1].readthrough)
        # Move and remove polymerase
        self.move_polymerase_n(polymer, self.pol1, 20)
        polymer.terminate(self.pol1, "gene1")
        # Run polymerase through terminator again, this time it should terminate
        polymer.bind_polymerase(self.pol1, "p1")
        self.move_polymerase_n(polymer, self.pol1, 36)
        self.assertFalse(polymer.elements[1].readthrough)
        with self.assertRaises(ValueError):
            polymer.polymerases.index(self.pol1)
        # Run the polymerase off the end of the polymer
        polymer.bind_polymerase(self.pol1, "p1")
        random.seed(19)
        self.shift_mask_n(polymer, 1000)
        self.move_polymerase_n(polymer, self.pol1, 42)
        self.assertTrue(polymer.elements[1].readthrough)
        self.move_polymerase_n(polymer,
                               self.pol1,
                               polymer.stop - self.pol1.stop + 1)
        self.assertFalse(polymer.elements[1].readthrough)
        with self.assertRaises(ValueError):
            polymer.polymerases.index(self.pol1)

    def test_uncover_elements(self):
        pro1, pro2, pro3, ter, polymer = self.construct_polymer_multi_promoter()
        # Define and bind polymerase
        pol = feature.Polymerase("ecolipol", 10, 30)
        polymer.bind_polymerase(pol, "p1")
        # Move polymerase manually
        pol.start = 12
        pol.stop = 12 + pol.footprint
        # Manually cover elements that should be covered
        pro2.cover()
        pro3.cover()
        # Test uncovered elements
        polymer._uncover_elements(pol)
        self.assertFalse(pro1.is_covered())
        self.assertFalse(pro2.is_covered())
        self.assertFalse(pro3.is_covered())
        # Move polymerase manually again
        pol.start = 16
        pol.stop = 16 + pol.footprint
        pro1.cover()
        pro2.cover()
        pro3.cover()
        # Test for uncovered elements
        polymer._uncover_elements(pol)
        self.assertTrue(pro1.is_covered())
        self.assertFalse(pro2.is_covered())
        self.assertFalse(pro3.is_covered())

    def test_recover_elements(self):
        # Set up
        pro1, pro2, pro3, ter, polymer = self.construct_polymer_multi_promoter()
        # Sanity check
        self.assertFalse(pro1.is_covered())
        self.assertFalse(pro2.is_covered())
        self.assertFalse(pro3.is_covered())
        # Add polymerase
        pol = feature.Polymerase("ecolipol", 10, 30)
        polymer.bind_polymerase(pol, "p1")
        # Uncover promoter after binding
        pro1.uncover()
        # Manually move polymerase
        pol.start = 12
        pol.stop = 12 + pol.footprint
        # Make sure everything is recovered appropriately
        polymer._recover_elements(pol)
        self.assertTrue(pro1.is_covered())
        self.assertTrue(pro2.is_covered())
        self.assertTrue(pro3.is_covered())
        # Uncover everything again
        pro1.uncover()
        pro2.uncover()
        pro3.uncover()
        # Manually move polymerase
        pol.start = 16
        pol.stop = 16 + pol.footprint
        # Again make sure that everything is covered appropriately and that
        # left_most_element is updated
        old_left_most_element = pol.left_most_element
        polymer._recover_elements(pol)
        self.assertFalse(pro1.is_covered())
        self.assertTrue(pro2.is_covered())
        self.assertTrue(pro3.is_covered())
        self.assertNotEqual(old_left_most_element, pol.left_most_element)
        self.assertEqual(pol.left_most_element, 1)
        # Test that everything is uncovered upon termination
        pro2.uncover()
        pro3.uncover()
        pol.start = 21
        pol.stop = 21 + pol.footprint
        polymer._recover_elements(pol)
        self.assertFalse(pro1.is_covered())
        self.assertFalse(pro2.is_covered())
        self.assertFalse(pro3.is_covered())

    def test_resolve_termination(self):
        self.setUp()
        polymer = self.construct_polymer()
        random.seed(22)
        # set up terminator; we don't actually have to place it on the polymer
        term = feature.Terminator("t1", 23, 60,
                                  {"rnapol": {"efficiency": 1.0},
                                   "ecolipol": {"efficiency": 0.6}
                                   }
                                  )
        self.pol1.start = 1
        self.pol1.stop = 23
        term.gene = "gene1"
        term.reading_frame = 1
        self.pol1.reading_frame = 1
        # Create temp termination signal
        self.pol2.release_signal.connect(lambda x: self.assertEqual(x, 60))
        # Terminator will enter readthrough
        term.readthrough = False
        polymer._resolve_termination(self.pol1, term)
        polymer._resolve_termination(self.pol1, term)
        self.assertTrue(term.readthrough)
        # Polymerase and terminator are in different reading frames
        term.readthrough = False
        self.pol2.start = 1
        self.pol2.stop = 23
        self.pol2.reading_frame = 2
        result = polymer._resolve_termination(self.pol2, term)
        self.assertFalse(result)

    def test_resolve_mask_collisions(self):
        self.setUp()
        polymer = self.construct_polymer()
        # Expose promoter
        self.shift_mask_n(polymer, 10)
        # This polymerase is capable of shifting mask backwards
        polymer.bind_polymerase(self.pol1, "p1")
        self.pol1.start += 6
        self.pol1.stop += 6
        old_mask_start = polymer.mask.start
        self.assertFalse(polymer._resolve_mask_collisions(self.pol1))
        self.assertEqual(polymer.mask.start, old_mask_start + 1)
        # The polymerase and mask should never overlap by more than 1 position
        self.pol1.start += 5
        self.pol1.stop += 5
        with self.assertRaises(RuntimeError):
            polymer._resolve_mask_collisions(self.pol1)
        # Let's do this again with a polymerase that cant shift the mask
        self.setUp()
        polymer = self.construct_polymer()
        # Expose promoter
        self.shift_mask_n(polymer, 10)
        polymer.bind_polymerase(self.pol2, "p1")
        old_mask_start = polymer.mask.start
        old_start = self.pol2.start
        old_stop = self.pol2.stop
        # Move pol2 6 spaces
        self.pol2.start += 6
        self.pol2.stop += 6
        self.assertTrue(polymer._resolve_mask_collisions(self.pol2))
        # pol2 should get bumped back to start of mask
        self.assertEqual(self.pol2.start, old_start + 5)
        self.assertEqual(self.pol2.stop, old_stop + 5)
        self.assertEqual(self.pol2.stop, polymer.mask.start - 1)
        self.assertEqual(polymer.mask.start, old_mask_start)

    def test_resolve_collisions(self):
        self.setUp()
        polymer = self.construct_polymer()
        self.shift_mask_n(polymer, 10)
        # Arrange polymerases
        polymer.bind_polymerase(self.pol1, "p1")
        self.move_polymerase_n(polymer, self.pol1, 20)
        polymer.bind_polymerase(self.pol2, "p1")
        # This should result in collision
        for i in range(11):
            self.pol2.move()
        # Record old start position
        old_start = self.pol2.start
        self.assertTrue(polymer._resolve_collisions(self.pol2))
        # pol2 should be bumped back
        self.assertEqual(old_start - 1, self.pol2.start)
        # Make sure polymerases can never overlap
        for i in range(3):
            self.pol2.move()
        with self.assertRaises(RuntimeError):
            polymer._resolve_collisions(self.pol2)

    def test_elements_intersect(self):
        polymer = self.construct_polymer()
        element1 = feature.Promoter("promoter1",
                                    5,
                                    15,
                                    ["ecolipol", "rnapol"]
                                    )
        element2 = feature.Promoter("promoter1",
                                    5,
                                    15,
                                    ["ecolipol", "rnapol"]
                                    )
        self.assertTrue(polymer.elements_intersect(element1, element2))
        self.assertTrue(polymer.elements_intersect(element2, element1))

        element1.start = 15
        element1.stop = 25
        self.assertTrue(polymer.elements_intersect(element1, element2))
        self.assertTrue(polymer.elements_intersect(element2, element1))

        element1.start = 16
        self.assertFalse(polymer.elements_intersect(element1, element2))
        self.assertFalse(polymer.elements_intersect(element2, element1))

        element1.start = 3
        element1.stop = 17
        self.assertTrue(polymer.elements_intersect(element1, element2))
        self.assertTrue(polymer.elements_intersect(element2, element1))

        element1.start = 7
        element1.stop = 10
        self.assertTrue(polymer.elements_intersect(element1, element2))
        self.assertTrue(polymer.elements_intersect(element2, element1))


class TestGenomeMethods(TestPolymerMethods):

    def setUp(self):
        super().setUp()
        # Set up transcript template
        self.transcript_template = [{'type': 'transcript',
                                     'name': 'rnapol',
                                     'length': 200,
                                     'start': 10,
                                     'stop': 210,
                                     'rbs': -15},
                                    {'type': 'transcript',
                                     'name': 'proteinX',
                                     'length': 40,
                                     'start': 230,
                                     'stop': 270,
                                     'rbs': -15},
                                    {'type': 'transcript',
                                     'name': 'proteinY',
                                     'length': 300,
                                     'start': 300,
                                     'stop': 600,
                                     'rbs': -15}]
        # Additional attributes for recording signal firings
        self.fired_transcript = False

    def construct_polymer(self):
        self.mask.stop = 700
        poly = polymer.Genome("g1", 700, [self.promoter, self.terminator],
                              self.transcript_template, self.mask)
        poly.transcript_signal.connect(self.fire_transcript)
        poly.termination_signal.connect(self.fire_termination)
        return poly

    def fire_transcript(self, transcript):
        self.fired_transcript = True
        self.transcript = transcript

    def test_bind_polymerase(self):
        self.setUp()
        polymer = self.construct_polymer()
        with self.assertRaises(RuntimeError):
            polymer.bind_polymerase(self.pol1, "p1")

        self.shift_mask_n(polymer, 20)

        self.assertFalse(self.fired_transcript)
        polymer.bind_polymerase(self.pol1, "p1")
        self.assertTrue(self.fired_transcript)

        self.assertEqual(self.transcript.stop, polymer.stop)
        self.assertTrue(
            self.transcript.shift_mask in self.pol1.move_signal._handlers
        )

        # Test that mask in transcript gets uncovered appropriately as
        # polymerase moves
        old_mask_start = self.transcript.mask.start
        self.move_polymerase_n(polymer, self.pol1, 10)

        self.assertEqual(self.transcript.mask.start, old_mask_start + 10)

    def test_build_transcript(self):
        # TODO: update to test that overlapping genes are ordered correctly
        self.setUp()
        polymer = self.construct_polymer()

        with self.assertRaises(AssertionError):
            polymer._build_transcript(0, 0)

        with self.assertRaises(RuntimeError):
            polymer._build_transcript(600, 600)

        transcript = polymer._build_transcript(0, 230)
        self.assertEqual(len(transcript.elements), 2)

        transcript = polymer._build_transcript(0, 300)
        self.assertEqual(len(transcript.elements), 4)

        transcript = polymer._build_transcript(0, 600)
        self.assertEqual(len(transcript.elements), 6)

        transcript = polymer._build_transcript(200, 600)
        self.assertEqual(len(transcript.elements), 4)

        self.assertEqual(transcript.mask.start, 200)
        self.assertEqual(transcript.mask.stop, 600)

        # Check positions of RBSs and tranlation stop sites
        self.assertEqual(transcript.elements[0].start, 215)
        self.assertEqual(transcript.elements[0].stop, 230)
        self.assertEqual(transcript.elements[1].start, 269)
        self.assertEqual(transcript.elements[1].stop, 270)
        self.assertEqual(transcript.elements[2].start, 285)
        self.assertEqual(transcript.elements[2].stop, 300)
        self.assertEqual(transcript.elements[3].start, 599)
        self.assertEqual(transcript.elements[3].stop, 600)


class TestTranscriptMethods(unittest.TestCase):
    """
    No transcript methods to test right now...
    """
    pass
