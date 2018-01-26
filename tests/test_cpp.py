# Test simulation
import unittest
import subprocess


class MainTest(unittest.TestCase):

    def test_cpp(self):
        print("\n\nRunning C++ tests...\n")
        subprocess.call(['pinetree_test'])
