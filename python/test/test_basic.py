import os
import unittest
from contextlib import contextmanager
from pathlib import Path
from opm2.simulators import BlackOilSimulator

@contextmanager
def pushd(path):
    cwd = os.getcwd()
    if not os.path.isdir(path):
        os.makedirs(path)
    os.chdir(path)
    yield
    os.chdir(cwd)


class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # NOTE: Loading the below extension module involves loading a
        #    a shared library like "simulators.cpython-37m-x86_64-linux-gnu.so"
        #    It turns out Python cannot unload this module, see:
        #
        #  https://stackoverflow.com/a/8295590/2173773
        #  https://bugs.python.org/issue34309
        #
        #    This is a problem when we want to create a new instance for each unit
        #    test. For example, when creating the first instance, static variables in
        #    in the shared object are initialized. However, when the
        #    corresponding Python object is later deleted (when the test finishes),
        #    the shared object is not unloaded and its static variables
        #    stays the same. So when a second Python instance is created,
        #    the same address is used for the static variables in the shared library
        #    i.e. the static variables are referring to the same memory location
        #    as for the first object (and they are not reinitialized).
        #
        #    Unfortunatly, this leads to undefined behavior since the C++ code
        #    for flow simulation uses static variable to keep state information
        #    and since it was not built under the assumption that it would
        #    used as a shared library. It was assumed (?) that a flow simulation
        #    was executed from an executable file (not library file) and only
        #    executed once. To execute another simulation, it was assumed that
        #    the executable would be restarted from a controlling program like
        #    the Shell (which would reload and initialize the object into fresh memory).
        #
        #  TODO: Fix the C++ code such that it allows multiple runs whith the same
        #     object file.
        #
        #  NOTE:  The result of the above is that we can only instantiate a
        #    single simulator object during the unit tests.
        #    This is not how the unittest module was supposed to be used. Usually one
        #    would write multiple test_xxx() methods that are independent and
        #    each method receives a new simulator object (also note that the order
        #    in which each test_xxx() method is called by unittest is not defined).
        #    However, as noted above this is not currently possible.
        #
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir = test_dir.parent.joinpath("test_data/SPE1CASE1a")


    def test_all(self):
        with pushd(self.data_dir):
            sim = BlackOilSimulator("SPE1CASE1.DATA")
            sim.step_init()
            sim.step()

            poro = sim.get_porosity()
            self.assertEqual(len(poro), 300, 'length of porosity vector')
            self.assertAlmostEqual(poro[0], 0.3, places=7, msg='value of porosity')
            poro = poro *.95
            sim.set_porosity(poro)
            sim.step()
            poro2 = sim.get_porosity()
            self.assertAlmostEqual(poro2[0], 0.285, places=7, msg='value of porosity 2')

