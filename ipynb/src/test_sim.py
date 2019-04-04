import unittest
import numpy as np
import numpy.testing as npt
from sklearn.utils import check_random_state
from sim import (
    partition_microbes, partition_metabolites,
    cystic_fibrosis_simulation
)


class TestSim(unittest.TestCase):

    def setUp(self):
        pass

    def test_partition_microbes(self):
        microbe_in = np.array([1, 2, 3, 4, 5])
        sigmaQ = 0.01
        state = check_random_state(0)
        num_microbes = 2
        res = partition_microbes(num_microbes, sigmaQ, microbe_in, state)
        npt.assert_allclose(res.sum(axis=1), microbe_in)

    def test_partition_metabolites(self):
        microbe_partition = np.array([[1, 2, 3, 4, 5],
                                      [2, 4, 6, 8, 10]]).T
        metabolite_in = np.array([1, 3, 5, 7, 9])
        uU, uV = 0, 0
        sigmaU, sigmaV = 1, 1
        num_metabolites = 3
        latent_dim = 2
        state = check_random_state(0)
        U, V, res = partition_metabolites(
            uU, sigmaU, uV, sigmaV,
            num_metabolites,
            latent_dim, microbe_partition,
            metabolite_in, state)
        npt.assert_allclose(res.sum(axis=1), metabolite_in)


if __name__ == "__main__":
    unittest.main()
