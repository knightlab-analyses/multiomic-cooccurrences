import unittest
import pandas as pd
import numpy as np
from evaluate import _edge_roc_curve
import pandas.util.testing as pdt


class TestEvaluate(unittest.TestCase):
    def setUp(self):
        pass

    def test_top_rank_accuracy(self):
        pass


class TestEvaluateBiofilm(unittest.TestCase):

    def setUp(self):
        self.ranks = pd.DataFrame(
            {
                'A': [1, -4, 4, 5],
                'B': [4, -2, 3, 1],
                'C': [4, -1, -3, 2]
            },
            index=['a', 'b', 'c', 'd']
        )

        self.edges = pd.DataFrame(
            {
                'microbe': ['A', 'A', 'B', 'C', 'B'],
                'metabolite': ['a', 'b', 'a', 'd', 'd'],
                'direction': [0, 'I', 'R', 'C', 'C']
            }
        )

    def test_edge_roc_curve(self):

        resR, resC, resI = _edge_roc_curve(self.ranks.T, self.edges, k_max=2)

        expR = pd.DataFrame(
            [[0, 5, 5, 1],
             [0, 5, 5, 1]],
            columns=['FN', 'FP', 'TN', 'TP'],
            index=[1, 2]
        )
        pdt.assert_frame_equal(resR, expR)

        expC = pd.DataFrame(
            [[1, 5, 5, 1],
             [1, 5, 5, 1]],
            columns=['FN', 'FP', 'TN', 'TP'],
            index=[1, 2]
        )
        pdt.assert_frame_equal(resC, expC)

        expI = pd.DataFrame(
            [[1, 6, 5, 0],
             [1, 6, 5, 0]],
            columns=['FN', 'FP', 'TN', 'TP'],
            index=[1, 2]
        )
        pdt.assert_frame_equal(resI, expI)


if __name__ == "__main__":
    unittest.main()
