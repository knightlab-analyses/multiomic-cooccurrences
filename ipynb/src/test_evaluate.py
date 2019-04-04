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
                'microbe': ['A', 'A', 'B', 'C'],
                'metabolite': ['a', 'b', 'a', 'd'],
                'direction': [0, -1, 1, 1]
            }
        )

    def test_edge_roc_curve(self):
        res_pos, res_neg = _edge_roc_curve(self.ranks.T, self.edges, k_max=2)

        exp_pos = pd.DataFrame(
            [[1, 2, 4, 1],
             [0, 4, 4, 2]],
            columns=['FN', 'FP', 'TN', 'TP'],
            index=[1, 2]
        )
        pdt.assert_frame_equal(res_pos, exp_pos)

        exp_neg = pd.DataFrame(
            [[0, 2, 4, 1],
             [0, 5, 4, 1]],
            columns=['FN', 'FP', 'TN', 'TP'],
            index=[1, 2]
        )
        pdt.assert_frame_equal(res_neg, exp_neg)


if __name__ == "__main__":
    unittest.main()
