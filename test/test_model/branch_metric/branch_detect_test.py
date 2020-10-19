import unittest
import sys
from pyneval.metric.branch_leaf_metric import branch_leaf_metric
from pyneval.model.swc_node import SwcTree, SwcNode
from pyneval.io.read_json import read_json


class BranchDetectTest(unittest.TestCase):
    def test_branch1(self):
        sys.setrecursionlimit(1000000)
        gold_swc_tree = SwcTree()
        test_swc_tree = SwcTree()
        test_swc_tree.load("..\\..\\..\\data\\branch_metric_data\\test\\194444.swc")
        gold_swc_tree.load("..\\..\\..\\data\\branch_metric_data\\gold\\194444.swc")
        config = read_json("..\\..\\..\\config\\branch_metric.json")
        config['threshold_dis'] = 5
        config['mode'] = 1

        branch_result, leaf_result = \
            branch_leaf_metric(test_swc_tree=gold_swc_tree, gold_swc_tree=test_swc_tree, config=config)
        print(branch_result[3])
        print(branch_result[4])
        print(branch_result[5])
        self.assertEqual(61, branch_result[3])
        self.assertEqual(73, branch_result[4])
        self.assertEqual(1.6772653444540893, branch_result[5])

    def test_branch2(self):
        sys.setrecursionlimit(1000000)
        gold_swc_tree = SwcTree()
        test_swc_tree = SwcTree()
        test_swc_tree.load("..\\..\\..\\data\\branch_metric_data\\test\\194444.swc")
        gold_swc_tree.load("..\\..\\..\\data\\branch_metric_data\\gold\\194444.swc")
        config = read_json("..\\..\\..\\config\\branch_metric.json")
        config['threshold_dis'] = 10
        branch_result, leaf_result = \
            branch_leaf_metric(test_swc_tree=gold_swc_tree, gold_swc_tree=test_swc_tree, config=config)
        print(branch_result[3])
        print(branch_result[4])
        print(branch_result[5])
        self.assertEqual(branch_result[3], 40)
        self.assertEqual(branch_result[4], 52)
        self.assertEqual(branch_result[5], 2.335216925858577)


if __name__ == '__main__':
    unittest.main()
