import sys

from pyneval.model.swc_node import SwcTree
from pyneval.metric.utils.km_utils import KM, get_dis_graph
from pyneval.io.read_json import read_json
from pyneval.metric.utils.point_match_utils import get_swc2swc_dicts
from pyneval.io.read_config import *
import jsonschema


def debug_out_list(swc_list, _str):
    print("[debug_out_list ]" + _str + str(len(swc_list)))
    for node in swc_list:
        print(node.get_id(), end=" ")
    print("")


def get_result(test_len, gold_len, switch, km, threshold_dis):
    false_pos_num, false_neg_num, true_pos_num = 0, 0, 0
    # count numer of nodes which are matched, calculate FP, TN, TP
    for i in range(0, gold_len):
        if km.match[i] != -1 and km.G[km.match[i]][i] != -0x3f3f3f3f / 2:
            true_pos_num += 1
    false_neg_num = gold_len - true_pos_num
    false_pos_num = test_len - true_pos_num

    # definition of swich is in function "get_dis_graph"
    if switch:
        false_neg_num, false_pos_num = false_pos_num, false_neg_num

    if true_pos_num != 0:
        mean_dis = -km.get_max_dis() / true_pos_num
    else:
        mean_dis = 0.0
    if mean_dis == -0.0:
        mean_dis = 0.0

    pt_cost = -km.get_max_dis() + threshold_dis * (false_neg_num + false_pos_num) / (
                false_neg_num + false_pos_num + true_pos_num)

    # debug:
    # print("output")
    # print(false_pos_num)
    # print(true_neg_num)
    # print(mean_dis)
    # print(pt_cost)
    return true_pos_num, false_neg_num, false_pos_num, mean_dis, mean_dis * true_pos_num, pt_cost


def get_colored_tree(test_node_list, gold_node_list, switch, km, color):
    '''
    color[0] = tp's color
    color[1] = fp's color
    color[2] = fn's color
    '''
    tp_set = set()
    for i in range(0, len(gold_node_list)):
        if km.match[i] != -1 and km.G[km.match[i]][i] != -0x3f3f3f3f / 2:
            gold_node_list[i]._type = color[0]
            test_node_list[km.match[i]]._type = color[0]
            tp_set.add(test_node_list[km.match[i]])
        else:
            if switch:
                gold_node_list[i]._type = color[1]
            else:
                gold_node_list[i]._type = color[2]
    for node in test_node_list:
        if node not in tp_set:
            if switch:
                node._type = color[2]
            else:
                node._type = color[1]


def score_point_distance(gold_tree, test_tree, test_node_list, gold_node_list,
                         test_gold_dict, threshold_dis, color, metric_mode):
    # disgraph is a 2D ndarray store the distance of nodes in gold and test
    # test_node_list contains only branch or leaf nodes
    dis_graph, switch, test_len, gold_len = get_dis_graph(gold_tree=gold_tree,
                                                          test_tree=test_tree,
                                                          test_node_list=test_node_list,
                                                          gold_node_list=gold_node_list,
                                                          test_gold_dict=test_gold_dict,
                                                          threshold_dis=threshold_dis,
                                                          metric_mode=metric_mode)
    # create a KM object and calculate the minimum match
    km = KM(maxn=max(test_len, gold_len)+10, nx=test_len, ny=gold_len, G=dis_graph)
    km.solve()
    # calculate the result
    true_pos_num, false_neg_num, false_pos_num, \
    mean_dis, tot_dis, pt_cost = get_result(test_len=test_len,
                                            gold_len=gold_len,
                                            switch=switch,
                                            km=km,
                                            threshold_dis=threshold_dis)
    # calculate the number of isolated nodes
    iso_node_num = 0
    for node in test_tree.get_node_list():
        if node.is_isolated():
            iso_node_num += 1
    # get a colored tree with
    get_colored_tree(test_node_list=test_node_list, gold_node_list=gold_node_list,
                     switch=switch, km=km, color=color)
    return gold_len, test_len, true_pos_num, false_neg_num, false_pos_num, mean_dis, tot_dis, pt_cost, iso_node_num


def branch_leaf_metric(gold_swc_tree, test_swc_tree, config):
    """
    calculate branch metric value of two swc trees
    :param gold_swc_tree(Swc Tree) gold standard tree
    :param test_swc_tree(Swc Tree) reconstructed Swc Tree object
    :param config(dict) dict read by json object

    :return branch_result(tuple) a tuple of 9 metric results
    """
    threshold_dis = read_float_config(config=config, config_name="threshold_dis", default=2)
    metric_mode = read_int_config(config=config, config_name="metric_mode", default=1)
    threshold_mode = read_int_config(config=config, config_name="threshold_mode", default=1)

    if threshold_mode == 2:
        # length of the entire gold swc forest
        tot_dis = gold_swc_tree.length()
        # number of edges in the forest
        edge_num = len(gold_swc_tree.get_node_list())-1-len(gold_swc_tree.root().children)
        threshold_dis = threshold_dis * tot_dis / edge_num
    color = [
        read_int_config(config=config, config_name="true_positive_type", default=1),
        read_int_config(config=config, config_name="false_negative_type", default=1),
        read_int_config(config=config, config_name="false_positive_type", default=1)
    ]
    gold_swc_tree.type_clear(0, 0)
    test_swc_tree.type_clear(0, 0)
    test_branch_swc_list = test_swc_tree.get_branch_swc_list()
    gold_branch_swc_list = gold_swc_tree.get_branch_swc_list()

    test_gold_dict = get_swc2swc_dicts(src_node_list=test_swc_tree.get_node_list(),
                                       tar_node_list=gold_swc_tree.get_node_list())

    # debug output
    # debug_out_list(test_branch_swc_list, "test_branch_swc_list")
    # debug_out_list(gold_branch_swc_list, "gold_branch_swc_list")
    # debug_out_list(test_leaf_swc_list, "test_leaf_swc_list")
    # debug_out_list(gold_leaf_swc_list, "gold_leaf_swc_list")

    branch_result_tuple = score_point_distance(gold_tree=gold_swc_tree,
                                               test_tree=test_swc_tree,
                                               test_node_list=test_branch_swc_list,
                                               gold_node_list=gold_branch_swc_list,
                                               test_gold_dict=test_gold_dict,
                                               threshold_dis=threshold_dis,
                                               color=color,
                                               metric_mode=metric_mode)

    res = {
        "gold_len": branch_result_tuple[0],
        "test_len": branch_result_tuple[1],
        "true_pos_num": branch_result_tuple[2],
        "false_neg_num": branch_result_tuple[3],
        "false_pos_num": branch_result_tuple[4],
        "mean_dis": branch_result_tuple[5],
        "tot_dis": branch_result_tuple[6],
        "pt_cost": branch_result_tuple[7],
        "iso_node_num": branch_result_tuple[8]
    }
    return res


if __name__ == "__main__":
    sys.setrecursionlimit(1000000)
    file_name = "fake_data11"
    gold_swc_tree = SwcTree()
    test_swc_tree = SwcTree()

    test_swc_tree.load("../../data/test_data/topo_metric_data/gold_fake_data4.swc")
    gold_swc_tree.load("../../data/test_data/topo_metric_data/test_fake_data4.swc")

    config = read_json("..\\..\\config\\branch_metric.json")
    config_schema = read_json("..\\..\\config\\schemas\\branch_metric_schema.json")
    try:
        jsonschema.validate(config, config_schema)
    except Exception as e:
        raise Exception("[Error: ]Error in analyzing config json file")

    config["metric_mode"] = 2
    config["threshold_dis"] = 1
    config["threshold_mode"] = 2
    branch_result = \
        branch_leaf_metric(test_swc_tree=test_swc_tree, gold_swc_tree=gold_swc_tree, config=config)
    print("---------------Result---------------")
    print("gole_branch_num = {}, test_branch_num = {}\n"
          "true_positive_number  = {}\n"
          "false_negative_num    = {}\n"
          "false_positive_num    = {}\n"
          "matched_mean_distance = {}\n"
          "matched_sum_distance  = {}\n"
          "pt_score              = {}\n"
          "isolated node number  = {}".
          format(branch_result["gold_len"], branch_result["test_len"], branch_result["true_pos_num"],
                 branch_result["false_neg_num"], branch_result["false_pos_num"], branch_result["mean_dis"],
                 branch_result["tot_dis"], branch_result["pt_cost"], branch_result["iso_node_num"]))
    print("----------------End-----------------")
    # with open("../../output/branch_metric/{}_gold.swc".format(file_name), 'w') as f:
    #     f.write(gold_swc_tree.to_str_list())
    # with open("../../output/branch_metric/{}_test.swc".format(file_name), 'w') as f:
    #     f.write(test_swc_tree.to_str_list())

