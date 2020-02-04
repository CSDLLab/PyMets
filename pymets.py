import argparse
import sys,os,platform
from pymets.io.read_swc import read_swc_trees
from pymets.io.read_json import read_json
from pymets.metric.diadem_metric import diadem_metric
from pymets.metric.length_metric import length_metric

metric_list = [
    "diadem_metric",
    "overall_length",
    "matched_length",
    "DM","OL","ML"
]

def read_parameters():
    parser = argparse.ArgumentParser(
        description="pymet 1.0"
    )

    parser.add_argument(
        "--test",
        "-T",
        help="the route of the test file",
        required=True,
        nargs='*',
    )
    parser.add_argument(
        "--gold",
        "-G",
        help="the route of the gold file",
        required=True
    )
    parser.add_argument(
        "--metric",
        "-M",
        help="choose a metric method",
        required=True
    )
    parser.add_argument(
        "--output",
        "-O",
        help="the route of the output file.\nif not specified, output to screen",
        required=False
    )
    parser.add_argument(
        "--config",
        "-C",
        help="special config for different metric method",
        required=False
    )
    parser.add_argument(
        "--reverse",
        "-R",
        help="output the answer when we switch the gold and test tree",
        required=False
    )
    return parser.parse_args()

def pymets(DEBUG=True):
    abs_dir = os.path.abspath("")

    sys.path.append(abs_dir)
    sys.path.append(os.path.join(abs_dir, "src"))
    sys.path.append(os.path.join(abs_dir, "test"))

    args = read_parameters()

    test_swc_files = [os.path.join(abs_dir, path) for path in args.test]
    gold_swc_file = os.path.join(abs_dir, args.gold)

    reverse = args.reverse
    if reverse is None:
        reverse = False
    metric = args.metric
    if metric not in metric_list:
        raise Exception("[Error: ] Unknown metric method {}".format(
            metric
        ))
    output_dest = args.output
    if output_dest is not None:
        output_dest = os.path.join(abs_dir, output_dest)
    config = args.config
    if config is None:
        if platform.system() == "Windows":
            if metric == "diadem_metric" or metric == "DM":
                config = os.path.join(abs_dir, "config\\diadem_metric.json")
            if metric in ["overall_length", "matched_length", "OL", "ML"]:
                config = os.path.join(abs_dir, "config\\length_metric.json")
        elif platform.system() == "Linux":
            if metric == "diadem_metric" or metric == "DM":
                config = os.path.join(abs_dir, "config/diadem_metric.json")
            if metric in ["overall_length", "matched_length", "OL", "ML"]:
                config = os.path.join(abs_dir, "config/length_metric.json")

    if DEBUG:
        print("Config = {}".format(config))

    test_swc_trees = []
    for test_swc_file in test_swc_files:
        test_swc_trees += read_swc_trees(test_swc_file)
    gold_swc_trees = read_swc_trees(gold_swc_file)
    config = read_json(config)

    print("There are {} test image(s) and {} gold image(s)".format(len(test_swc_trees), len(gold_swc_trees)))
    if len(gold_swc_trees) == 0:
        raise Exception("[Error:  ] No gold image detected")
    if len(gold_swc_trees) > 1:
        print("[Warning:  ] More than one gold image detected, only the first one will be used")

    gold_swc_treeroot = gold_swc_trees[0]
    for test_swc_treeroot in test_swc_trees:
        if metric == "diadem_metric" or metric == "DM":
            diadem_metric(swc_test_tree=test_swc_treeroot,
                          swc_gold_tree=gold_swc_treeroot,
                          config=read_json("D:\gitProject\mine\PyMets\config\diadem_metric.json"))
            if reverse:
                diadem_metric(swc_test_tree=gold_swc_treeroot,
                              swc_gold_tree=test_swc_treeroot,
                              config=read_json("D:\gitProject\mine\PyMets\config\diadem_metric.json"))
        if metric == "overall_length" or metric == "OL":
            config["method"] = 1
            length_metric(gold_swc_treeroot, test_swc_treeroot,
                          abs_dir, config)
            if reverse:
                length_metric(test_swc_treeroot, gold_swc_treeroot ,
                              abs_dir, config)
        if metric == "matched_length" or metric == "ML":
            config["method"] = 2
            length_metric(gold_swc_treeroot, test_swc_treeroot,
                          abs_dir, config)
            if reverse:
                length_metric(test_swc_treeroot, gold_swc_treeroot ,
                              abs_dir, config)

if __name__ == "__main__":
    pymets()
# python ./pymets.py --test D:\gitProject\mine\PyMets\test\data_example\test\30_18_10_test.swc --gold D:\gitProject\mine\PyMets\test\data_example\gold\30_18_10_gold.swc --metric matched_length
# python ./pymets.py --test D:\gitProject\mine\PyMets\test\data_example\test\ExampleTest.swc --gold D:\gitProject\mine\PyMets\test\data_example\gold\ExampleGoldStandard.swc --metric diadem_metric
