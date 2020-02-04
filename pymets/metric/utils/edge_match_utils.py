from pymets.model.euclidean_point import EuclideanPoint,Line
from pymets.metric.utils.config_utils import DINF
from pymets.model.swc_node import add_mid_node,SwcNode,get_bounds,get_lca

import os
import time
import numpy as np
import math
from anytree import PreOrderIter
from rtree import index

MIN_SIZE = 0.8


# get a dict, id_edge_dict[id] = tuple(swc_node1, swc_node2)
def get_idedge_dict(swc_tree=None):
    id_edge_dict = {}
    swc_tree_list = [node for node in PreOrderIter(swc_tree.root())]
    for node in swc_tree_list:
        if node.is_virtual() or node.parent.is_virtual():
            continue
        id_edge_dict[node.get_id()] = tuple([node, node.parent])
    return id_edge_dict


# construct rtree
def get_edge_rtree(swc_tree=None):
    swc_tree_list = [node for node in PreOrderIter(swc_tree.root())]
    p = index.Property()
    p.dimension = 3
    idx3d = index.Index(properties=p)
    for node in swc_tree_list:
        if node.is_virtual() or node.parent.is_virtual():
            continue

        idx3d.insert(node.get_id(), get_bounds(node, node.parent, extra = node.radius()))
    return idx3d


# find the closest edge base on rtree
def get_nearest_edge_fast(idx3d, point, id_edge_dict, DEBUG=False):
    point_box = (point._pos[0] - MIN_SIZE, point._pos[1] - MIN_SIZE, point._pos[2] - MIN_SIZE,
                 point._pos[0] + MIN_SIZE, point._pos[1] + MIN_SIZE, point._pos[2] + MIN_SIZE)
    hits = list(idx3d.intersection(point_box))

    confirm_list = []
    for h in hits:
        line_tuple = id_edge_dict[h]
        if DEBUG:
            print("\npoint = {}, line_a = {}, line_b = {}".format(
                point._pos, line_tuple[0]._pos, line_tuple[1]._pos)
            )
        e_point = EuclideanPoint(center=point._pos)
        d = e_point.distance(Line(swc_node_1=line_tuple[0], swc_node_2=line_tuple[1]))
        if d < point.radius():
            confirm_list.append(line_tuple)

    return confirm_list


# find successful matched edge
def get_match_edges_e_fast(gold_swc_tree=None, test_swc_tree=None,
                           dis_threshold=0.1, detail_path=None, DEBUG=False):
    match_edge = set()
    used_edge = set()
    idx3d = get_edge_rtree(test_swc_tree)
    id_edge_dict = get_idedge_dict(test_swc_tree)
    gold_node_list = [node for node in PreOrderIter(gold_swc_tree.root())]

    if detail_path is not None:
        with open(detail_path, 'a') as f:
            f.write("# [Detail: ]List of unmatched edges")

    for node in gold_node_list:
        if node.is_virtual() or node.parent.is_virtual():
            continue

        e_node = EuclideanPoint(node._pos)
        e_parent = EuclideanPoint(node.parent._pos)

        conf_list_a = get_nearest_edge_fast(idx3d, node, id_edge_dict)
        conf_list_b = get_nearest_edge_fast(idx3d, node.parent, id_edge_dict)

        test_length = 0.0
        gold_length = 0.0

        done = False
        for line_tuple_a in conf_list_a:
            if done:
                break
            for line_tuple_b in conf_list_b:
                if done:
                    break
                if DEBUG:
                    print("\nnode = {} {}\nnode.parent = {} {}\nnode_line = {},{}\nnode_p_line = {},{}\n".format(
                        node.get_id(), node._pos,
                        node.parent.get_id(), node.parent._pos,
                        line_tuple_a[0]._pos, line_tuple_a[1]._pos,
                        line_tuple_b[0]._pos, line_tuple_b[1]._pos,
                    ))

                test_length = get_lca_length(test_swc_tree, \
                                             line_tuple_a, \
                                             line_tuple_b, \
                                             Line([e_node._pos, e_parent._pos]))
                gold_length = node.parent_distance()
                if math.fabs(test_length - gold_length) > gold_length/10:
                    continue

                if line_tuple_a[0].get_id() == 6:
                    print("!!")

                if not LCA_clean(idx3d, id_edge_dict, test_swc_tree, line_tuple_a, line_tuple_b, Line([e_node._pos, e_parent._pos]), used_edge):
                    continue
                match_edge.add(tuple([node, node.parent]))
                done = True

        if not done and detail_path is not None:
            localtime = time.asctime(time.localtime(time.time()))
            with open(detail_path, 'a') as f:
                # f.write("\n------------------------------------\n")
                # f.write(localtime)
                # f.write("\nedge:\npoint_a = {} {}\npoint_b = {} {}\n".format(
                #     node.get_id(), node._pos,
                #     node.parent.get_id(),node.parent._pos
                # ))
                f.write("{} {} {} {} {} {} {}\n".format(
                    node.get_id(), 10, node.get_x(), node.get_y(), node.get_z(), node.radius(), node.parent.get_id())
                )
                f.write("{} {} {} {} {} {} {}\n".format(
                    node.parent.get_id(), 10, node.parent.get_x(), node.parent.get_y(), node.parent.get_z(), node.parent.radius(), node.parent.parent.get_id())
                )
    if detail_path is not None:
        with open(detail_path, 'a') as f:
            f.write("# -----END-----")
    return match_edge


# get all node from current node to the LCA node
def get_route_node(current_node, lca_id):
    res_list = []
    while not current_node.is_virtual() and not current_node.get_id() == lca_id:
        res_list.append(current_node)
        current_node = current_node.parent
    # if current_node.is_virtual():
    #
    #     raise Exception("[Error: ] something wrong in LCA process")
    res_list.append(current_node)
    return res_list


def is_same_poi(node_a, node_b):
    if math.fabs(node_a.get_x() - node_b.get_x()) < 1e-8 and \
            math.fabs(node_a.get_x() - node_b.get_x()) < 1e-8 and \
            math.fabs(node_a.get_x() - node_b.get_x()) < 1e-8:
        return True
    return False


# get the distance of two matched closest edges
def get_lca_length(gold_swc_tree, gold_line_tuple_a, gold_line_tuple_b, test_line):
    point_a, point_b = test_line.get_points()
    gold_line_a = Line(coords=[gold_line_tuple_a[0]._pos, gold_line_tuple_a[1]._pos])
    gold_line_b = Line(coords=[gold_line_tuple_b[0]._pos, gold_line_tuple_b[1]._pos])

    foot_a = point_a.get_foot_point(gold_line_a)
    foot_b = point_b.get_foot_point(gold_line_b)

    if gold_line_tuple_a[0].get_id() == gold_line_tuple_b[0].get_id() and \
        gold_line_tuple_a[1].get_id() == gold_line_tuple_b[1].get_id():
            return foot_a.distance(foot_b)

    lca_id = get_lca(gold_line_tuple_a[0], gold_line_tuple_b[0])

    route_list_a = get_route_node(gold_line_tuple_a[0], lca_id)
    route_list_b = get_route_node(gold_line_tuple_b[0], lca_id)

    lca_length = 0.0
    if gold_line_tuple_a[1] in route_list_a:
        route_list_a.remove(gold_line_tuple_a[0])
        lca_length += foot_a.distance(EuclideanPoint(gold_line_tuple_a[1]._pos))
    else:
        lca_length += foot_a.distance(EuclideanPoint(gold_line_tuple_a[0]._pos))
    route_list_a.pop()
    if gold_line_tuple_b[1] in route_list_b:
        route_list_b.remove(gold_line_tuple_b[0])
        lca_length += foot_b.distance(EuclideanPoint(gold_line_tuple_b[1]._pos))
    else:
        lca_length += foot_b.distance(EuclideanPoint(gold_line_tuple_b[0]._pos))
    route_list_b.pop()

    route_list = route_list_a + route_list_b
    for node in route_list:
        lca_length += node.parent_distance()
    return lca_length


def LCA_clean(idx3d, id_edge_dict, gold_swc_tree, gold_line_tuple_a, gold_line_tuple_b, test_line, used_edge):
    point_a, point_b = test_line.get_points()
    gold_line_a = Line(coords=[gold_line_tuple_a[0]._pos, gold_line_tuple_a[1]._pos])
    gold_line_b = Line(coords=[gold_line_tuple_b[0]._pos, gold_line_tuple_b[1]._pos])

    foot_a = point_a.get_foot_point(gold_line_a)
    foot_b = point_b.get_foot_point(gold_line_b)
    swc_foot_a = SwcNode(nid=gold_swc_tree._size + 1, ntype=3, center=foot_a._pos,
                         radius=(gold_line_tuple_a[0].radius() + gold_line_tuple_a[1].radius())/2,
                         )
    swc_foot_b = SwcNode(nid=gold_swc_tree._size + 2, ntype=3, center=foot_b._pos,
                         radius=(gold_line_tuple_b[0].radius() + gold_line_tuple_b[1].radius())/2,
                         )
    if gold_line_tuple_a[0].get_id() == gold_line_tuple_b[0].get_id() and \
            gold_line_tuple_a[1].get_id() == gold_line_tuple_b[1].get_id():
        if tuple([gold_line_tuple_a[0].get_id(), gold_line_tuple_a[1].get_id()]) not in used_edge:
            used_edge.add(tuple([gold_line_tuple_a[0].get_id(), gold_line_tuple_a[1].get_id()]))
            dis1 = swc_foot_a.distance(gold_line_tuple_a[1])
            dis2 = swc_foot_b.distance(gold_line_tuple_a[1])
            dis3 = swc_foot_a.distance(gold_line_tuple_a[0])
            dis4 = swc_foot_b.distance(gold_line_tuple_a[0])
            if (dis1 < 1e-8 or dis3 < 1e-8) and (dis2 < 1e-8 or dis4 < 1e-8):
                return True
            if dis1 < 1e-8 or dis3 < 1e-8:
                gold_swc_tree._size += 1
                swc_foot_b._id = gold_swc_tree._size
                add_mid_node(idx3d,id_edge_dict,gold_line_tuple_a[0],gold_line_tuple_a[1],swc_foot_b)
            elif dis2 < 1e-8 or dis4 < 1e-8:
                gold_swc_tree._size += 1
                swc_foot_a._id = gold_swc_tree._size
                add_mid_node(idx3d,id_edge_dict,gold_line_tuple_a[0],gold_line_tuple_a[1],swc_foot_a)
            elif dis1 > dis2:
                gold_swc_tree._size += 1
                swc_foot_a._id = gold_swc_tree._size
                add_mid_node(idx3d,id_edge_dict,gold_line_tuple_a[0],gold_line_tuple_a[1],swc_foot_a)
                gold_swc_tree._size += 1
                swc_foot_b._id = gold_swc_tree._size
                add_mid_node(idx3d,id_edge_dict,swc_foot_a,gold_line_tuple_a[1],swc_foot_b)
            else:
                gold_swc_tree._size += 1
                swc_foot_b._id = gold_swc_tree._size
                add_mid_node(idx3d,id_edge_dict,gold_line_tuple_a[0],gold_line_tuple_a[1],swc_foot_b)
                gold_swc_tree._size += 1
                swc_foot_a._id = gold_swc_tree._size
                add_mid_node(idx3d,id_edge_dict,swc_foot_b,gold_line_tuple_a[1],swc_foot_a)
            return True
        else:
            return False

    lca_id = get_lca(gold_line_tuple_a[0], gold_line_tuple_b[0])

    route_list_a = get_route_node(gold_line_tuple_a[0], lca_id)
    route_list_b = get_route_node(gold_line_tuple_b[0], lca_id)

    route_list_a.pop()
    route_list_b.pop()

    if gold_line_tuple_a[1] in route_list_a and is_same_poi(gold_line_tuple_a[1], foot_a):
        route_list_a.remove(gold_line_tuple_a[0])
    if gold_line_tuple_b[1] in route_list_b and is_same_poi(gold_line_tuple_b[1], foot_b):
        route_list_b.remove(gold_line_tuple_b[0])
    if gold_line_tuple_a[1] not in route_list_a and not is_same_poi(gold_line_tuple_a[0], foot_a):
        route_list_a.append(gold_line_tuple_a[1])
    if gold_line_tuple_b[1] not in route_list_b and not is_same_poi(gold_line_tuple_b[0], foot_b):
        route_list_b.append(gold_line_tuple_b[1])

    route_list = route_list_a + route_list_b

    for node in route_list:
        if tuple([node.get_id(), node.parent.get_id()]) in used_edge:
            return False

    if not is_same_poi(swc_foot_a, gold_line_tuple_a[0]) and not is_same_poi(swc_foot_a,gold_line_tuple_a[1]):
        gold_swc_tree._size += 1
        swc_foot_a._id = gold_swc_tree._size
        add_mid_node(idx3d, id_edge_dict, gold_line_tuple_a[0], gold_line_tuple_a[1], swc_foot_a)
        route_list.append(swc_foot_a)
    if not is_same_poi(swc_foot_b, gold_line_tuple_b[0]) and not is_same_poi(swc_foot_b,gold_line_tuple_b[1]):
        gold_swc_tree._size += 1
        swc_foot_b._id = gold_swc_tree._size
        add_mid_node(idx3d, id_edge_dict, gold_line_tuple_b[0], gold_line_tuple_b[1], swc_foot_b)
        route_list.append(swc_foot_b)

    for node in route_list:
        used_edge.add(tuple([node.get_id(), node.parent.get_id()]))
    return True


if __name__ == "__main__":
    pass

