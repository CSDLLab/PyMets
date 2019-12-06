# bennieHan 2019-11-12 16:01
# all right reserved

from anytree import NodeMixin, iterators, RenderTree
from src.model.euclidean_point import EuclideanPoint,Line
from anytree import PreOrderIter
from rtree import index
import kdtree
import math
import queue
import numpy as np

dis_threshold = 0.1
_3D = "3d"
_2D = "2d"

def Make_Virtual():
    return SwcNode(nid=-1)


def compute_platform_area(r1, r2, h):
    return (r1 + r2) * h * math.pi


# to test
def compute_two_node_area(tn1, tn2, remain_dist):
    """Returns the surface area formed by two nodes
    """
    r1 = tn1.radius()
    r2 = tn2.radius()
    d = tn1.distance(tn2)
    print(remain_dist)

    if remain_dist >= d:
        h = d
    else:
        h = remain_dist
        a = remain_dist / d
        r2 = r1 * (1 - a) + r2 * a

    area = compute_platform_area(r1, r2, h)
    return area


# to test
def compute_surface_area(tn, range_radius):
    area = 0

    # backtrace
    currentDist = 0
    parent = tn.parent
    while parent and currentDist < range_radius:
        remainDist = range_radius - currentDist
        area += compute_two_node_area(tn, parent, remainDist)
        currentDist += tn.distance(parent)
        tn = parent
        parent = tn.parent

    # forwardtrace
    currentDist = 0
    childList = tn.children
    while len(childList) == 1 and currentDist < range_radius:
        child = childList[0]
        remainDist = range_radius - currentDist
        area += compute_two_node_area(tn, child, remainDist)
        currentDist += tn.distance(child)
        tn = child
        childList = tn.children

    return area


class SwcNode(NodeMixin):
    """
        this is a class that temporarily store SWC file
        Attributes:
        id: id of the node,
        type: leaf = 1,continuation = 2, bifurcation = 3,
        parent: pa node id,
        son=[]: son list,
        x: x coordinate,
        y: y coordinate,
        z: z coordinate,
        radius: radius of the node
        surface_area: surface area of the cylinder, radius is current radius, length is the distance to its parent
        volume: volume of the cylinder. radious is the same as above

        parent_trajectory: distance to root
        left_trajectory: distance to the farthest son of left_son
        right_trajectory: distance to the farthest son of right_son

        path_length: distance to parent
        xy_path_length: distance to parent regardless z coordinate
        z_path_lenth: distance to parent
    """

    def __init__(self,
                 nid=-1,
                 ntype=0,
                 radius=1,
                 center=[0, 0, 0],
                 parent=None,
                 depth=0,

                 surface_area=0.0,
                 volume=0.0,

                 parent_trajectory=None,
                 left_trajectory=None,
                 right_trajectory=None,

                 path_length=0.0,
                 xy_path_length=0.0,
                 z_path_lenth=0.0):
        self._id = nid
        self._type = ntype
        self.parent = parent
        self._pos = center
        self._radius = radius
        self.surface_area = surface_area
        self.volume = volume
        self._depth = depth

        self.parent_trajectory=parent_trajectory
        self.left_trajectory=left_trajectory
        self.right_trajectory=right_trajectory

        self.path_length=path_length
        self.xy_path_length=xy_path_length
        self.z_path_length=z_path_lenth

    def add_length(self, swc_node):
        self.path_length += swc_node.path_length
        self.xy_path_length += swc_node.xy_path_length
        self.z_path_length += swc_node.z_path_length

    def add_data(self, swc_node):
        self.path_length += swc_node.path_length
        self.xy_path_length += swc_node.xy_path_length
        self.z_path_length += swc_node.z_path_length
        self.volume += swc_node.volume
        self.surface_area += swc_node.surface_area

    def is_virtual(self):
        """Returns True iff the node is virtual.
        """
        return self._id < 0

    def depth(self):
        return self._depth

    def is_regular(self):
        """Returns True iff the node is NOT virtual.
        """
        return self._id >= 0

    def get_id(self):
        """Returns the ID of the node.
        """
        return self._id

    def distance(self, tn = None, mode = _3D):
        """ Returns the distance to another node.
        It returns 0 if either of the nodes is not regular.

        Args:
          tn : the target node for distance measurement
        """
        if tn is None:
            return 0.0
        if type(tn) == type([]):
            tn = SwcNode(nid=1,center=tn)
        if tn and self.is_regular() and tn.is_regular():
            dx = self._pos[0] - tn._pos[0]
            dy = self._pos[1] - tn._pos[1]
            dz = self._pos[2] - tn._pos[2]
            if mode == _2D:
                dz = 0.0
            d2 = dx * dx + dy * dy + dz * dz

            return math.sqrt(d2)

        return 0.0

    def parent_distance(self):
        """ Returns the distance to it parent.
        """
        return self.distance(self.parent)

    def radius(self):
        return self._radius

    def scale(self, sx, sy, sz, adjusting_radius=True):
        """Transform a node by scaling
        """

        self._pos[0] *= sx
        self._pos[1] *= sy
        self._pos[2] *= sz

        if adjusting_radius:
            self._radius *= math.sqrt(sx * sy)

    def to_swc_str(self):
        return '%d %d %g %g %g %g' % (self._id, self._type, self._pos[0], self._pos[1], self._pos[2], self._radius)

    def get_parent_id(self):
        return -2 if self.is_root else self.parent.get_id()

    def __str__(self):
        return '%d (%d): %s, %g' % (self._id, self._type, str(self._pos), self._radius)


class SwcTree:
    """A class for representing one or more SWC trees.
    For simplicity, we always assume that the root is a virtual node.
    """

    def __init__(self):
        self._root = Make_Virtual()
        self._size = None
        self._total_length = None

        self.depth_array = None
        self.LOG_NODE_NUM=None
        self.lca_parent=None

    def _print(self):
        print(RenderTree(self._root).by_attr("_id"))

    def clear(self):
        self._root = Make_Virtual()

    def is_comment(self, line):
        return line.strip().startswith('#')

    def root(self):
        return self._root

    def regular_root(self):
        return self._root.children

    def node_from_id(self, nid):
        niter = iterators.PreOrderIter(self._root)
        for tn in niter:
            if tn.get_id() == nid:
                return tn
        return None

    def parent_id(self, nid):
        tn = self.node_from_id(nid)
        if tn:
            return tn.get_parent_id()

    def parent_node(self, nid):
        tn = self.node_from_id(nid)
        if tn:
            return tn.parent

    def child_list(self, nid):
        tn = self.node_from_id(nid)
        if tn:
            return tn.children

    def load(self, path):
        self.clear()
        with open(path, 'r') as fp:
            lines = fp.readlines()
            nodeDict = dict()
            for line in lines:
                if not self.is_comment(line):
                    #                     print line
                    data = list(map(float, line.split()))
                    #                     print(data)
                    if len(data) == 7:
                        nid = int(data[0])
                        ntype = int(data[1])
                        pos = data[2:5]
                        radius = data[5]
                        parentId = data[6]
                        tn = SwcNode(nid=nid, ntype=ntype, radius=radius, center=pos)
                        nodeDict[nid] = (tn, parentId)
            fp.close()

            for _, value in nodeDict.items():
                tn = value[0]
                parentId = value[1]
                if parentId == -1:
                    tn.parent = self._root
                    tn._depth = 0
                else:
                    parentNode = nodeDict.get(parentId)
                    if parentNode:
                        tn.parent = parentNode[0]
                        tn._depth = tn.parent._depth+1

    def save(self, path):
        with open(path, 'w') as fp:
            niter = iterators.PreOrderIter(self._root)
            for tn in niter:
                if tn.is_regular():
                    fp.write('%s %d\n' % (tn.to_swc_str(), tn.get_parent_id()))
            fp.close()

    def has_regular_node(self):
        return len(self.regular_root()) > 0

    def node_count(self, regular=True, force_update=False):
        if force_update == False and self._size is not None:
            return self._size

        count = 0
        niter = iterators.PreOrderIter(self._root)
        for tn in niter:
            if regular:
                if tn.is_regular():
                    count += 1
            else:
                count += 1
        self._size = count
        return count

    def parent_distance(self, nid):
        d = 0
        tn = self.node(nid)
        if tn:
            parent_tn = tn.parent
            if parent_tn:
                d = tn.distance(parent_tn)

        return d

    def scale(self, sx, sy, sz, adjusting_radius=True):
        niter = iterators.PreOrderIter(self._root)
        for tn in niter:
            tn.scale(sx, sy, sz, adjusting_radius)

    def length(self, force_update=False):
        if self._total_length is not None and force_update == False:
            return self._total_length

        niter = iterators.PreOrderIter(self._root)
        result = 0
        for tn in niter:
            result += tn.parent_distance()

        return result

    def radius(self, nid):
        return self.node(nid).radius()

    def get_depth_array(self):
        self.depth_array = [0] * (self.node_count()+10)
        for node in PreOrderIter(self.root()):
            self.depth_array[node.get_id()] = node.depth()

    def get_lca_preprocess(self):
        self.get_depth_array()
        self.LOG_NODE_NUM = math.ceil(math.log(self.node_count(), 2))
        self.lca_parent = np.zeros(shape=(self.node_count()+10, self.LOG_NODE_NUM),dtype=int)
        tree_node_list = [node for node in PreOrderIter(self.root())]

        for node in tree_node_list:
            if node.is_virtual():
                continue
            self.lca_parent[node.get_id()][0] = node.parent.get_id()
        for k in range(self.LOG_NODE_NUM - 1):
            for v in range(1,self.node_count()):
                if self.lca_parent[v][k] < 0:
                    self.lca_parent[v][k + 1] = -1
                else:
                    self.lca_parent[v][k + 1] = self.lca_parent[int(self.lca_parent[v][k])][k]
        return True

    def get_lca(self, u, v):
        lca_parent = self.lca_parent
        LOG_NODE_NUM = self.LOG_NODE_NUM
        depth_array = self.depth_array

        if depth_array[u] > depth_array[v]:
            u,v = v,u
        for k in range(LOG_NODE_NUM):
            if depth_array[v] - depth_array[u] >> k & 1:
                v = lca_parent[v][k]
        if u == v:
            return u
        for k in range(LOG_NODE_NUM -1,-1,-1):
            if lca_parent[u][k] != lca_parent[v][k]:
                u = lca_parent[u][k]
                v = lca_parent[v][k]
        return lca_parent[u][0]

    def align_roots(self, gold_tree, mode="average", DEBUG=False):
        offset = EuclideanPoint()
        stack = queue.LifoQueue()
        gold_anchor = np.zeros(3)
        test_anchor = np.zeros(3)
        if mode == "average":
            gold_tree_list = [node for node in PreOrderIter(gold_tree.root())]
            test_tree_list = [node for node in PreOrderIter(self.root())]
            for node in gold_tree_list:
                gold_anchor += np.array(node._pos)
            for node in test_tree_list:
                test_anchor += np.array(node._pos)
            gold_anchor /= len(gold_tree_list) - 1
            test_anchor /= len(test_tree_list) - 1
        elif mode == "root":
            test_anchor = list(self.root().children)[0]
            gold_anchor = list(gold_tree.root().children)[0]

        offset._pos = (gold_anchor - test_anchor).tolist()
        if DEBUG:
            print("off_set:x = {}, y = {}, z = {}".format(offset._pos[0], offset._pos[1], offset._pos[2]))

        stack.put(self.root().children[0])
        while not stack.empty():
            node = stack.get()
            if node.is_virtual():
                continue

            node._pos[0] += offset._pos[0]
            node._pos[1] += offset._pos[1]
            node._pos[2] += offset._pos[2]

            for son in node.children:
                stack.put(son)


def get_default_threshold(gold_swc_tree):
    global dis_threshold
    total_length = gold_swc_tree.length()
    total_node = gold_swc_tree.node_count()
    if total_node <= 1:
        dis_threshold = 0.1
    else:
        dis_threshold = (total_length/total_node)/10


def get_bounds(point_a, point_b):
    point_a = np.array(point_a._pos)
    point_b = np.array(point_b._pos)
    res = np.where(point_a>point_b,point_b,point_a).tolist() + np.where(point_a>point_b,point_a,point_b).tolist()
    return tuple(res)


def get_idedge_dict(swc_tree=None):
    id_edge_dict = {}
    swc_tree_list = [node for node in PreOrderIter(swc_tree.root())]
    for node in swc_tree_list:
        if node.is_virtual() or node.parent.is_virtual():
            continue
        id_edge_dict[node.get_id()] = tuple([node, node.parent])
    return id_edge_dict


def get_edge_rtree(swc_tree=None):
    swc_tree_list = [node for node in PreOrderIter(swc_tree.root())]
    p = index.Property()
    p.dimension = 3
    idx3d = index.Index(properties=p)
    for node in swc_tree_list:
        if node.is_virtual() or node.parent.is_virtual():
            continue
        if node.get_id() == 7:
            print("---")
        idx3d.insert(node.get_id(), get_bounds(node, node.parent))
    return idx3d


def get_nearest_edge(idx3d, point,id_edge_dict):
    nearest_line_id = list(idx3d.nearest(get_bounds(point,point)))[0]
    line_tuple = id_edge_dict[nearest_line_id]
    line = Line(coords=[line_tuple[0]._pos, line_tuple[1]._pos], is_segment=True)
    # print("point = {}, line_a = {}, line_b = {}".format(point._pos, line.coords[0], line.coords[1]))
    dis = point.distance(line)
    return line_tuple, dis


def get_match_edges_e(gold_swc_tree=None, test_swc_tree=None, DEBUG=False):
    match_edge = set()
    idx3d = get_edge_rtree(test_swc_tree)
    id_edge_dict = get_idedge_dict(test_swc_tree)
    gold_node_list = [node for node in PreOrderIter(gold_swc_tree.root())]
    global dis_threshold

    for node in gold_node_list:
        if node.is_virtual() or node.parent.is_virtual():
            continue

        e_node = EuclideanPoint(node._pos)
        e_parent = EuclideanPoint(node.parent._pos)

        line_tuple_a, dis_a = get_nearest_edge(idx3d, e_node, id_edge_dict)
        line_tuple_b, dis_b = get_nearest_edge(idx3d, e_parent, id_edge_dict)

        # get_max_dis(test_swc_tree, line_tuple_a, line_tuple_b, Line([e_node._pos, e_parent._pos]))
        if dis_a <= dis_threshold and dis_b <= dis_threshold:
            match_edge.add(tuple([node,node.parent]))
    return match_edge


def get_kdtree_data(kd_node):
    return kd_node[0].data


def check_match(gold_node_knn, son_node_knn, edge_set, id_center_dict):
    for pa in gold_node_knn:
        test_pa_node = id_center_dict[tuple(get_kdtree_data(pa))]
        for tpn in test_pa_node:
            for sn in son_node_knn:
                test_son_node = id_center_dict[tuple(get_kdtree_data(sn))]
                for tsn in test_son_node:
                    if tuple([tpn, tsn]) in edge_set:
                        edge_set.remove(tuple([tpn, tsn]))
                        edge_set.remove(tuple([tsn, tpn]))
                        return tuple([tpn, tsn])
    return None


def get_match_edges_p(gold_swc_tree=None, test_swc_tree=None, knn=3, DEBUG=False):
    match_edge = {}
    test_swc_list = [node for node in PreOrderIter(test_swc_tree.root())]
    if DEBUG:
        for item in test_swc_list:
            print("---{} {}".format(item.get_id(), item._pos))
    id_center_dict = {}
    center_list = []
    edge_set = set()
    global dis_threshold
    for node in test_swc_list:
        if node.is_virtual():
            continue
        if tuple(node._pos) not in id_center_dict.keys():
            id_center_dict[tuple(node._pos)] = []
        id_center_dict[tuple(node._pos)].append(node)
        center_list.append(node._pos)
        for son in node.children:
            edge_set.add(tuple([node, son]))
            edge_set.add(tuple([son, node]))

    test_kdtree = kdtree.create(center_list)

    stack = queue.LifoQueue()
    stack.put(gold_swc_tree.root())
    while not stack.empty():
        gold_node = stack.get()
        for son in gold_node.children:
            stack.put(son)

        if gold_node.is_virtual():
            continue

        gold_node_knn = test_kdtree.search_knn(gold_node._pos, knn)

        for node in gold_node_knn:
            if gold_node.distance(get_kdtree_data(node)) > dis_threshold:
                gold_node_knn.remove(node)
        if DEBUG:
            print("knn of gold = {}".format(gold_node_knn))
        for son in gold_node.children:
            son_node_knn = test_kdtree.search_knn(son._pos, knn)
            if DEBUG:
                print("son of gold = {}".format(son_node_knn))
            match = check_match(gold_node_knn, son_node_knn, edge_set, id_center_dict)
            if match is not None:
                match_edge[tuple([gold_node, son])] = match

    return match_edge


def get_route_node(current_node, lca_id):
    res_list = []
    while not current_node.is_virtual() and not current_node.get_id() == lca_id:
        res_list.append(current_node)
    if current_node.is_virtual():
        raise Exception("[Error: ] something wrong in LCA process")
    res_list.append(current_node)
    return res_list


def get_max_dis(gold_swc_tree, gold_line_a, gold_line_b, test_line):
    gold_swc_tree.get_lca_preprocess()
    lca_id = gold_swc_tree.get_lca(gold_line_a[0].get_id(), gold_line_b[0].get_id())
    route_list = []
    route_list += get_route_node(gold_line_a[0], lca_id)
    route_list += get_route_node(gold_line_b[0], lca_id)
    # root节点可能有问题
    for node in route_list:
        print(node.get_id())


if __name__ == '__main__':
    print('testing ...')
    tree = SwcTree()
    tree.load("D:\gitProject\mine\PyMets\\test\data_example\gold\gold.swc")
    tree.get_lca_preprocess()
    print(tree.get_lca(2,6))
