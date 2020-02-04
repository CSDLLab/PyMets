# bennieHan 2019-11-12 16:01
# all right reserved

from anytree import NodeMixin, iterators, RenderTree
from pymets.model.euclidean_point import EuclideanPoint
from anytree import PreOrderIter
import math
import queue
import numpy as np

_3D = "3d"
_2D = "2d"


# input node id of two swc_node，calculate LCA
def get_lca(u, v):
    tmp_set = set()
    tmp_u = u
    tmp_v = v
    while u.get_id() != -1:
        tmp_set.add(u.get_id())
        u = u.parent

    while v.get_id() != -1:
        if v.get_id() in tmp_set:
            return v.get_id()
        v = v.parent
    return None


# get bounding box of a segment
def get_bounds(point_a, point_b, extra = 0):
    point_a = np.array(point_a._pos)
    point_b = np.array(point_b._pos)
    res = (np.where(point_a>point_b,point_b,point_a) - extra).tolist() + (np.where(point_a>point_b,point_a,point_b) + extra).tolist()

    return tuple(res)


def add_mid_node(rtree, id_edge_dict, son, parent, new_node):
    rtree.delete(son.get_id(), get_bounds(son, son.parent, extra=son.radius()))

    son.parent = new_node
    new_node.parent = parent

    rtree.insert(son.get_id(), get_bounds(son, son.parent, extra=son.radius()))
    rtree.insert(new_node.get_id(), get_bounds(new_node, new_node.parent, extra=new_node.radius()))

    id_edge_dict[son.get_id()] = tuple([son, son.parent])
    id_edge_dict[new_node.get_id()] = tuple([new_node, new_node.parent])


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

    def get_x(self):
        return self._pos[0]

    def get_y(self):
        return self._pos[1]

    def get_z(self):
        return self._pos[2]

    def set_x(self, x):
        self._pos[0] = x

    def set_y(self, y):
        self._pos[1] = y

    def set_z(self, z):
        self._pos[2] = z

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
        if tn and self.is_regular() and (isinstance(tn, EuclideanPoint) or tn.is_regular()):
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
        return '%d %d %g %g %g %g %d\n' % (self._id, self._type, self._pos[0], self._pos[1], self._pos[2], self._radius, self.parent.get_id())

    def get_parent_id(self):
        return -2 if self.is_root else self.parent.get_id()

    def add_child(self, swc_node):
        if not isinstance(swc_node, SwcNode):
            return False
        children = list(self.children)
        children.append(swc_node)
        self.children = tuple(children)

    def remove_child(self, swc_node):
        if not isinstance(swc_node, SwcNode):
            return False
        children = list(self.children)
        children.remove(swc_node)
        self.children = tuple(children)

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
            #dededebug
            # if tn.parent is None:
            #     print("tn = {} tn.pa = None dis = {}".format(
            #         tn.get_id(), tn.parent_distance()
            #     ))
            # else:
            #     print("tn = {} tn.pa = {} dis = {}".format(
            #         tn.get_id(),tn.parent.get_id(),tn.parent_distance()
            #     ))
            if tn.parent is not None and tn.get_id() > tn.parent.get_id():
                result += tn.parent_distance()

        return result

    def radius(self, nid):
        return self.node(nid).radius()

    def get_depth_array(self):
        self.depth_array = [0] * (self.node_count()+10)
        for node in PreOrderIter(self.root()):
            self.depth_array[node.get_id()] = node.depth()

    # initialize LCA data structure in swc_tree
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


if __name__ == '__main__':
    # print('testing ...')
    # tree = SwcTree()
    # tree.load("D:\gitProject\mine\PyMets\\test\data_example\gold\gold.swc")
    # tree.get_lca_preprocess()
    # print(tree.get_lca(2,6))

    # tree = SwcTree()
    # tree.load("D:\gitProject\mine\PyMets\\test\data_example\gold\gold.swc")
    # node = tree.root()
    # tmp = list(node.children)
    # # tmp.pop()
    #
    # node.children = tuple(tmp)
    # tmp = [node for node in PreOrderIter(tree.root())]
    # for node in tmp:
    #     print(node.get_id())
    pass