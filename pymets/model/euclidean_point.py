import numpy as np
import math


class EuclideanPoint(object):
    def __init__(self,
                 center=[0, 0, 0]):
        self._pos=center

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

    def to_str(self):
        print("EuclideanPoint: {}".format(self._pos))

    def add_coord(self, point_a):
        self._pos[0] += point_a.get_x()
        self._pos[1] += point_a.get_y()
        self._pos[2] += point_a.get_z()

    def get_foot_point(self, line):
        if len(line.coords) != 2:
            raise Exception("[Error: ]in function get_foot_point. read line error")

        p = np.array(self._pos)
        a = np.array(line.coords[0])
        b = np.array(line.coords[1])
        if (np.abs(a-b) < 1e-8).all():
            print(a-b)
            return EuclideanPoint(a.tolist())
        k = -((a - p).dot(b - a)) / (((b - a) ** 2).sum())
        foot = k * (b - a) + a
        return EuclideanPoint(foot.tolist())

    def get_foot_point_fast(self, line):
        if len(line.coords) != 2:
            raise Exception("[Error: ]in function get_foot_point. read line error")

        a_p = [line.coords[0][0] - self._pos[0], line.coords[0][1] - self._pos[1], line.coords[0][2] - self._pos[2]]
        b_a = [line.coords[1][0] - line.coords[0][0], line.coords[1][1] - line.coords[0][1], line.coords[1][2] - line.coords[0][2]]
        k_up = a_p[0]*b_a[0] + a_p[1]*b_a[1] + a_p[2]*b_a[2]
        k_down = b_a[0]**2+b_a[1]**2+b_a[2]**2
        k = -k_up/k_down
        foot = [k*b_a[0]+line.coords[0][0], k*b_a[1]+line.coords[0][1], k*b_a[2]+line.coords[0][2]]
        return EuclideanPoint(foot)

    def on_line(self, line):
        p = np.array(self._pos)
        a = np.array(line.coords[0])
        b = np.array(line.coords[1])

        if min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and \
                min(a[1], b[1]) <= p[1] <= max(a[1], b[1]) and \
                min(a[2], b[2]) <= p[2] <= max(a[2], b[2]):
            return True
        return False

    def distance_to_point(self, point):
        if isinstance(point, list) and len(point) == 3:
            point = EuclideanPoint(point)
        if type(point) != type(EuclideanPoint()):
            raise Exception("[Error:  ] expect point. got {}".format(point))

        sub = [self._pos[0] - point._pos[0],
               self._pos[1] - point._pos[1],
               self._pos[2] - point._pos[2]]
        return math.sqrt(sub[0]*sub[0] + sub[1]*sub[1] + sub[2]*sub[2])

    def distance_to_line(self, line):
        foot = self.get_foot_point(line)
        return self.distance_to_point(foot)

    def distance_to_segment(self, line):
        foot = self.get_foot_point(line)
        # foot = EuclideanPoint([0, 0, 0])
        if foot.on_line(line):
            return self.distance_to_point(foot)
        else:
            dis1 = self.distance_to_point(line.coords[0])
            dis2 = self.distance_to_point(line.coords[1])
            return min(dis1,dis2)

    def distance(self, object):
        if isinstance(object, EuclideanPoint):
            return self.distance_to_point(object)
        elif isinstance(object, Line):
            if object.is_segment:
                return self.distance_to_segment(object)
            else:
                return self.distance_to_line(object)
        else:
            raise Exception("[Error: ] unexpected object type {}".format(type(object)))


class Line:
    def __init__(self,
                 coords=None,
                 swc_node_1=None,
                 swc_node_2=None,
                 is_segment=True):
        if coords is not None:
            self.coords = coords
        else:
            self.coords = [[],[]]
            self.coords[0] = swc_node_1._pos
            self.coords[1] = swc_node_2._pos
        self.is_segment = is_segment

    def to_str(self):
        print("Line: {}, {}".format(self.coords[0],self.coords[1]))

    def get_points(self):
        point_a = EuclideanPoint(self.coords[0])
        point_b = EuclideanPoint(self.coords[1])
        return point_a, point_b

    def distance(self, object):
        if type(object) == type(EuclideanPoint()):
            return object.disdance(self)
        elif type(object) == type(Line()):
            raise Exception("[Error: ] distance between lines is not support temporarily")
        else:
            raise Exception("[Error: ] unexpected object type" + type(object))


if __name__ == "__main__":
    points = [
        EuclideanPoint([0,0,0]),
        EuclideanPoint([0,0,1]),
        EuclideanPoint([0,1,0]),
        EuclideanPoint([1,0,0]),
        EuclideanPoint([1,2,3]),
        EuclideanPoint([2, 2, 2])
    ]
    for i in range(0,len(points)):
        for j in range(i+1, len(points)):
            print(points[i].distance(points[j]))
    print("----")
    line_a = Line(coords=[[0,0,0],[1,1,1]],is_segment=False)
    for point in points:
        print(point.distance(line_a))



