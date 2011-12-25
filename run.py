#!/usr/bin/python

import os, sys
src_folder  = os.path.abspath('./') + '/py_lib'
swig_folder = os.path.abspath('./') + '/swig'
if src_folder not in sys.path:
    sys.path.insert(0, src_folder)
if swig_folder not in sys.path:
    sys.path.insert(0, swig_folder)

import gdspy
import point
import segment
import polygon
import api

profiler = api.Profiler()
profiler.start("python flow")

point1 = point.Point()
point1.set_point(-10, -20)
print "point1 = ", point1.get_gx(), point1.get_gy()

point2 = point.Point()
point2.set_point(30, -40)
print "point2 = ", point2.get_gx(), point2.get_gy()

point3 = point.Point()
point3.set_point(50, 60)
print "point3 = ", point3.get_gx(), point3.get_gy()

point4 = point.Point()
point4.set_point(-70, 80)
print "point4 = ", point4.get_gx(), point4.get_gy()

segment1 = segment.Segment()
segment1.set_segment(point1, point2)
print "dump segment1:\n", segment1.dump()

segment2 = segment.Segment()
segment2.set_segment(point2, point3)
print "dump segment2:\n", segment2.dump()

segment3 = segment.Segment()
segment3.set_segment(point3, point4)
print "dump segment3:\n", segment3.dump()

segment4 = segment.Segment()
segment4.set_segment(point4, point1)
print "dump segment4:\n", segment4.dump()

segments1 = segment.Segments()
segments1.append(segment1)
segments1.append(segment2)
segments1.append(segment3)
segments1.append(segment4)
print "dump segments1:\n", segments1.dump()

polygon_1 = polygon.Polygon()
polygon_1.set_polygon(segments1)
print "dump polygon_1\n", polygon_1.dump()

print "before init another polygon"
polygon_copy = polygon.Polygon()
print "after init another polygon"
polygon_copy.copy(polygon_1)
print "dump polygon_copy\n", polygon_copy.dump()

polygon_1.get_segment(0).get_head().set_point(100, 100)
print "dump polygon_1 after modified\n", polygon_1.dump()
print "dump polygon_copy after polygon_1 modified", polygon_copy.dump()

# this should be wrapped from c
def add_polygon(polygon):
    points_array = []
    length = polygon.get_segment_number()
    for index in range(length):
       point = (polygon.get_segment(index).get_head().get_gx(), \
                polygon.get_segment(index).get_head().get_gy())
       points_array.append(point)
    return points_array

points_array = add_polygon(polygon_copy)
print "points array: ", points_array

g_cell = gdspy.Cell('polygon')
points = [(0, 0), (2, 2), (2, 6), (-6, 6), (-6, -6), (-4, -4), (-4, 4), (0, 4)]
polygon1 = gdspy.Polygon(1, points_array)
polygon2 = gdspy.Polygon(2, points)
g_cell.add(polygon1)
g_cell.add(polygon2)

name = os.path.abspath(os.path.dirname(os.sys.argv[0])) + os.sep + 'sample'
gdspy.gds_print(name + '.gds', unit = 1.0e-6, precision = 1.0e-9)
print 'Sample gds file saved: ' + name + '.gds'

polygon_copy.free()
profiler.end()
