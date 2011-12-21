#!/usr/bin/python

import os, sys
src_folder  = os.path.abspath('./') + '/py_lib'
swig_folder = os.path.abspath('./') + '/swig'
if src_folder not in sys.path:
    sys.path.insert(0, src_folder)
if swig_folder not in sys.path:
    sys.path.insert(0, swig_folder)

import gdspy
import api
import polygon
import segment
import point

point1 = api.Point()
point1.set_point(10, 20)
print "point1 = ", point1.get_gx(), point1.get_gy()

g_cell = gdspy.Cell('polygon')
points = [(0, 0), (2, 2), (2, 6), (-6, 6), (-6, -6), (-4, -4), (-4, 4), (0, 4)]
polygon1 = gdspy.Polygon(1, points)
g_cell.add(polygon1)

name = os.path.abspath(os.path.dirname(os.sys.argv[0])) + os.sep + 'sample'
gdspy.gds_print(name + '.gds', unit = 1.0e-6, precision = 1.0e-9)
print 'Sample gds file saved: ' + name + '.gds'


