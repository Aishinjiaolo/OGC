import os, sys
src_folder = os.path.abspath('./') + '/src'
if src_folder not in sys.path:
    sys.path.insert(0, src_folder)

import gdspy

g_cell = gdspy.Cell('polygon')
points = [(0, 0), (2, 2), (2, 6), (-6, 6), (-6, -6), (-4, -4), (-4, 4), (0, 4)]
polygon1 = gdspy.Polygon(1, points)
g_cell.add(polygon1)

name = os.path.abspath(os.path.dirname(os.sys.argv[0])) + os.sep + 'gdspy_sample'
gdspy.gds_print(name + '.gds', unit = 1.0e-6, precision = 1.0e-9)
print 'Sample gds file saved: ' + name + '.gds'


