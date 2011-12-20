########################################################################
##                                                                    ##
##  Copyright 2009-2011 Lucas Heitzmann Gabrielli                     ##
##                                                                    ##
##  This file is part of gdspy.                                       ##
##                                                                    ##
##  gdspy is free software: you can redistribute it and/or modify it  ##
##  under the terms of the GNU General Public License as published    ##
##  by the Free Software Foundation, either version 3 of the          ##
##  License, or any later version.                                    ##
##                                                                    ##
##  gdspy is distributed in the hope that it will be useful, but      ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of        ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
##  GNU General Public License for more details.                      ##
##                                                                    ##
##  You should have received a copy of the GNU General Public         ##
##  License along with gdspy.  If not, see                            ##
##  <http://www.gnu.org/licenses/>.                                   ##
##                                                                    ##
########################################################################

import numpy

__doc__ = """
Boolean is a Python module that implements boolean operations between
polygons (polygon clipping).

It is distributed as part of the *gdspy* module to help the creation of
complex structures in the GDSII stream format, but the function
interfaces should allow its usage in more general situations.
"""


class _Point:
    """
    Holds information about the position of a point in the graph
    structure, as well as to which segments it belongs.

    Parameters
    ----------
    xy : array-like[2]
        Position of the point.
    segments : list
        List of segments which this point belongs to.
    """
    def __init__(self, xy):
        self.xy = numpy.array(xy, dtype=long)
        self.segments = []
    
    def __repr__(self):
        """
        Return a string representation of this instance.

        Returns
        -------
        str : string
            String representation of this instance.
        """
        return '(' + str(self.xy[0]) + ', ' + str(self.xy[1]) + ')'

def compare_points(point1, point2):
    """
    Compare points.
    
    Points are ordered by first coordinate and, if equal, by second
    coordinate.
    
    Parameters
    ----------
    point1, point2 : _Point
        Points to be compared.
    
    Returns
    -------
    out : int
        Returns -1 if ``point1 < point2``, +1 if ``point1 > point2``,
        and 0 if they are equal.
    """
    if point1.xy[0] < point2.xy[0]:
        return -1
    elif point1.xy[0] > point2.xy[0]:
        return 1
    elif point1.xy[1] < point2.xy[1]:
        return -1
    elif point1.xy[1] > point2.xy[1]:
        return 1
    else:
        return 0


def compare_pair_0(point1, point2):
    """
    Compare points.
    
    Points are ordered by first coordinate.
    
    Parameters
    ----------
    point1, point2 : array-like[2]
        Points to be compared.
    
    Returns
    -------
    out : int
        Returns -1 if ``point1[0] < point2[0]``, +1 if ``point1[0] > 
        point2[0]``, and 0 if they are equal.
    """
    if point1[0] < point2[0]:
        return -1
    elif point1[0] > point2[0]:
        return 1
    else:
        return 0


def compare_pair_1(point1, point2):
    """
    Compare points.
    
    Points are ordered by second coordinate.
    
    Parameters
    ----------
    point1, point2 : array-like[2]
        Points to be compared.
    
    Returns
    -------
    out : int
        Returns -1 if ``point1[1] < point2[1]``, +1 if ``point1[1] >
        point2[1]``, and 0 if they are equal.
    """
    if point1[1] < point2[1]:
        return -1
    elif point1[1] > point2[1]:
        return 1
    else:
        return 0


class _Segment:
    """
    Holds information about a segment in the graph structure.

    Parameters
    ----------
    polygon : array of integers
        The value ``polygon[n]`` identifies whether this segment belongs to
        polygon ``n`` in the graph.
        If 0, the segment does not belong to polygon ``n``. If 1, the segment
        belongs and is oriented from ``points[0]`` to ``points[1]``. If -1, it
        belongs with the reverse orientation.
    points : list[2] of _Point
        Begin and end points of this segment.
    in_phase : bool
        Must be set to ``True`` if ``points[0] < points[1]``, ``False``
        otherwise.
    """
    def __init__(self, polygon, points, in_phase=False):
        self.polygon = numpy.array(polygon, dtype=int)
        self.points = points
        self.in_phase = in_phase

    def __repr__(self):
        """
        Return a string representation of this instance.

        Returns
        -------
        str : string
            String representation of this instance.
        """
        return str(self.points[0]) + '->' + str(self.points[1]) + ' @ ' + str(self.polygon)
    
def higher_segment(segment1, segment2):
    """
    Compare the y component of the segments at a given x coordinate.

    The x coordinate is chosen as the center of the interval where
    the segments overlap in the x coordinate.
    If such an interval does not exist, 0 is returned.
    
    Parameters
    ----------
    segment1, segment2 : _Segment
        Segment to be compared.
    
    Returns
    -------
    out : int
        Returns -1 if ``segment1`` is lower than ``segment2``, +1 if
        ``segment1`` is higher than ``segment2``, and 0 if they are
        at same heights.
    """
    x_list = [segment1.points[0].xy[0], segment1.points[1].xy[0], segment2.points[0].xy[0], segment2.points[1].xy[0]]
    x_list.sort()
    x = (x_list[1] + x_list[2]) * 0.5
    dx1 = long(segment1.points[1].xy[0] - segment1.points[0].xy[0])
    dx2 = long(segment2.points[1].xy[0] - segment2.points[0].xy[0])
    if dx1 == 0:
        y1 = long(segment1.points[0].xy[1] + segment1.points[1].xy[1]) // 2L
    else:
        y1 = segment1.points[0].xy[1] + (segment1.points[1].xy[1] - segment1.points[0].xy[1]) * (x - segment1.points[0].xy[0]) / (segment1.points[1].xy[0] - segment1.points[0].xy[0])
    if dx2 == 0:
        y2 = long(segment2.points[0].xy[1] + segment2.points[1].xy[1]) // 2L
    else:
        y2 = segment2.points[0].xy[1] + (segment2.points[1].xy[1] - segment2.points[0].xy[1]) * (x - segment2.points[0].xy[0]) / (segment2.points[1].xy[0] - segment2.points[0].xy[0])
    if y1 > y2:
        return 1
    elif y1 < y2:
        return -1
    return 0


def _merge_overlaps(points, segments=None):
    """
    Merge overlapping points and segments in the graph structure.

    Parameters
    ----------
    points : list of _Point
        List of points in the graph structure.
    segments : list of _Segment or None
        If not ``None`` overlapping segments will also be merged.
    """
    ii = 1
    while ii < len(points):
        if compare_points(points[ii-1], points[ii]) == 0:
            points[ii-1].segments += points[ii].segments
            for seg in points[ii].segments:
                if seg.points[0] is points[ii]:
                    seg.points[0] = points[ii-1]
                else:
                    seg.points[1] = points[ii-1]
            points.remove(points[ii])
        else:
            if segments is not None:
                if points[ii].segments[0].points[0] is points[ii]:
                    out_nodes = [points[ii].segments[0].points[1]]
                else:
                    out_nodes = [points[ii].segments[0].points[0]]
                jj = 1
                while jj < len(points[ii].segments):
                    if points[ii].segments[jj].points[0] is points[ii]:
                        cur_node = points[ii].segments[jj].points[1]
                    else:
                        cur_node = points[ii].segments[jj].points[0]
                    kk = 0
                    while kk < len(out_nodes):
                        if cur_node is out_nodes[kk]:
                            if points[ii].segments[kk].points[0] is points[ii].segments[jj].points[0]:
                                points[ii].segments[kk].polygon += points[ii].segments[jj].polygon
                            else:
                                points[ii].segments[kk].polygon -= points[ii].segments[jj].polygon
                            segments.remove(points[ii].segments[jj])
                            cur_node.segments.remove(points[ii].segments[jj])
                            points[ii].segments.remove(points[ii].segments[jj])
                            kk = len(out_nodes) + 1
                        else:
                            kk += 1                                                                    
                    if kk == len(out_nodes):
                        out_nodes.append(cur_node)
                        jj += 1
            ii += 1


def _intersect_point_segment(points, segments):
    """
    Calculates all intersections of segments with points.

    Parameters
    ----------
    points : list of _Point
        List of points in the graph structure.
    segments : list of _Segment
        List of segments in the graph structure.

    Returns
    -------
    out : bool
        ``True`` if any intersections are found. ``False`` otherwise.
    """
    changed = False
    seg_check = list(points[0].segments)
    for ii in range(1, len(points)):
        ss = 0
        while ss < len(seg_check):
            dx0 = long(seg_check[ss].points[1].xy[0] - seg_check[ss].points[0].xy[0])
            dy0 = long(seg_check[ss].points[1].xy[1] - seg_check[ss].points[0].xy[1]) 
            dx1 = long(points[ii].xy[0] - seg_check[ss].points[0].xy[0])
            dy1 = long(points[ii].xy[1] - seg_check[ss].points[0].xy[1])
            seg_len2 = dx0 * dx0 + dy0 * dy0
            cross2 = (dx1 * dy0 - dx0 * dy1)**2
            dot = dx1 * dx0 + dy1 * dy0

            if (points[ii] not in seg_check[ss].points) and (cross2 < seg_len2) and (dot > 0) and (dot < seg_len2):
                changed = True
                segments.append(_Segment(seg_check[ss].polygon, [seg_check[ss].points[0], points[ii]], compare_points(points[ii], seg_check[ss].points[0]) > 0))
                seg_check[ss].points[0].segments.remove(seg_check[ss])
                seg_check[ss].points[0].segments.append(segments[-1])
                seg_check[ss].points[0] = points[ii]
                seg_check[ss].in_phase = (compare_points(seg_check[ss].points[1], points[ii]) > 0)
                points[ii].segments += [seg_check[ss], segments[-1]]
                if (segments[-1].in_phase and (compare_points(segments[-1].points[1], points[ii]) >= 0) and (compare_points(segments[-1].points[0], points[ii-1]) <= 0)) or ((not segments[-1].in_phase) and (compare_points(segments[-1].points[0], points[ii]) >= 0) and (compare_points(segments[-1].points[1], points[ii-1]) <= 0)):
                    seg_check.append(segments[-1])
                if not ((seg_check[ss].in_phase and (compare_points(seg_check[ss].points[1], points[ii]) >= 0) and (compare_points(seg_check[ss].points[0], points[ii-1]) <= 0)) or ((not seg_check[ss].in_phase) and (compare_points(seg_check[ss].points[0], points[ii]) >= 0) and (compare_points(seg_check[ss].points[1], points[ii-1]) <= 0))):
                    seg_check.remove(seg_check[ss])
            else:
                ss += 1
        ## Update the scanbeam ``seg_check``
        for seg in points[ii].segments:
            if ((seg.points[0] is points[ii]) and seg.in_phase) or ((seg.points[1] is points[ii]) and not seg.in_phase):
                seg_check.append(seg)
            elif ((seg.points[0] is points[ii]) and not seg.in_phase) or ((seg.points[1] is points[ii]) and seg.in_phase):
                seg_check.remove(seg)
    return changed
    

def _polygon_area(polygon):
    """
    Calculate the (counter-clockwise) area of a polygon.

    Parameters
    ----------
    polygon : array-like[N][2]
        List of the polygon vertices.

    Returns
    -------
    out : number
        Area of the polygon.
    """
    poly_area = 0
    for ii in range(1, len(polygon) - 1):
        poly_area += long(polygon[ii][0] - polygon[0][0]) * long(polygon[ii+1][1] - polygon[0][1]) - long(polygon[ii][1] - polygon[0][1]) * long(polygon[ii+1][0] - polygon[0][0])
    return poly_area // 2


def _inside(point, polygon):
    """
    Determine if a point is located inside a polygon.

    Parameters
    ----------
    point : array-like[2]
        Coordinates of the point to be checked.
    polygon : array-like[N][2]
        List of the polygon vertices.

    Returns
    -------
    out : bool
        ``True`` if the point is inside the polygon, ``False`` otherwise.
    """
    odd = False
    for ii in range(len(polygon)):
        if ((polygon[ii][1] < point[1]) and (polygon[ii-1][1] >= point[1])) or ((polygon[ii-1][1] < point[1]) and (polygon[ii][1] >= point[1])):
            if (polygon[ii-1][0] - polygon[ii][0]) * (point[1] - polygon[ii][1]) / (polygon[ii-1][1] - polygon[ii][1]) - point[0] + polygon[ii][0] < 0:
                odd = (not odd)
    return odd


def _add_polygon(new_polygon, out_polygons, bounding_boxes):
    """
    Append a new polygon to a list of polygons. If the polygon is in
    clockwise direction, it is reversed.

    Parameters
    ----------
    new_polygon : array-like[N][2]
        List of vertices of the polygon to be added to the list.
    out_polygons : list of array-like[N][2]
        List of polygons where ``new_polygon`` will be inserted.
    bounding_boxes : list of array-like[4]
        List of bounding boxes of the polygons. Each bounding box is
        an array with [min x, max x, min y, max y] coordinates.
    """
    if _polygon_area(new_polygon) < 0:
        new_polygon.reverse()
    min_x = new_polygon[0][0]
    max_x = new_polygon[0][0]
    min_y = new_polygon[0][1]
    max_y = new_polygon[0][1]
    for ii in range(1, len(new_polygon)):
        if new_polygon[ii][1] < min_y:
            min_y = new_polygon[ii][1]
        elif new_polygon[ii][1] > max_y:
            max_y = new_polygon[ii][1]
        if new_polygon[ii][0] < min_x:
            min_x = new_polygon[ii][0]
        elif new_polygon[ii][0] > max_x:
            max_x = new_polygon[ii][0]
    out_polygons.append(new_polygon)
    bounding_boxes.append([min_x, max_x, min_y, max_y])


def chop(polygon, position, axis):
    """
    Slice polygon at a given position along a given axis.
    
    Parameters
    ----------
    polygon : array-like[N][2]
        Coordinates of the vertices of the polygon.
    position : number
        Position to perform the slicing operation along the specified
        axis.
    axis : 0 or 1
        Axis along which the polygon will be sliced.
    
    Returns
    -------
    out : tuple[2]
        Each element is a list of polygons (array-like[N][2]).  The first
        list contains the polygons left before the slicing position, and
        the second, the polygons left after that position.
    """
    out_polygons = ([], [])
    polygon = list(polygon)
    while polygon[-1][axis] == position:
        polygon = [polygon[-1]] + polygon[:-1]
    cross = list(numpy.sign(numpy.array(polygon)[:, axis] - position))
    bnd = ([], [])
    i = 0
    while i < len(cross):
        if cross[i - 1] * cross[i] < 0:
            if axis == 0:
                polygon.insert(i, [position, polygon[i - 1][1] + (position - polygon[i - 1][0]) * (polygon[i][1] - polygon[i - 1][1]) / float(polygon[i][0] - polygon[i - 1][0])])
            else:
                polygon.insert(i, [polygon[i - 1][0] + (position - polygon[i - 1][1]) * (polygon[i][0] - polygon[i - 1][0]) / float(polygon[i][1] - polygon[i - 1][1]), position])
            cross.insert(i, 0)
            bnd[1 * (cross[i + 1] > 0)].append(i)
            i += 2
        elif cross[i] == 0:
            j = i + 1
            while cross[j] == 0:
                j += 1
            if cross[i - 1] * cross[j] < 0:
                bnd[1 * (cross[j] > 0)].append(j - 1)
            i = j + 1
        else:
            i += 1
    if len(bnd[0]) == 0:
        out_polygons[1 * (numpy.sum(cross) > 0)].append(polygon)
        return out_polygons
    bnd = (numpy.array(bnd[0]), numpy.array(bnd[1]))
    bnd = (list(bnd[0][numpy.argsort(numpy.array(polygon)[bnd[0], 1 - axis])]),
           list(bnd[1][numpy.argsort(numpy.array(polygon)[bnd[1], 1 - axis])]))
    cross = numpy.ones(len(polygon), dtype=int)
    cross[bnd[0]] = -2
    cross[bnd[1]] = -1
    i = 0
    while i < len(polygon):
        if cross[i] > 0 and polygon[i][axis] != position:
            start = i
            side = 1 * (polygon[i][axis] > position)
            out_polygons[side].append([polygon[i]])
            cross[i] = 0
            nxt = i + 1
            if nxt == len(polygon):
                nxt = 0
            boundary = True
            while nxt != start:
                out_polygons[side][-1].append(polygon[nxt])
                if cross[nxt] > 0:
                    cross[nxt] = 0
                if cross[nxt] < 0 and boundary:
                    j = bnd[cross[nxt] + 2].index(nxt)
                    nxt = bnd[-cross[nxt] - 1][j]
                    boundary = False
                else:
                    nxt += 1
                    if nxt == len(polygon):
                        nxt = 0
                    boundary = True
        i += 1
    return out_polygons


def clip(polygons, operation, precision, max_points=0):
    """
    Perform a boolean operation (clipping) on a set of polygons.
    
    Parameters
    ----------
    polygons : list of array-like[N][2]
        List of polygons. Each polygon is an array-like[N][2] object
        with the coordinates of the vertices of the polygon.
    operation : function
        This function should accept a list[N] of integer variables as
        argument and output an integer or boolean representing the
        desired operation to be performed on the polygons. Each input
        variable corresponds to one polygon of ``polygons``.
    precision : number
        The snapping distance for the vertices of the polygons (grid
        spacing).
    max_points : integer
        If greater than 4, fractures the resulting polygons to ensure
        they have at most ``max_points`` vertices. This is not a
        tessellating function, so this number should be as high as
        possible. For example, it should be set to 199 for polygons
        being drawn in GDSII files.
    
    Returns
    -------
    out : list of array-like[N][2]
        List of polygons resulting from the boolean operation.
    """
    ## Create a graph structure with all polygons
    points = []
    segments = []
    
    for pp in range(len(polygons)):
        pol_mask = numpy.zeros(len(polygons))
        pol_mask[pp] = 1
        if polygons[pp].__class__ == numpy.ndarray:
            polygons[pp] = polygons[pp].tolist()
        
        ## Round points to ``precision``
        for ii in range(len(polygons[pp])):
            polygons[pp][ii] = [long(round(polygons[pp][ii][0] / precision)),  long(round(polygons[pp][ii][1] / precision))]
        
        ## Remove 0-area triangles formed by three consective points
        removed = True
        while removed:
            removed = False
            ii = 0
            while ii < len(polygons[pp]):
                if (polygons[pp][ii-2][0] - polygons[pp][ii-1][0]) * (polygons[pp][ii][1] - polygons[pp][ii-1][1]) - (polygons[pp][ii][0] - polygons[pp][ii-1][0]) * (polygons[pp][ii-2][1] - polygons[pp][ii-1][1]) == 0:
                    polygons[pp].pop(ii-1)
                    removed = True
                else:
                    ii += 1
            if len(polygons[pp]) < 3:
                removed = False
        
        ## Add polygon to the graph structure
        if len(polygons[pp]) > 2:
            for ii in range(len(polygons[pp])):
                points.append(_Point(polygons[pp][ii]))
                if ii != 0:
                    segments.append(_Segment(pol_mask, [points[-2], points[-1]]))
                    points[-1].segments.append(segments[-1])
                    points[-2].segments.append(segments[-1])
            segments.append(_Segment(pol_mask, [points[-1], points[-len(polygons[pp])]]))
            points[-1].segments.append(segments[-1])
            points[-len(polygons[pp])].segments.append(segments[-1])
    
    ## Sort points by x-coordinate (y if same x), and orient segments
    points.sort(compare_points)
    for seg in segments:
        if compare_points(seg.points[0], seg.points[1]) < 0:
            seg.in_phase = True

    ## Merge overlaping points
    _merge_overlaps(points)

    ## Calculate intersections between segments and points
    changed = True
    while changed:
        changed = False
        for ii in range(2):
            changed = _intersect_point_segment(points, segments)
            for point in points:
                point.xy[1], point.xy[0] = point.xy
            points.sort(compare_points)
            for seg in segments:
                if compare_points(seg.points[0], seg.points[1]) < 0:
                    seg.in_phase = True
                else:
                    seg.in_phase = False

    ## Add points to segment crossings
    new_points = []
    seg_check = list(points[0].segments)
    for point in points:
        ## Update seg_check
        for seg in point.segments:
            if ((seg.points[0] is point) and seg.in_phase) or ((seg.points[1] is point) and not seg.in_phase):
                seg_check.append(seg)
            elif ((seg.points[0] is point) and not seg.in_phase) or ((seg.points[1] is point) and seg.in_phase):
                seg_check.remove(seg)
                ## Check for intersections
                dx1 = long(seg.points[1].xy[0] - seg.points[0].xy[0])
                dy1 = long(seg.points[1].xy[1] - seg.points[0].xy[1])
                for seg0 in seg_check:
                    dx0 = long(seg0.points[1].xy[0] - seg0.points[0].xy[0])
                    dy0 = long(seg0.points[1].xy[1] - seg0.points[0].xy[1])
                    den = long(dx0 * dy1 - dx1 * dy0)
                    if den != 0:
                        dx = long(seg.points[1].xy[0] - seg0.points[0].xy[0])
                        dy = long(seg.points[1].xy[1] - seg0.points[0].xy[1])
                        alpha = float(dy1 * dx - dx1 * dy) / den
                        beta = float(dx0 * dy - dy0 * dx) / den
                        if alpha > 0 and beta > 0 and alpha < 1 and beta < 1:
                            ## Insert crossing
                            x = dx0 * dx1 * long(seg0.points[0].xy[1] - seg.points[0].xy[1]) + dy1 * dx0 * long(seg.points[0].xy[0]) - dy0 * dx1 * long(seg0.points[0].xy[0])
                            y = -(dy0 * dy1 * long(seg0.points[0].xy[0] - seg.points[0].xy[0]) + dx1 * dy0 * long(seg.points[0].xy[1]) - dx0 * dy1 * long(seg0.points[0].xy[1]))
                            if abs(x - (x // den) * den) > abs(x - (x // den + 1) * den):
                                x += den
                            elif abs(x - (x // den) * den) > abs(x - (x // den - 1) * den):
                                x -= den
                            if abs(y - (y // den) * den) > abs(y - (y // den + 1) * den):
                                y += den
                            elif abs(y - (y // den) * den) > abs(y - (y // den - 1) * den):
                                y -= den
                            new_points.append(_Point((x // den, y // den)))
    
    ## Add crossings to graph
    for new_point in new_points:
        if compare_points(new_point, points[-1]) > 0:
            points.append(new_point)
        else:
            jj = 0
            while compare_points(new_point, points[jj]) > 0:
                jj += 1
            if compare_points(new_point, points[jj]) < 0:
                points.insert(jj, new_point)
    new_points = None

    ## Recalculate intersections between segments and points
    changed = True
    while changed:
        changed = False
        for ii in range(2):
            changed = _intersect_point_segment(points, segments)
            for point in points:
                point.xy[1], point.xy[0] = point.xy
            points.sort(compare_points)
            for seg in segments:
                if compare_points(seg.points[0], seg.points[1]) < 0:
                    seg.in_phase = True
                else:
                    seg.in_phase = False

    ## Merge overlaping points and segments
    _merge_overlaps(points, segments)
    
    ## Perform operation
    ii = 1
    seg_check = list(points[0].segments)
    seg_check.sort(higher_segment)
    seg_output = [0] * len(seg_check)
    while ii < len(points):
        phase = numpy.zeros(len(polygons), dtype=bool)
        counter = numpy.zeros(len(polygons), dtype=long)
        last_result = False
        for ss in range(len(seg_check) - 1):
            for jj in range(len(polygons)):
                if seg_check[ss].polygon[jj] != 0:
                    if counter[jj] == 0:
                        if seg_check[ss].polygon[jj] > 0:
                            phase[jj] = seg_check[ss].in_phase
                        else:
                            phase[jj] = (not seg_check[ss].in_phase)
                        counter[jj] += 1
                    else:
                        if seg_check[ss].polygon[jj] > 0:
                            if seg_check[ss].in_phase == phase[jj]:
                                counter[jj] += 1
                            else:
                                counter[jj] -= 1
                        else:
                            if seg_check[ss].in_phase == phase[jj]:
                                counter[jj] -= 1
                            else:
                                counter[jj] += 1
            result = (operation(counter) > 0)
            if result != last_result:
                if (seg_check[ss].in_phase and result) or ((not seg_check[ss].in_phase) and (not result)):
                    seg_output[ss] = 1
                else:
                    seg_output[ss] = -1
                last_result = result
        if (len(seg_output) > 0) and last_result:
            if seg_check[-1].in_phase:
                seg_output[-1] = -1
            else:
                seg_output[-1] = 1
        ss = 0
        while ss < len(points[ii].segments):
            seg = points[ii].segments[ss]
            if ((seg.points[0] is points[ii]) and seg.in_phase) or ((seg.points[1] is points[ii]) and (not seg.in_phase)):
                if len(seg_check) == 0:
                    seg_check.append(seg)
                    seg_output.append(0)
                elif higher_segment(seg, seg_check[-1]) > 0:
                    seg_check.append(seg)
                    seg_output.append(0)
                else:
                    jj = 0
                    while higher_segment(seg, seg_check[jj]) > 0:
                        jj += 1
                    seg_check.insert(jj, seg)
                    seg_output.insert(jj, 0)
                ss += 1
            else:
                jj = seg_check.index(seg)
                seg_check.remove(seg)
                if seg_output[jj] != 0:
                    seg.in_phase = (seg_output[jj] > 0)
                    ss += 1
                else:
                    segments.remove(points[ii].segments[ss])
                    seg.points[0].segments.remove(seg)
                    seg.points[1].segments.remove(seg)
                if jj == 0:
                    seg_output = seg_output[1:]
                elif jj == len(seg_output) - 1:
                    seg_output = seg_output[:-1]
                else:
                    seg_output = seg_output[:jj] + seg_output[jj+1:]
        ii += 1

    ## Find resulting polygons
    out_polygons = []
    bounding_boxes = []
    while len(segments) > 0:
        seg = segments.pop()
        seg.points[0].segments.remove(seg)
        seg.points[1].segments.remove(seg)
        new_polygon = [seg.points[0]]                     
        new_point = seg.points[1]
        true_left = seg.in_phase
        cross = []
        while new_point is not new_polygon[0]:
            if new_point in cross:
                ii = new_polygon.index(new_point)
                for pt in new_polygon[ii:]:
                    if (pt in cross) and (len(pt.segments) <= 2):
                        cross.remove(pt)
                _add_polygon([pt.xy.copy() for pt in new_polygon[ii:]], out_polygons, bounding_boxes)
                new_polygon = new_polygon[:ii]

            new_polygon.append(new_point)
            for seg in new_point.segments:
                if ((seg.points[0] is new_point) and (seg.in_phase == true_left)) or ((seg.points[1] is new_point) and (seg.in_phase != true_left)):
                    break

            segments.remove(seg)
            seg.points[0].segments.remove(seg)
            seg.points[1].segments.remove(seg)
            if len(new_point.segments) > 0:
                cross.append(new_point)
            if seg.points[0] is new_point:
                new_point = seg.points[1]
            else:
                new_point = seg.points[0]

        _add_polygon([pt.xy.copy() for pt in new_polygon], out_polygons, bounding_boxes)
    
    ## Create a hole information tree
    hole_tree = [[] for ii in range(len(out_polygons))]
    for ii in range(len(out_polygons) - 1):
        for jj in range(ii + 1, len(out_polygons)):
            if (bounding_boxes[jj][0] == bounding_boxes[ii][0]) and (bounding_boxes[jj][1] == bounding_boxes[ii][1]) and (bounding_boxes[jj][2] == bounding_boxes[ii][2]) and (bounding_boxes[jj][3] == bounding_boxes[ii][3]):
                if _inside((out_polygons[ii][0] + out_polygons[ii][1]) * 0.5, out_polygons[jj]):
                    hole_tree[jj].append(ii)
                elif _inside((out_polygons[jj][0] + out_polygons[jj][1]) * 0.5, out_polygons[ii]):
                    hole_tree[ii].append(jj)
            elif (bounding_boxes[jj][0] <= bounding_boxes[ii][0]) and (bounding_boxes[jj][1] >= bounding_boxes[ii][1]) and (bounding_boxes[jj][2] <= bounding_boxes[ii][2]) and (bounding_boxes[jj][3] >= bounding_boxes[ii][3]):
                if _inside((out_polygons[ii][0] + out_polygons[ii][1]) * 0.5, out_polygons[jj]):
                    hole_tree[jj].append(ii)
            elif (bounding_boxes[jj][0] >= bounding_boxes[ii][0]) and (bounding_boxes[jj][1] <= bounding_boxes[ii][1]) and (bounding_boxes[jj][2] >= bounding_boxes[ii][2]) and (bounding_boxes[jj][3] <= bounding_boxes[ii][3]):
                if _inside((out_polygons[jj][0] + out_polygons[jj][1]) * 0.5, out_polygons[ii]):
                    hole_tree[ii].append(jj)
    is_hole = [False] * len(hole_tree)
    for ii in range(len(hole_tree)):
        for jj in range(len(hole_tree[ii])):
            is_hole[hole_tree[ii][jj]] = (not is_hole[hole_tree[ii][jj]])
    for ii in range(len(hole_tree)):
        if not is_hole[ii]:
            jj = len(hole_tree[ii]) - 1
            while jj >= 0:
                if is_hole[hole_tree[ii][jj]]:
                    kk = len(hole_tree[ii]) - 1
                    while kk >= 0:
                        if is_hole[hole_tree[ii][kk]] and (hole_tree[ii][jj] in hole_tree[hole_tree[ii][kk]]):
                            hole_tree[ii].remove(hole_tree[ii][jj])
                            kk = -1
                        kk -= 1
                else: 
                    hole_tree[ii].remove(hole_tree[ii][jj])
                jj -= 1

    ## Connect holes to their outer polygons
    for ii in range(len(out_polygons)):
        if (not is_hole[ii]):
            ## Connect any holes that share a node with the outer polygon
            jj = len(hole_tree[ii]) - 1
            while jj >= 0:
                kk = hole_tree[ii][jj]

                ## Check for connection with outer polygon
                ip = len(out_polygons[ii]) - 1
                while ip >= 0:
                    kp = len(out_polygons[kk]) - 1
                    while kp >= 0:
                        if (out_polygons[ii][ip][0] == out_polygons[kk][kp][0]) and (out_polygons[ii][ip][1] == out_polygons[kk][kp][1]):
                            out_polygons[kk].reverse()
                            kp = len(out_polygons[kk]) - 1 - kp
                            out_polygons[ii] = out_polygons[ii][:ip] + out_polygons[kk][kp:] + out_polygons[kk][:kp] + out_polygons[ii][ip:]
                            out_polygons[kk] = []
                            hole_tree[ii].remove(kk)
                            kp = -1
                            ip = -1
                        kp -= 1
                    ip -= 1
                if ip > -2:
                    ## Check for connection with other holes
                    ll = jj - 1
                    while ll >= 0:
                        mm = hole_tree[ii][ll]
                        mp = len(out_polygons[mm]) - 1
                        while mp >= 0:
                            kp = len(out_polygons[kk]) - 1
                            while kp >= 0:
                                if (out_polygons[mm][mp][0] == out_polygons[kk][kp][0]) and (out_polygons[mm][mp][1] == out_polygons[kk][kp][1]):
                                    out_polygons[mm] = out_polygons[mm][:mp] + out_polygons[kk][kp:] + out_polygons[kk][:kp] + out_polygons[mm][mp:]
                                    out_polygons[kk] = []
                                    hole_tree[ii].remove(kk)
                                    kp = -1
                                    mp = -1
                                    ll = -1
                                kp -= 1
                            mp -= 1
                        ll -= 1
                jj -= 1
            if len(hole_tree[ii]) > 0:
                ## Sort holes by their bounding boxes
                sort_funtion = (lambda i,j: int(2 * (bounding_boxes[i][0] > bounding_boxes[j][0]) - 1))
                hole_tree[ii].sort(sort_funtion)
                ## Connect holes to the boundary
                for jj in hole_tree[ii]:
                    out_polygons[jj].reverse()
                    jp = 0
                    for pp in range(1, len(out_polygons[jj])):
                        if (out_polygons[jj][jp][0] > out_polygons[jj][pp][0]) or ((out_polygons[jj][jp][0] == out_polygons[jj][pp][0]) and (out_polygons[jj][jp][1] > out_polygons[jj][pp][1])):
                            jp = pp
                    ## Intersect outer polygon
                    xp = -1e30
                    ip = -1
                    for pp in range(len(out_polygons[ii])):
                        x = xp
                        if out_polygons[ii][pp][1] == out_polygons[jj][jp][1]:
                            x = out_polygons[ii][pp][0]
                        elif ((out_polygons[ii][pp][1] > out_polygons[jj][jp][1]) and (out_polygons[ii][pp-1][1] < out_polygons[jj][jp][1])) or ((out_polygons[ii][pp][1] < out_polygons[jj][jp][1]) and (out_polygons[ii][pp-1][1] > out_polygons[jj][jp][1])):
                            x = out_polygons[ii][pp][0] + (out_polygons[ii][pp-1][0] - out_polygons[ii][pp][0]) * (out_polygons[jj][jp][1] - out_polygons[ii][pp][1]) / (out_polygons[ii][pp-1][1] - out_polygons[ii][pp][1])
                        if (x > xp) and (x < out_polygons[jj][jp][0]):
                            ip = pp
                            xp = x
                    if out_polygons[ii][ip][1] == out_polygons[jj][jp][1]:
                        ## Connect to existing vertex
                        out_polygons[ii] = out_polygons[ii][:ip+1] + out_polygons[jj][jp:] + out_polygons[jj][:jp+1] + out_polygons[ii][ip:]
                    else:
                        ## Connect to new vertex
                        out_polygons[ii] = out_polygons[ii][:ip] + [numpy.array([long(round(xp)), out_polygons[jj][jp][1]], dtype=long)] + out_polygons[jj][jp:] + out_polygons[jj][:jp+1] + [numpy.array([long(round(xp)), out_polygons[jj][jp][1]], dtype=long)] + out_polygons[ii][ip:]
                    out_polygons[jj] = []

    ## Cut polygons with more than ``max_points`` points
    if max_points > 4:
        ii = 0
        while ii < len(out_polygons):
            if len(out_polygons[ii]) > max_points:
                #divisions = len(out_polygons[ii]) // max_points
                #min_x = out_polygons[ii][0][0]
                #max_x = out_polygons[ii][0][0]
                #min_y = out_polygons[ii][0][1]
                #max_y = out_polygons[ii][0][1]
                #for jj in range(1, len(out_polygons[ii])):
                #    if out_polygons[ii][jj][1] < min_y:
                #        min_y = out_polygons[ii][jj][1]
                #    elif out_polygons[ii][jj][1] > max_y:
                #        max_y = out_polygons[ii][jj][1]
                #    if out_polygons[ii][jj][0] < min_x:
                #        min_x = out_polygons[ii][jj][0]
                #    elif out_polygons[ii][jj][0] > max_x:
                #        max_x = out_polygons[ii][jj][0]
                #if max_y - min_y >= max_x - min_x:
                #    ## Cut in the x direction
                #    size = (max_y - min_y) / (divisions + 1.0)
                #    boundaries = [[(min_x, min_y + jj * size), (max_x, min_y + jj * size), (max_x, min_y + (jj + 1) * size), (min_x, min_y + (jj + 1) * size)] for jj in range(divisions + 1)]
                #else:
                #    ## Cut in the y direction
                #    size = (max_x - min_x) / (divisions + 1.0)
                #    boundaries = [[(min_x + jj * size, min_y), (min_x + jj * size, max_y), (min_x + (jj + 1) * size, max_y), (min_x + (jj + 1) * size, min_y)] for jj in range(divisions + 1)]
                #for jj in range(divisions + 1):
                #    out_polygons += clip([boundaries[jj], out_polygons[ii]], lambda p: p[0] and p[1], 1, max_points)
                pts0 = list(out_polygons[ii])
                pts1 = list(out_polygons[ii])
                pts0.sort(compare_pair_0)
                pts1.sort(compare_pair_1)
                if pts0[-1][0] - pts0[0][0] > pts1[-1][1] - pts1[0][1]:
                    ## Vertical cuts
                    chopped = chop(out_polygons[ii], (pts0[len(pts0) // 2][0] + pts0[len(pts0) // 2 + 1][0]) // 2, 0)
                else:
                    ## Horizontal cuts
                    chopped = chop(out_polygons[ii], (pts1[len(pts1) // 2][1] + pts1[len(pts1) // 2 + 1][1]) // 2, 1)
                out_polygons.pop(ii) #remove(out_polygons[ii])
                out_polygons += chopped[0]
                out_polygons += chopped[1]
            else:
                ii += 1

    ii = len(out_polygons) - 1
    while ii >= 0:
        if len(out_polygons[ii]) < 3:
            out_polygons.remove(out_polygons[ii])
        else:
            out_polygons[ii] = numpy.array(out_polygons[ii]) * precision
        ii -= 1
    return out_polygons
