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

import os
import struct
import datetime
import numpy
from boolean import clip, chop, compare_pair_0, compare_pair_1
from viewer import LayoutViewer, gds_image#, gds_view

__version__ = '0.2.8'
__doc__ = """
gdspy is a Python module that allows the creation of GDSII stream files.

Many features of the GDSII format are implemented, such as cell references
and arrays, but the support for fonts is quite limited. Text is only
available through polygonal objects.

If the Python Imaging Library is installed, it can be used to output the
geometry created to an image file.
"""

def _eight_byte_real(value):
    """
    Convert a number into the GDSII 8 byte real format.
    
    Parameters
    ----------
    value : number
        The number to be converted.
    
    Returns
    -------
    out : string
        The GDSII binary string that represents ``value``.
    """
    byte1 = 0
    byte2 = 0
    short3 = 0
    long4 = 0
    if value != 0:
        if value < 0:
            byte1 = 0x80
            value = -value
        exponent = int(numpy.floor(numpy.log2(value) * 0.25))
        mantissa = long(value * 16L**(14 - exponent))
        while mantissa >= 72057594037927936L:
            exponent += 1
            mantissa = long(value * 16L**(14 - exponent))
        byte1 += exponent + 64
        byte2 = (mantissa // 281474976710656L)
        short3 = (mantissa % 281474976710656L) // 4294967296L
        long4 = mantissa % 4294967296L
    return struct.pack(">HHL", byte1 * 256 + byte2, short3, long4)


def _eight_byte_real_to_float(value):
    """
    Convert a number from GDSII 8 byte real format to float.

    Parameters
    ----------
    value : string
        The GDSII binary string representation of the number.

    Returns
    -------
    out : float
        The number represented by ``value``.
    """
    short1, short2, long3 = struct.unpack('>HHL', value)
    exponent = (short1 & 0x7f00) // 256
    mantissa = (((short1 & 0x00ff) * 65536L + short2) * 4294967296L + long3) / 72057594037927936.0
    return (-1 if (short1 & 0x8000) else 1) * mantissa * 16L ** (exponent - 64)


class Polygon:
    """
    Polygonal geometric object.

    Parameters
    ----------
    layer : integer
        The GDSII layer number for this element.
    points : array-like[N][2]
        Coordinates of the vertices of the polygon.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).
    
    Notes
    -----
    The last point should not be equal to the first (polygons are
    automatically closed).
    
    The GDSII specification supports only a maximum of 199 vertices per
    polygon.

    Examples
    --------
    >>> triangle_pts = [(0, 40), (15, 40), (10, 50)]
    >>> triangle = gdspy.Polygon(1, triangle_pts)
    >>> myCell.add(triangle)
    """
    
    def __init__(self, layer, points, datatype=0) :
        if len(points) > 199:
            print "gdspy - WARNING: a polygon with more than 199 points was created (not officially supported by the GDSII format)."
        self.layer = layer
        self.points = numpy.array(points)
        self.datatype = datatype

    def to_gds(self, multiplier): 
        """
        Convert this object to a GDSII element.
        
        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.
        
        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        data = struct.pack('>10h', 4, 0x0800, 6, 0x0D02, self.layer, 6, 0x0E02, self.datatype, 12 + 8 * len(self.points), 0x1003)
        for point in self.points:
            data += struct.pack('>2l', int(round(point[0] * multiplier)), int(round(point[1] * multiplier)))
        return data + struct.pack('>2l2h', int(round(self.points[0][0] * multiplier)), int(round(self.points[0][1] * multiplier)), 4, 0x1100)
    
    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this object.
        
        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.
        
        Returns
        -------
        out : ``Polygon``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        self.points = (self.points - c0) * ca + (self.points - c0)[:,::-1] * sa + c0
        return self

    def area(self, by_layer=False):
        """
        Calculate the total area of this object.
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary
            ``{layer: area}``.
        
        Returns
        -------
        out : number, dictionary
            Area of this object.
        """
        poly_area = 0
        for ii in range(1, len(self.points) - 1):
            poly_area += (self.points[0][0] - self.points[ii + 1][0]) * (self.points[ii][1] - self.points[0][1]) - (self.points[0][1] - self.points[ii + 1][1]) * (self.points[ii][0] - self.points[0][0])
        if by_layer:
            return {self.layer: 0.5 * abs(poly_area)}
        else:
            return 0.5 * abs(poly_area)

    def fracture(self, max_points=199):
        """
        Slice this polygon in the horizontal and vertical directions so
        that each resulting piece has at most ``max_points``.
        
        Parameters
        ----------
        max_points : integer
            Maximal number of points in each resulting polygon (must be
            greater than 4).
        
        Returns
        -------
        out : ``PolygonSet``
            Resulting polygons from the fracture operation.
        """
        out_polygons = [self.points]
        if max_points > 4:
            ii = 0
            while ii < len(out_polygons):
                if len(out_polygons[ii]) > max_points:
                    pts0 = list(out_polygons[ii])
                    pts1 = list(out_polygons[ii])
                    pts0.sort(compare_pair_0)
                    pts1.sort(compare_pair_1)
                    if pts0[-1][0] - pts0[0][0] > pts1[-1][1] - pts1[0][1]:
                        ## Vertical cuts
                        chopped = chop(out_polygons[ii], (pts0[len(pts0) // 2][0] + pts0[len(pts0) // 2 + 1][0]) / 2, 0)
                    else:
                        ## Horizontal cuts
                        chopped = chop(out_polygons[ii], (pts1[len(pts1) // 2][1] + pts1[len(pts1) // 2 + 1][1]) / 2, 1)
                    out_polygons.pop(ii)
                    out_polygons += chopped[0]
                    out_polygons += chopped[1]
                else:
                    ii += 1
        return PolygonSet(self.layer, out_polygons, self.datatype)

    def fillet(self, radius, points_per_2pi=128, max_points=199):
        """
        Round the corners of this polygon and fractures it into polygons
        with less vertices if necessary.
        
        Parameters
        ----------
        radius : number
            Radius of the corners.
        points_per_2pi : integer
            Number of vertices used to approximate a full circle.  The
            number of vertices in each corner of the polygon will be the
            fraction of this number corresponding to the angle encompassed
            by that corner with respect to 2 pi.
        max_points : integer
            Maximal number of points in each resulting polygon (must be
            greater than 4).
        
        Returns
        -------
        out : ``Polygon`` or ``PolygonSet``
            If no fracturing occurs, return this object; otherwise return a
            ``PolygonSet`` with the fractured result (this object will have
            more than ``max_points`` vertices).
        """
        two_pi = 2 * numpy.pi
        vec = self.points.astype(float) - numpy.roll(self.points, 1, 0)
        length = numpy.sqrt(numpy.sum(vec ** 2, 1))
        ii = numpy.flatnonzero(length)
        if len(ii) < len(length):
            self.points = self.points[ii]
            vec = self.points - numpy.roll(self.points, 1, 0)
            length = numpy.sqrt(numpy.sum(vec ** 2, 1))
        vec[:,0] /= length
        vec[:,1] /= length
        dvec = numpy.roll(vec, -1, 0) - vec
        norm = numpy.sqrt(numpy.sum(dvec ** 2, 1))
        ii = numpy.flatnonzero(norm)
        dvec[ii,0] /= norm[ii]
        dvec[ii,1] /= norm[ii]
        theta = numpy.arccos(numpy.sum(numpy.roll(vec, -1, 0) * vec, 1))
        ct = numpy.cos(theta * 0.5)
        tt = numpy.tan(theta * 0.5)

        new_points = []
        for ii in range(-1, len(self.points) - 1):
            if theta[ii] > 0:
                a0 = -vec[ii] * tt[ii] - dvec[ii] / ct[ii]
                a0 = numpy.arctan2(a0[1], a0[0])
                a1 = vec[ii + 1] * tt[ii] - dvec[ii] / ct[ii]
                a1 = numpy.arctan2(a1[1], a1[0])
                if a1 - a0 > numpy.pi:
                    a1 -= two_pi
                elif a1 - a0 < -numpy.pi:
                    a1 += two_pi
                n = max(int(numpy.ceil(abs(a1 - a0) / two_pi * points_per_2pi) + 0.5), 2)
                a = numpy.linspace(a0, a1, n)
                l = radius * tt[ii]
                if l > 0.49 * length[ii]:
                    r = 0.49 * length[ii] / tt[ii]
                    l = 0.49 * length[ii]
                else:
                    r = radius
                if l > 0.49 * length[ii + 1]:
                    r = 0.49 * length[ii + 1] / tt[ii]
                new_points += list(r * dvec[ii] / ct[ii] + self.points[ii] + numpy.vstack((r * numpy.cos(a), r * numpy.sin(a))).transpose())
            else:
                new_points.append(self.points[ii])
        
        self.points = numpy.array(new_points)
        if len(self.points) > max_points:
            return self.fracture()
        else:
            return self


class PolygonSet:
    """ 
    Set of polygonal objects.

    Parameters
    ----------
    layer : integer
        The GDSII layer number for this element.
    polygons : list of array-like[N][2]
        List containing the coordinates of the vertices of each polygon.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).
    
    Notes
    -----
    The last point should not be equal to the first (polygons are
    automatically closed).
    
    The GDSII specification supports only a maximum of 199 vertices per
    polygon.
    """

    def __init__(self, layer, polygons, datatype=0):
        self.layers = [layer] * len(polygons)
        self.datatypes = [datatype] * len(polygons)
        self.polygons = [None] * len(polygons)
        for i in range(len(polygons)):
            self.polygons[i] = numpy.array(polygons[i])
            if len(polygons[i]) > 199:
                print "gdspy - WARNING: a polygon with more than 199 points was created (not officially supported by the GDSII format)."

    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this object.
        
        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.
        
        Returns
        -------
        out : ``PolygonSet``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        self.polygons = [(points - c0) * ca + (points - c0)[:,::-1] * sa + c0 for points in self.polygons]
        return self

    def to_gds(self, multiplier):
        """
        Convert this object to a series of GDSII elements.
        
        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            elements.
        
        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        data = ''
        for ii in range(len(self.polygons)):
            data += struct.pack('>10h', 4, 0x0800, 6, 0x0D02, self.layers[ii], 6, 0x0E02, self.datatypes[ii], 12 + 8 * len(self.polygons[ii]), 0x1003)
            for point in self.polygons[ii]:
                data += struct.pack('>2l', int(round(point[0] * multiplier)), int(round(point[1] * multiplier)))
            data += struct.pack('>2l2h', int(round(self.polygons[ii][0][0] * multiplier)), int(round(self.polygons[ii][0][1] * multiplier)), 4, 0x1100)
        return data

    def area(self, by_layer=False):
        """
        Calculate the total area of the path(s).
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary
            C{{layer: area}}.
        
        Returns
        -------
        out : number, dictionary
            Area of this object.
        """
        if by_layer:
            path_area = {}
            for jj in range(len(self.polygons)):
                poly_area = 0
                for ii in range(1, len(self.polygons[jj]) - 1):
                    poly_area += (self.polygons[jj][0][0] - self.polygons[jj][ii + 1][0]) * (self.polygons[jj][ii][1] - self.polygons[jj][0][1]) - (self.polygons[jj][0][1] - self.polygons[jj][ii + 1][1]) * (self.polygons[jj][ii][0] - self.polygons[jj][0][0])
                if path_area.has_key(self.layers[jj]):
                    path_area[self.layers[jj]] += 0.5 * abs(poly_area)
                else:
                    path_area[self.layers[jj]] = 0.5 * abs(poly_area)
        else:
            path_area = 0
            for points in self.polygons:
                poly_area = 0
                for ii in range(1, len(points) - 1):
                    poly_area += (points[0][0] - points[ii + 1][0]) * (points[ii][1] - points[0][1]) - (points[0][1] - points[ii + 1][1]) * (points[ii][0] - points[0][0])
                path_area += 0.5 * abs(poly_area)
        return path_area

    def fracture(self, max_points=199):
        """
        Slice these polygons in the horizontal and vertical directions so
        that each resulting piece has at most ``max_points``.  This
        opertaion occur in place.
        
        Parameters
        ----------
        max_points : integer
            Maximal number of points in each resulting polygon (must be
            greater than 4).

        Returns
        -------
        out : ``PolygonSet``
            This object.
        """
        if max_points > 4:
            ii = 0
            while ii < len(self.polygons):
                if len(self.polygons[ii]) > max_points:
                    pts0 = list(self.polygons[ii])
                    pts1 = list(self.polygons[ii])
                    pts0.sort(compare_pair_0)
                    pts1.sort(compare_pair_1)
                    if pts0[-1][0] - pts0[0][0] > pts1[-1][1] - pts1[0][1]:
                        ## Vertical cuts
                        chopped = chop(self.polygons[ii], (pts0[len(pts0) // 2][0] + pts0[len(pts0) // 2 + 1][0]) / 2, 0)
                    else:
                        ## Horizontal cuts
                        chopped = chop(self.polygons[ii], (pts1[len(pts1) // 2][1] + pts1[len(pts1) // 2 + 1][1]) / 2, 1)
                    self.polygons.pop(ii)
                    layer = self.layers.pop(ii)
                    datatype = self.datatypes.pop(ii)
                    self.polygons += [numpy.array(x) for x in chopped[0] + chopped[1]]
                    self.layers += [layer] * (len(chopped[0]) + len(chopped[1]))
                    self.datatypes += [datatype] * (len(chopped[0]) + len(chopped[1]))
                else:
                    ii += 1
        return self

    def fillet(self, radius, points_per_2pi=128, max_points=199):
        """
        Round the corners of these polygons and fractures them into
        polygons with less vertices if necessary.
        
        Parameters
        ----------
        radius : number
            Radius of the corners.
        points_per_2pi : integer
            Number of vertices used to approximate a full circle.  The
            number of vertices in each corner of the polygon will be the
            fraction of this number corresponding to the angle encompassed
            by that corner with respect to 2 pi.
        max_points : integer
            Maximal number of points in each resulting polygon (must be
            greater than 4).
        
        Returns
        -------
        out : ``PolygonSet``
            This object.
        """
        two_pi = 2 * numpy.pi
        fracture = False

        for jj in range(len(self.polygons)):
            vec = self.polygons[jj].astype(float) - numpy.roll(self.polygons[jj], 1, 0)
            length = numpy.sqrt(numpy.sum(vec ** 2, 1))
            ii = numpy.flatnonzero(length)
            if len(ii) < len(length):
                self.polygons[jj] = self.polygons[jj][ii]
                vec = self.polygons[jj] - numpy.roll(self.polygons[jj], 1, 0)
                length = numpy.sqrt(numpy.sum(vec ** 2, 1))
            vec[:,0] /= length
            vec[:,1] /= length
            dvec = numpy.roll(vec, -1, 0) - vec
            norm = numpy.sqrt(numpy.sum(dvec ** 2, 1))
            ii = numpy.flatnonzero(norm)
            dvec[ii,0] /= norm[ii]
            dvec[ii,1] /= norm[ii]
            theta = numpy.arccos(numpy.sum(numpy.roll(vec, -1, 0) * vec, 1))
            ct = numpy.cos(theta * 0.5)
            tt = numpy.tan(theta * 0.5)

            new_points = []
            for ii in range(-1, len(self.polygons[jj]) - 1):
                if theta[ii] > 0:
                    a0 = -vec[ii] * tt[ii] - dvec[ii] / ct[ii]
                    a0 = numpy.arctan2(a0[1], a0[0])
                    a1 = vec[ii + 1] * tt[ii] - dvec[ii] / ct[ii]
                    a1 = numpy.arctan2(a1[1], a1[0])
                    if a1 - a0 > numpy.pi:
                        a1 -= two_pi
                    elif a1 - a0 < -numpy.pi:
                        a1 += two_pi
                    n = max(int(numpy.ceil(abs(a1 - a0) / two_pi * points_per_2pi) + 0.5), 2)
                    a = numpy.linspace(a0, a1, n)
                    l = radius * tt[ii]
                    if l > 0.49 * length[ii]:
                        r = 0.49 * length[ii] / tt[ii]
                        l = 0.49 * length[ii]
                    else:
                        r = radius
                    if l > 0.49 * length[ii + 1]:
                        r = 0.49 * length[ii + 1] / tt[ii]
                    new_points += list(r * dvec[ii] / ct[ii] + self.polygons[jj][ii] + numpy.vstack((r * numpy.cos(a), r * numpy.sin(a))).transpose())
                else:
                    new_points.append(self.polygons[jj][ii])
            self.polygons[jj] = numpy.array(new_points)
            if len(new_points) > max_points:
                fracture = True

        if fracture:
            self.fracture(max_points)
        return self


class Rectangle(Polygon):
    """
    Rectangular geometric object.

    Parameters
    ----------
    layer : integer
        The GDSII layer number for this element.
    point1 : array-like[2]
        Coordinates of a corner of the rectangle.
    point2 : array-like[2]
        Coordinates of the corner of the rectangle opposite to ``point1``.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Examples
    --------
    >>> rectangle = gdspy.Rectangle(1, (0, 0), (10, 20))
    >>> myCell.add(rectangle)
    """
    
    def __init__(self, layer, point1, point2, datatype=0):
        self.layer = layer
        self.points = numpy.array([[point1[0], point1[1]], [point1[0], point2[1]], [point2[0], point2[1]], [point2[0], point1[1]]])
        self.datatype = datatype


class Round(PolygonSet):
    """
    Circular geometric object.
    Represent a circle, a circular section, a ring or a ring section.

    Parameters
    ----------
    layer : integer
        The GDSII layer number for this element.
    center : array-like[2]
        Coordinates of the center of the circle/ring.
    radius : number
        Radius of the circle/outer radius of the ring.
    inner_radius : number
        Inner radius of the ring.
    initial_angle : number
        Initial angle of the circular/ring section (in *radians*).
    final_angle : number
        Final angle of the circular/ring section (in *radians*).
    number_of_points : integer
        number of vertices that form the object (polygonal
        approximation).
    max_points : integer
        if ``number_of_points > max_points``, the element will be
        fractured in smaller polygons with at most ``max_points`` each.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).
    
    Notes
    -----
    The GDSII specification supports only a maximum of 199 vertices per
    polygon.

    Examples
    --------
    >>> circle = gdspy.Round(2, (30, 5), 8)
    >>> ring = gdspy.Round(2, (50, 5), 8, inner_radius=5)
    >>> pie_slice = gdspy.Round(2, (30, 25), 8, initial_angle=0,
    ...                         final_angle=-5.0*numpy.pi/6.0)
    >>> arc = gdspy.Round(2, (50, 25), 8, inner_radius=5,
    ...                   initial_angle=-5.0*numpy.pi/6.0,
    ...                   final_angle=0)
    """
    
    def __init__(self, layer, center, radius, inner_radius=0, initial_angle=0, final_angle=0, number_of_points=199, max_points=199, datatype=0):
        pieces = int(numpy.ceil(number_of_points / float(max_points)))
        number_of_points = number_of_points // pieces
        self.layers = [layer] * pieces
        self.datatypes = [datatype] * pieces
        self.polygons = [numpy.zeros((number_of_points, 2)) for _ in range(pieces)]
        if final_angle == initial_angle and pieces > 1:
            final_angle += 2 * numpy.pi
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        for ii in range(pieces):
            if angles[ii+1] == angles[ii]:
                if inner_radius <= 0:
                    angle = numpy.arange(number_of_points) * 2.0 * numpy.pi / number_of_points
                    self.polygons[ii][:,0] = numpy.cos(angle)
                    self.polygons[ii][:,1] = numpy.sin(angle)
                    self.polygons[ii] = self.polygons[ii] * radius + numpy.array(center)
                else:
                    n2 = number_of_points // 2
                    n1 = number_of_points - n2
                    angle = numpy.arange(n1) * 2.0 * numpy.pi / (n1 - 1.0)
                    self.polygons[ii][:n1,0] = numpy.cos(angle) * radius + center[0]
                    self.polygons[ii][:n1,1] = numpy.sin(angle) * radius + center[1]
                    angle = numpy.arange(n2) * -2.0 * numpy.pi / (n2 - 1.0)
                    self.polygons[ii][n1:,0] = numpy.cos(angle) * inner_radius + center[0]
                    self.polygons[ii][n1:,1] = numpy.sin(angle) * inner_radius + center[1]
            else:
                if inner_radius <= 0:
                    angle = numpy.linspace(angles[ii], angles[ii+1], number_of_points - 1)
                    self.polygons[ii][1:,0] = numpy.cos(angle)
                    self.polygons[ii][1:,1] = numpy.sin(angle)
                    self.polygons[ii] = self.polygons[ii] * radius + numpy.array(center)
                else:
                    n2 = number_of_points // 2
                    n1 = number_of_points - n2
                    angle = numpy.linspace(angles[ii], angles[ii+1], n1)
                    self.polygons[ii][:n1,0] = numpy.cos(angle) * radius + center[0]
                    self.polygons[ii][:n1,1] = numpy.sin(angle) * radius + center[1]
                    angle = numpy.linspace(angles[ii+1], angles[ii], n2)
                    self.polygons[ii][n1:,0] = numpy.cos(angle) * inner_radius + center[0]
                    self.polygons[ii][n1:,1] = numpy.sin(angle) * inner_radius + center[1]


class Text(PolygonSet):
    """
    Polygonal text object.
    
    Each letter is formed by a series of polygons.

    Parameters
    ----------
    layer : integer
        The GDSII layer number for these elements.
    text : string
        The text to be converted in geometric objects.
    size : number
        Base size of each character.
    position : array-like[2]
        Text position (lower left corner).
    horizontal : bool
        If ``True``, the text is written from left to right; if
        ``False``, from top to bottom.
    angle : number
        The angle of rotation of the text.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Examples
    --------
    >>> text = gdspy.Text(8, 'Sample text', 20, (-10, -100))
    >>> myCell.add(text)
    """
    __font = {
        '!':[[(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 4), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (2, 6), (2, 5)]],
        '"':[[(1, 7), (2, 7), (2, 8), (2, 9), (1, 9), (1, 8)], [(3, 7), (4, 7), (4, 8), (4, 9), (3, 9), (3, 8)]],
        '#':[[(0, 3), (1, 3), (1, 2), (2, 2), (2, 3), (2, 4), (2, 5), (3, 5), (3, 4), (2, 4), (2, 3), (3, 3), (3, 2), (4, 2), (4, 3), (5, 3), (5, 4), (4, 4), (4, 5), (5, 5), (5, 6), (4, 6), (4, 7), (3, 7), (3, 6), (2, 6), (2, 7), (1, 7), (1, 6), (0, 6), (0, 5), (1, 5), (1, 4), (0, 4)]],
        '$':[[(0, 2), (1, 2), (2, 2), (2, 1), (3, 1), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (4, 4), (4, 5), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (5, 8), (4, 8), (3, 8), (3, 9), (2, 9), (2, 8), (1, 8), (1, 7), (2, 7), (2, 6), (1, 6), (1, 5), (2, 5), (2, 4), (2, 3), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(4, 3), (5, 3), (5, 4), (4, 4)]],
        '%':[[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (2, 7), (2, 8), (2, 9), (1, 9), (0, 9), (0, 8)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (3, 4), (3, 3)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        '&':[[(0, 3), (1, 3), (1, 4), (1, 5), (0, 5), (0, 4)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (1, 6)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 5), (5, 5), (5, 6), (4, 6)]],
        '\'':[[(2, 7), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8)]],
        '(':[[(1, 4), (2, 4), (2, 5), (2, 6), (2, 7), (1, 7), (1, 6), (1, 5)], [(2, 3), (3, 3), (3, 4), (2, 4)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(3, 8), (4, 8), (4, 9), (3, 9)]],
        ')':[[(1, 2), (2, 2), (2, 3), (1, 3)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 3), (3, 3), (3, 4), (2, 4)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 4), (4, 4), (4, 5), (4, 6), (4, 7), (3, 7), (3, 6), (3, 5)]],
        '*':[[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 4), (1, 4), (1, 3), (2, 3), (2, 2), (3, 2), (3, 3), (4, 3), (4, 4), (5, 4), (5, 5), (4, 5), (4, 6), (3, 6), (3, 7), (2, 7), (2, 6), (1, 6), (1, 5), (0, 5)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        '+':[[(0, 4), (1, 4), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (4, 4), (5, 4), (5, 5), (4, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6), (2, 5), (1, 5), (0, 5)]],
        ',':[[(1, 0), (2, 0), (2, 1), (1, 1)], [(2, 1), (3, 1), (3, 2), (3, 3), (2, 3), (2, 2)]],
        '-':[[(0, 4), (1, 4), (2, 4), (3, 4), (4, 4), (5, 4), (5, 5), (4, 5), (3, 5), (2, 5), (1, 5), (0, 5)]],
        '.':[[(2, 2), (3, 2), (3, 3), (2, 3)]],
        '/':[[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        '0':[[(0, 3), (1, 3), (1, 4), (2, 4), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 5), (4, 4), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (3, 7)]],
        '1':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (1, 8), (1, 7), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (1, 3)]],
        '2':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (1, 8), (0, 8)], [(1, 4), (2, 4), (3, 4), (3, 5), (2, 5), (1, 5)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '3':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '4':[[(0, 4), (1, 4), (2, 4), (3, 4), (3, 3), (3, 2), (4, 2), (4, 3), (4, 4), (5, 4), (5, 5), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (3, 9), (2, 9), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (2, 5), (1, 5), (1, 6), (0, 6), (0, 5)], [(1, 6), (2, 6), (2, 7), (2, 8), (1, 8), (1, 7)]],
        '5':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        '6':[[(0, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        '7':[[(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 7), (4, 6), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (2, 5), (2, 4), (2, 3)], [(3, 5), (4, 5), (4, 6), (3, 6)]],
        '8':[[(0, 3), (1, 3), (1, 4), (1, 5), (0, 5), (0, 4)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '9':[[(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)]],
        ':':[[(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 5), (3, 5), (3, 6), (2, 6)]],
        ';':[[(1, 0), (2, 0), (2, 1), (1, 1)], [(2, 1), (3, 1), (3, 2), (3, 3), (2, 3), (2, 2)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        '<':[[(0, 5), (1, 5), (1, 6), (0, 6)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 3), (3, 3), (4, 3), (4, 4), (3, 4), (2, 4)], [(2, 7), (3, 7), (4, 7), (4, 8), (3, 8), (2, 8)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 8), (5, 8), (5, 9), (4, 9)]],
        '=':[[(0, 3), (1, 3), (2, 3), (3, 3), (4, 3), (5, 3), (5, 4), (4, 4), (3, 4), (2, 4), (1, 4), (0, 4)], [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5), (5, 5), (5, 6), (4, 6), (3, 6), (2, 6), (1, 6), (0, 6)]],
        '>':[[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 8), (1, 8), (1, 9), (0, 9)], [(1, 3), (2, 3), (3, 3), (3, 4), (2, 4), (1, 4)], [(1, 7), (2, 7), (3, 7), (3, 8), (2, 8), (1, 8)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 5), (5, 5), (5, 6), (4, 6)]],
        '?':[[(0, 7), (1, 7), (1, 8), (0, 8)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '@':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 4), (3, 4), (4, 4), (4, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6), (2, 5)], [(4, 5), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6)]],
        'A':[[(0, 2), (1, 2), (1, 3), (1, 4), (2, 4), (3, 4), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'B':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'C':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'D':[[(0, 2), (1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 4), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5)]],
        'E':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'F':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'G':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (3, 6), (2, 6), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'H':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'I':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (1, 3)]],
        'J':[[(0, 3), (1, 3), (1, 4), (0, 4)], [(0, 8), (1, 8), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (3, 4), (3, 3), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3)]],
        'K':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (2, 6), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(2, 6), (3, 6), (3, 7), (2, 7)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 8), (5, 8), (5, 9), (4, 9)]],
        'L':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'M':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 7), (2, 8), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6)], [(3, 7), (4, 7), (4, 6), (4, 5), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (3, 8)]],
        'N':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 7), (2, 8), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6)], [(3, 4), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (3, 5)]],
        'O':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'P':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'Q':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4), (3, 4), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        'R':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (3, 4), (4, 4), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'S':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        'T':[[(0, 8), (1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)]],
        'U':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'V':[[(0, 5), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6)], [(1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 5), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6)]],
        'W':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (2, 3), (1, 3)], [(2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (2, 6), (2, 5), (2, 4)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'X':[[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        'Y':[[(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 5), (2, 5), (2, 6), (2, 7), (1, 7), (1, 6)], [(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (2, 5), (2, 4), (2, 3)], [(3, 5), (4, 5), (4, 6), (4, 7), (3, 7), (3, 6)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        'Z':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 7), (3, 7)]],
        '[':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (1, 3)]],
        '\\':[[(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)]],
        ']':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (3, 4), (3, 3), (2, 3), (1, 3)]],
        '^':[[(0, 6), (1, 6), (1, 7), (0, 7)], [(1, 7), (2, 7), (2, 8), (1, 8)], [(2, 8), (3, 8), (3, 9), (2, 9)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        '_':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)]],
        '`':[[(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 7), (3, 7), (3, 8), (2, 8)]],
        'a':[[(0, 3), (1, 3), (1, 4), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 4), (2, 4), (3, 4), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'b':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'c':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'd':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (3, 7), (2, 7), (1, 7), (1, 6), (2, 6), (3, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)]],
        'e':[[(0, 3), (1, 3), (1, 4), (2, 4), (3, 4), (4, 4), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'f':[[(0, 5), (1, 5), (1, 4), (1, 3), (1, 2), (2, 2), (2, 3), (2, 4), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (2, 7), (2, 8), (1, 8), (1, 7), (1, 6), (0, 6)], [(2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9)]],
        'g':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 0), (2, 0), (3, 0), (4, 0), (4, 1), (3, 1), (2, 1), (1, 1)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'h':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'i':[[(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (2, 7), (1, 7)], [(2, 8), (3, 8), (3, 9), (2, 9)]],
        'j':[[(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1)], [(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (2, 1), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (2, 7), (1, 7)], [(2, 8), (3, 8), (3, 9), (2, 9)]],
        'k':[[(0, 2), (1, 2), (1, 3), (1, 4), (2, 4), (3, 4), (3, 5), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        'l':[[(1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (1, 9)]],
        'm':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'n':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'o':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'p':[[(0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3), (0, 2), (0, 1)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'q':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (4, 0), (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'r':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (3, 6), (2, 6), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7)]],
        's':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 5), (1, 5), (1, 6), (0, 6)], [(1, 4), (2, 4), (3, 4), (4, 4), (4, 5), (3, 5), (2, 5), (1, 5)], [(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7)], [(4, 3), (5, 3), (5, 4), (4, 4)]],
        't':[[(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (1, 7)], [(3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3)]],
        'u':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'v':[[(0, 5), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6)], [(1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 5), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6)]],
        'w':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (2, 3), (1, 3)], [(2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (2, 6), (2, 5), (2, 4)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'x':[[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(1, 3), (2, 3), (2, 4), (1, 4)], [(1, 5), (2, 5), (2, 6), (1, 6)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        'y':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 0), (2, 0), (3, 0), (4, 0), (4, 1), (3, 1), (2, 1), (1, 1)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)]],
        'z':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (2, 4), (1, 4), (1, 3), (0, 3)], [(0, 6), (1, 6), (2, 6), (3, 6), (3, 5), (4, 5), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        '{':[[(1, 5), (2, 5), (2, 4), (2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (1, 6)], [(3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3)], [(3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9)]],
        '|':[[(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3)]],
        '}':[[(0, 2), (1, 2), (2, 2), (2, 3), (1, 3), (0, 3)], [(0, 8), (1, 8), (2, 8), (2, 9), (1, 9), (0, 9)], [(2, 3), (3, 3), (3, 4), (3, 5), (4, 5), (4, 6), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4)]],
        '~':[[(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]]}

    def __init__(self, layer, text, size, position=(0, 0), horizontal=True, angle=0, datatype=0) :
        self.polygons = []
        posX = 0
        posY = 0
        text_multiplier = size / 9.0
        if angle == 0:
            ca = 1
            sa = 0
        else:
            ca = numpy.cos(angle)
            sa = numpy.sin(angle)
        for jj in range(len(text)):
            if text[jj] == '\n':
                if horizontal:
                    posY -= 11
                    posX = 0
                else:
                    posX += 8
                    posY = 0
            elif text[jj] == '\t':
                if horizontal:
                    posX = posX + 32 - (posX + 8) % 32
                else:
                    posY = posY - 11 - (posY - 22) % 44
            else:
                if text[jj] in Text.__font.keys():
                    for p in Text.__font[text[jj]]:
                        polygon = p[:]
                        for ii in range(len(polygon)):
                            xp = text_multiplier * (posX + polygon[ii][0])
                            yp = text_multiplier * (posY + polygon[ii][1])
                            polygon[ii] = (position[0] + xp * ca - yp * sa, position[1] + xp * sa + yp * ca)
                        self.polygons.append(numpy.array(polygon))
                if horizontal:
                    posX += 8
                else:
                    posY -= 11
        self.layers = [layer] * len(self.polygons)
        self.datatypes = [datatype] * len(self.polygons) 


class Path(PolygonSet):
    """
    Series of geometric objects that form a path or a collection of
    parallel paths.

    Parameters
    ----------
    width : number
        The width of each path.
    initial_point : array-like[2]
        Starting position of the path.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number
        Distance between the centers of adjacent paths.

    Attributes
    ----------
    x : number
        Current position of the path in the x direction.
    y : number
        Current position of the path in the y direction.
    w : number
        *Half*-width of each path.
    n : integer
        Number of parallel paths.
    direction : {'+x', '-x', '+y', '-y'} or number
        Direction or angle (in *radians*) the path points to.
    distance : number
        Distance between the centers of adjacent paths.
    length : number
        Length of the central path axis. If only one path is created,
        this is the real length of the path.
    """
    
    def __init__(self, width, initial_point=(0, 0), number_of_paths=1, distance=0):
        self.x = initial_point[0]
        self.y = initial_point[1]
        self.w = width * 0.5
        self.n = number_of_paths
        self.direction = '+x'
        self.distance = distance
        self.length = 0.0
        self.polygons = []
        self.layers = []
        self.datatypes = []
        
    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this object.
        
        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.
        
        Returns
        -------
        out : ``Path``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        if self.direction.__class__ == ''.__class__:
            self.direction = {'+x': 0, '+y': 0.5, '-x': 1, '-y': -0.5}[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [(points - c0) * ca + (points - c0)[:,::-1] * sa + c0 for points in self.polygons]
        return self

    def segment(self, layer, length, direction=None, final_width=None, final_distance=None, datatype=0) :
        """
        Add a straight section to the path.
        
        Parameters
        ----------
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        length : number
            Length of the section to add.
        direction : {'+x', '-x', '+y', '-y'} or number
            Direction or angle (in *radians*) of rotation of the
            segment.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255). If the number of datatypes in the list is less
            than the number of paths, the list is repeated.
        
        Returns
        -------
        out : ``Path``
            This object.
        """
        if direction == None:
            direction = self.direction
        else:
            self.direction = direction
        if direction == '+x':
            ca = 1
            sa = 0
        elif direction == '-x':
            ca = -1
            sa = 0
        elif direction == '+y':
            ca = 0
            sa = 1
        elif direction == '-y':
            ca = 0
            sa = -1
        else:
            ca = numpy.cos(direction)
            sa = numpy.sin(direction)
        old_x = self.x
        old_y = self.y
        self.x += length*ca
        self.y += length*sa
        old_w = self.w
        old_distance = self.distance
        if not final_width == None:
            self.w = final_width * 0.5
        if not final_distance == None:
            self.distance = final_distance
        if (self.w != 0) or (old_w != 0):
            for ii in range(self.n):
                d0 = ii * self.distance - (self.n - 1) * self.distance * 0.5
                old_d0 = ii * old_distance - (self.n - 1) * old_distance * 0.5
                self.polygons.append(numpy.array([(old_x+(old_d0-old_w)*sa, old_y-(old_d0-old_w)*ca), (old_x+(old_d0+old_w)*sa, old_y-(old_d0+old_w)*ca), (self.x+(d0+self.w)*sa, self.y-(d0+self.w)*ca), (self.x+(d0-self.w)*sa, self.y-(d0-self.w)*ca)]))
                if self.w == 0:
                    self.polygons[-1] = self.polygons[-1][:-1,:]
                if old_w == 0:
                    self.polygons[-1] = self.polygons[-1][1:,:]
            self.length += abs(length)
            if (layer.__class__ == [].__class__):
                self.layers += (layer * (self.n // len(layer) + 1))[:self.n]
            else:
                self.layers += [layer] * self.n
            if (datatype.__class__ == [].__class__) :
                self.datatypes += (datatype * (self.n // len(datatype) + 1))[ :self.n]
            else:
                self.datatypes += [datatype] * self.n
        return self

    def arc(self, layer, radius, initial_angle, final_angle, number_of_points=199, max_points=199, final_width=None, final_distance=None, datatype=0) :
        """
        Add a curved section to the path.
        
        Parameters
        ----------
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        radius : number
            Central radius of the section.
        initial_angle : number
            Initial angle of the curve (in *radians*).
        final_angle : number
            Final angle of the curve (in *radians*).
        number_of_points : integer
            Number of vertices that form the object (polygonal
            approximation).
        max_points : integer
            if ``number_of_points > max_points``, the element will be
            fractured in smaller polygons with at most ``max_points`` each.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255). If the number of datatypes in the list is less
            than the number of paths, the list is repeated.
                    
        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
        """
        cx = self.x - radius*numpy.cos(initial_angle)
        cy = self.y - radius*numpy.sin(initial_angle)
        self.x = cx + radius*numpy.cos(final_angle)
        self.y = cy + radius*numpy.sin(final_angle)
        if final_angle > initial_angle:
            self.direction = final_angle + numpy.pi * 0.5
        else:
            self.direction = final_angle - numpy.pi * 0.5
        old_w = self.w
        old_distance = self.distance
        if not final_width == None:
            self.w = final_width * 0.5
        if not final_distance == None:
            self.distance = final_distance
        pieces = int(numpy.ceil(number_of_points / float(max_points)))
        number_of_points = number_of_points // pieces
        widths = numpy.linspace(old_w, self.w, pieces + 1)
        distances = numpy.linspace(old_distance, self.distance, pieces + 1)
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        if (self.w != 0) or (old_w != 0):
            for jj in range(pieces):
                for ii in range(self.n):
                    self.polygons.append(numpy.zeros((number_of_points, 2)))
                    r0 = radius + ii * distances[jj+1] - (self.n - 1) * distances[jj+1] * 0.5
                    old_r0 = radius + ii * distances[jj] - (self.n - 1) * distances[jj] * 0.5
                    pts2 = number_of_points // 2
                    pts1 = number_of_points - pts2
                    ang = numpy.linspace(angles[jj], angles[jj+1], pts1)
                    rad = numpy.linspace(old_r0 + widths[jj], r0 + widths[jj+1], pts1)
                    self.polygons[-1][:pts1,0] = numpy.cos(ang) * rad + cx
                    self.polygons[-1][:pts1,1] = numpy.sin(ang) * rad + cy
                    if widths[jj+1] == 0:
                        pts1 -= 1
                        pts2 += 1
                    if widths[jj] == 0:
                        self.polygons[-1][:pts1-1,:] = numpy.array(self.polygons[-1][1:pts1,:])
                        pts1 -= 1
                        pts2 += 1
                    ang = numpy.linspace(angles[jj+1], angles[jj], pts2)
                    rad = numpy.linspace(r0 - widths[jj+1], old_r0 - widths[jj], pts2)
                    self.polygons[-1][pts1:,0] = numpy.cos(ang) * rad + cx
                    self.polygons[-1][pts1:,1] = numpy.sin(ang) * rad + cy
                self.length += abs((angles[jj+1] - angles[jj])*radius)
                if (layer.__class__ == [].__class__):
                    self.layers += (layer * (self.n // len(layer) + 1))[:self.n]
                else:
                    self.layers += [layer] * self.n
                if (datatype.__class__ == [].__class__) :
                    self.datatypes += (datatype * (self.n // len(datatype) + 1))[ :self.n]
                else:
                    self.datatypes += [datatype] * self.n
        return self

    def turn(self, layer, radius, angle, number_of_points=199, max_points=199, final_width=None, final_distance=None, datatype=0):
        """
        Add a curved section to the path.
        
        Parameters
        ----------
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        radius : number
            Central radius of the section.
        angle : {'r', 'l', 'rr', 'll'} or number
            Angle (in *radians*) of rotation of the path. The values 'r'
            and 'l' represent 90-degree turns cw and ccw, respectively;
            the values 'rr' and 'll' represent analogous 180-degree
            turns.
        number_of_points : integer
            Number of vertices that form the object (polygonal
            approximation).
        max_points : integer
            if ``number_of_points > max_points``, the element will be
            fractured in smaller polygons with at most ``max_points`` each.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255). If the number of datatypes in the list is less
            than the number of paths, the list is repeated.
                    
        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
        """
        exact = True
        if angle == 'r':
            delta_i = 0.5 * numpy.pi
            delta_f = 0
        elif angle == 'rr':
            delta_i = 0.5 * numpy.pi
            delta_f = -delta_i
        elif angle == 'l':
            delta_i = -0.5 * numpy.pi
            delta_f = 0
        elif angle == 'll':
            delta_i = -0.5 * numpy.pi
            delta_f = -delta_i
        elif angle < 0:
            exact = False
            delta_i = 0.5 * numpy.pi
            delta_f = delta_i + angle
        else:
            exact = False
            delta_i = -0.5 * numpy.pi
            delta_f = delta_i + angle
        if self.direction == '+x':
            self.direction = 0
        elif self.direction == '-x':
            self.direction = numpy.pi
        elif self.direction == '+y':
            self.direction = 0.5 * numpy.pi
        elif self.direction == '-y':
            self.direction = -0.5 * numpy.pi
        elif exact:
            exact = False
        self.arc(layer, radius, self.direction + delta_i, self.direction + delta_f, number_of_points, max_points, final_width, final_distance, datatype)
        if exact:
            self.direction = ['+x', '+y', '-x', '-y'][int(round(self.direction / (0.5 * numpy.pi))) % 4]
        return self

    def parametric(self, layer, curve_function, curve_derivative=None, number_of_evaluations=99, max_points=199, final_width=None, final_distance=None, datatype=0):
        """
        Add a parametric curve to the path.
        
        Parameters
        ----------
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        curve_function : function
            Function that defines the curve. Must be a function of one
            argument (that varies from 0 to 1) that returns a 2-element
            list, tuple or array (x, y).
        curve_derivative : function
            If set, it should be the derivative of the curve function.
            Must be a function of one argument (that varies from 0 to 1)
            that returns a 2-element list, tuple or array (x,y). If
            ``None``, the derivative will be calculated numerically.
        number_of_evaluations : integer
            Number of points where the curve function will be evaluated.
            The final segment will have twice this number of points.
        max_points : integer
            if ``2 * number_of_evaluations > max_points``, the element
            will be fractured in smaller polygons with at most
            ``max_points`` each.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255). If the number of datatypes in the list is less
            than the number of paths, the list is repeated.
        
        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The norm of the vector returned by ``curve_derivative`` is not
        important. Only the direction is used.
        
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.

        Examples
        --------
        >>> def my_parametric_curve(t):
        ...     return (2**t, t**2)
        >>> def my_parametric_curve_derivative(t):
        ...     return (0.69315 * 2**t, 2 * t)
        >>> my_path.parametric(1, my_parametric_curve,
        ...                    my_parametric_curve_derivative)
        """
        old_w = self.w
        old_distance = self.distance
        if not final_width == None:
            self.w = final_width * 0.5
        if not final_distance == None:
            self.distance = final_distance
        if curve_derivative == None:
            def numerical_diff(t):
                delta = 0.5 / (number_of_evaluations - 1.0)
                if t == 0:
                    x0,y0 = curve_function(0)
                    x1,y1 = curve_function(delta)
                elif t == 1:
                    x0,y0 = curve_function(1 - delta)
                    x1,y1 = curve_function(1)
                else:
                    x0,y0 = curve_function(t - delta)
                    x1,y1 = curve_function(t + delta)
                return (x1 - x0, y1 - y0)
            curve_derivative = numerical_diff
        pieces = int(numpy.ceil(2 * number_of_evaluations / float(max_points)))
        number_of_evaluations = number_of_evaluations // pieces
        widths = numpy.linspace(old_w, self.w, pieces + 1)
        distances = numpy.linspace(old_distance, self.distance, pieces + 1)
        boundaries = numpy.linspace(0, 1, pieces + 1)
        if (self.w != 0) or (old_w != 0):
            for kk in range(pieces):
                x0 = numpy.zeros(number_of_evaluations)
                y0 = numpy.zeros(number_of_evaluations)
                dx = numpy.zeros(number_of_evaluations)
                dy = numpy.zeros(number_of_evaluations)
                for jj in range(number_of_evaluations):
                    u = boundaries[kk] + (boundaries[kk+1] - boundaries[kk]) * jj / (number_of_evaluations - 1.0)
                    v = curve_function(u)
                    x0[jj] = v[0] + self.x
                    y0[jj] = v[1] + self.y
                    if jj > 0:
                        self.length += ((x0[jj] - x0[jj-1])**2 + (y0[jj] - y0[jj-1])**2)**0.5
                    v = curve_derivative(u)
                    w = (v[0]**2 + v[1]**2)**0.5
                    dx[jj] = -v[1] / w
                    dy[jj] = v[0] / w
                for ii in range(self.n):
                    p1 = [None] * number_of_evaluations
                    p2 = [None] * number_of_evaluations
                    for jj in range(number_of_evaluations):
                        d0 = distances[kk] + (distances[kk+1] - distances[kk]) * jj / (number_of_evaluations - 1.0)
                        d0 = ii * d0 - (self.n - 1) * d0 * 0.5
                        w = widths[kk] + (widths[kk+1] - widths[kk]) * jj / (number_of_evaluations - 1.0)
                        p1[jj] = (x0[jj] + (d0 + w) * dx[jj], y0[jj] + (d0 + w) * dy[jj])
                        p2[number_of_evaluations - jj - 1] = (x0[jj] + (d0 - w) * dx[jj], y0[jj] + (d0 - w) * dy[jj])
                    if widths[kk+1] == 0:
                        p2 = p2[1:]
                    if widths[kk] == 0:
                        p1 = p1[1:]
                    self.polygons.append(numpy.array(p1 + p2))
                if (layer.__class__ == [].__class__):
                    self.layers += (layer * (self.n // len(layer) + 1))[:self.n]
                else:
                    self.layers += [layer] * self.n
                if (datatype.__class__ == [].__class__) :
                    self.datatypes += (datatype * (self.n // len(datatype) + 1))[:self.n]
                else:
                    self.datatypes += [datatype] * self.n
            self.x = x0[-1]
            self.y = y0[-1]
            self.direction = numpy.arctan2(-dx[-1], dy[-1])
        return self

        
class L1Path(PolygonSet):
    """
    Series of geometric objects that form a path or a collection of
    parallel paths with Manhattan geometry.

    Parameters
    ----------
    layer : integer, list
        The GDSII layer numbers for the elements of each path. If
        the number of layers in the list is less than the number of
        paths, the list is repeated.
    initial_point : array-like[2]
        Starting position of the path.
    direction : {'+x', '+y', '-x', '-y'}
        Starting direction of the path.
    width : number
        The initial width of each path.
    length : array-like
        Lengths of each section to add.
    turn : array-like
        Direction to turn before each section. The sign indicate the
        turn direction (ccw is positive), and the modulus is a
        multiplicative factor for the path width after each turn. Must
        have 1 element less then ``length``.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number
        Distance between the centers of adjacent paths.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0
        and 255). If the number of datatypes in the list is less
        than the number of paths, the list is repeated.
    
    Returns
    -------
    out : ``L1Path``
        This object.

    Attributes
    ----------
    x : number
        Final position of the path in the x direction.
    y : number
        Final position of the path in the y direction.
    direction : {'+x', '-x', '+y', '-y'} or number
        Direction or angle (in *radians*) the path points to. The
        numerical angle is returned only after a rotation of the object.

    Examples
    --------
    >>> length = [10, 30, 15, 15, 15, 15, 10]
    >>> turn = [1, -1, -1, 3, -1, 1]
    >>> l1path = gdspy.L1Path(5, (0, 0), '+x', 2, length, turn)
    >>> myCell.add(l1path)
    """
    def __init__(self, layer, initial_point, direction, width, length, turn, number_of_paths=1, distance=0, max_points=199, datatype=0):
        if (layer.__class__ != [].__class__):
            layer = [layer]
        if (datatype.__class__ != [].__class__) :
            datatype = [datatype]
        layer = (layer * (number_of_paths // len(layer) + 1))[:number_of_paths]
        datatype = (datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]
        w = width * 0.5
        points = max_points // 2 - 1
        paths = [[[], []] for ii in range(number_of_paths)]
        self.polygons = []
        self.layers = []
        self.datatypes = []
        self.x = initial_point[0]
        self.y = initial_point[1]
        if direction == '+x':
            direction = 0
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((initial_point[0], d0 + initial_point[1] - w))
                paths[ii][1].append((initial_point[0], d0 + initial_point[1] + w))
        elif direction == '+y':
            direction = 1
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                paths[ii][0].append((d0 + initial_point[0] + w, initial_point[1]))
                paths[ii][1].append((d0 + initial_point[0] - w, initial_point[1]))
        elif direction == '-x':
            direction = 2
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance 
                paths[ii][0].append((initial_point[0], d0 + initial_point[1] + w))
                paths[ii][1].append((initial_point[0], d0 + initial_point[1] - w))
        elif direction == '-y':
            direction = 3
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((d0 + initial_point[0] - w, initial_point[1]))
                paths[ii][1].append((d0 + initial_point[0] + w, initial_point[1]))
        for jj in range(len(turn)):
            points -= 1
            if direction == 0:
                for ii in range(number_of_paths):
                    d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                    paths[ii][0].append((self.x + length[jj] - (d0 - w) * turn[jj], paths[ii][0][-1][1]))
                    paths[ii][1].append((self.x + length[jj] - (d0 + w) * turn[jj], paths[ii][1][-1][1]))
                self.x += length[jj]
            elif direction == 1:
                for ii in range(number_of_paths):
                    d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                    paths[ii][0].append((paths[ii][0][-1][0], self.y + length[jj] - (d0 - w) * turn[jj]))
                    paths[ii][1].append((paths[ii][1][-1][0], self.y + length[jj] - (d0 + w) * turn[jj]))
                self.y += length[jj]
            elif direction == 2:
                for ii in range(number_of_paths):
                    d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance 
                    paths[ii][0].append((self.x - length[jj] - (d0 + w) * turn[jj], paths[ii][0][-1][1]))
                    paths[ii][1].append((self.x - length[jj] - (d0 - w) * turn[jj], paths[ii][1][-1][1]))
                self.x -= length[jj]
            elif direction == 3:
                for ii in range(number_of_paths):
                    d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                    paths[ii][0].append((paths[ii][0][-1][0], self.y - length[jj] - (d0 + w) * turn[jj]))
                    paths[ii][1].append((paths[ii][1][-1][0], self.y - length[jj] - (d0 - w) * turn[jj]))
                self.y -= length[jj]
            if points == 0:
                for p in paths:
                    if direction % 2 == 0:
                        min_dist = 1e300
                        for x1 in [p[0][-2][0], p[1][-2][0]]:
                            for x2 in [p[0][-1][0], p[1][-1][0]]:
                                if abs(x1 - x2) < min_dist:
                                    x0 = 0.5 * (x1 + x2)
                                    min_dist = abs(x1 - x2)
                        p0 = (x0, p[0][-1][1])
                        p1 = (x0, p[1][-1][1])
                    else:
                        min_dist = 1e300
                        for y1 in [p[0][-2][1], p[1][-2][1]]:
                            for y2 in [p[0][-1][1], p[1][-1][1]]:
                                if abs(y1 - y2) < min_dist:
                                    y0 = 0.5 * (y1 + y2)
                                    min_dist = abs(y1 - y2)
                        p0 = (p[0][-1][0], y0)
                        p1 = (p[1][-1][0], y0)
                    self.polygons.append(numpy.array(p[0][:-1] + [p0, p1] + p[1][-2::-1]))
                    p[0] = [p0, p[0][-1]]
                    p[1] = [p1, p[1][-1]]
                self.layers += layer
                self.datatypes += datatype
                points = max_points // 2 - 2
            if turn[jj] > 0:
                direction = (direction + 1) % 4
            else:
                direction = (direction - 1) % 4
        if direction == 0:
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((self.x + length[-1], paths[ii][0][-1][1]))
                paths[ii][1].append((self.x + length[-1], paths[ii][1][-1][1]))
            self.x += length[-1]
        elif direction == 1:
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((paths[ii][0][-1][0], self.y + length[-1]))
                paths[ii][1].append((paths[ii][1][-1][0], self.y + length[-1]))
            self.y += length[-1]
        elif direction == 2:
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance 
                paths[ii][0].append((self.x - length[-1], paths[ii][0][-1][1]))
                paths[ii][1].append((self.x - length[-1], paths[ii][1][-1][1]))
            self.x -= length[-1]
        elif direction == 3:
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance 
                paths[ii][0].append((paths[ii][0][-1][0], self.y - length[-1]))
                paths[ii][1].append((paths[ii][1][-1][0], self.y - length[-1]))
            self.y -= length[jj]
        self.direction = ['+x', '+y', '-x', '-y'][direction]
        self.polygons += [numpy.array(p[0] + p[1][::-1]) for p in paths]
        self.layers += layer
        self.datatypes += datatype

    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this object.
        
        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.
        
        Returns
        -------
        out : ``L1Path``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        self.direction = {'+x': 0, '+y': 0.5, '-x': 1, '-y': -0.5}[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [(points - c0) * ca + (points - c0)[:,::-1] * sa + c0 for points in self.polygons]
        return self


class PolyPath(PolygonSet):
    """
    Series of geometric objects that form a polygonal path or a collection of
    parallel polygonal paths.

    Parameters
    ----------
    layer : integer, list
        The GDSII layer numbers for the elements of each path. If
        the number of layers in the list is less than the number of
        paths, the list is repeated.
    points : array-like[N][2]
        Endpoints of each path segment.
    width : number or array-like[N]
        Width of the path. If an array is given, width at each endpoint.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number or array-like[N]
        Distance between the centers of adjacent paths. If na array is
        given, distance at each endpoint.
    corners : positive integer
        Type of the joins: 0 - miter join, 1 - bevel join
    max_points : integer
        The paths will be fractured in polygons with at most
        ``max_points``.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0
        and 255). If the number of datatypes in the list is less
        than the number of paths, the list is repeated.
    
    Returns
    -------
    out : ``PolyPath``
        This object.

    Notes
    -----
    The bevel joint will give  strange results if the number of paths is
    greater than 1.
    """
    def __init__(self, layer, points, width, number_of_paths=1, distance=0, corners=0, max_points=199, datatype=0):
        if (layer.__class__ != [].__class__):
            layer = [layer]
        if (datatype.__class__ != [].__class__) :
            datatype = [datatype]
        if width.__class__ is int or width.__class__ is long or width.__class__ is float:
            width = numpy.array([width * 0.5])
        else:
            width = numpy.array(width) * 0.5
        len_w = len(width)
        if distance.__class__ is int or distance.__class__ is long or distance.__class__ is float:
            distance = (distance,)
        len_d = len(distance)
        points = numpy.array(points)
        v = points[1,:] - points[0,:]
        v = numpy.array((-v[1], v[0])) / numpy.sqrt(numpy.sum(v * v))
        d0 = 0.5 * (number_of_paths - 1) * distance[0]
        d1 = 0.5 * (number_of_paths - 1) * distance[1 % len_d]
        self.polygons = []
        self.layers = []
        self.datatypes = []
        paths = [[[points[0,:] + (ii * distance[0] - d0 - width[0]) * v], [points[0,:] + (ii * distance[0] - d0 + width[0]) * v]] for ii in range(number_of_paths)]
        p1 = [(points[1,:] + (ii * distance[1 % len_d] - d1 - width[1 % len_w]) * v, points[1,:] + (ii * distance[1 % len_d] - d1 + width[1 % len_w]) * v) for ii in range(number_of_paths)]
        for jj in range(1, points.shape[0] - 1):
            j0d = jj % len_d
            j0w = jj % len_w
            j1d = (jj + 1) % len_d
            j1w = (jj + 1) % len_w
            v = points[jj + 1,:] - points[jj,:]
            v = numpy.array((-v[1], v[0])) / numpy.sqrt(numpy.sum(v * v))
            d0 = d1
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            p0 = p1
            p1 = []
            pp = []
            for ii in range(number_of_paths):
                pp.append((points[jj,:] + (ii * distance[j0d] - d0 - width[j0w]) * v, points[jj,:] + (ii * distance[j0d] - d0 + width[j0w]) * v))
                p1.append((points[jj + 1,:] + (ii * distance[j1d] - d1 - width[j1w]) * v, points[jj + 1,:] + (ii * distance[j1d] - d1 + width[j1w]) * v))
                for kk in (0, 1):
                    p0m = paths[ii][kk][-1] - p0[ii][kk]
                    p1p = pp[ii][kk] - p1[ii][kk]
                    vec = p0m[0] * p1p[1] - p1p[0] * p0m[1]
                    if abs(vec) > 1e-30:
                        p = numpy.array((1, -1)) * (p0m * p1p[::-1] * p1[ii][kk] - p1p * p0m[::-1] * p0[ii][kk] + p0m * p1p * (p0[ii][kk][::-1] - p1[ii][kk][::-1])) / vec
                        if corners > 0 and numpy.sum((p - pp[ii][kk]) * p1p) > 0 and numpy.sum((p - p0[ii][kk]) * p0m) < 0:
                            paths[ii][kk].append(p0[ii][kk])
                            paths[ii][kk].append(pp[ii][kk])
                        else:
                            paths[ii][kk].append(p)
                if len(paths[ii][0]) + len(paths[ii][1]) + 3 > max_points:
                    if numpy.sum((paths[ii][0][0] - paths[ii][1][0]) ** 2) == 0:
                        paths[ii][1] = paths[ii][1][1:]
                    if numpy.sum((paths[ii][0][-1] - paths[ii][1][-1]) ** 2) == 0:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][-2::-1]))
                    else:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][::-1]))
                    paths[ii][0] = paths[ii][0][-1:]
                    paths[ii][1] = paths[ii][1][-1:]
                    self.layers.append(layer[ii % len(layer)])
                    self.datatypes.append(datatype[ii % len(datatype)])
        for ii in range(number_of_paths):
            if numpy.sum((paths[ii][0][0] - paths[ii][1][0]) ** 2) == 0:
                paths[ii][1] = paths[ii][1][1:]
            if numpy.sum((p1[ii][0] - p1[ii][1]) ** 2) != 0:
                paths[ii][0].append(p1[ii][0])
            paths[ii][1].append(p1[ii][1])
        self.polygons += [numpy.array(p[0] + p[1][::-1]) for p in paths]
        self.layers += (layer * (number_of_paths // len(layer) + 1))[:number_of_paths]
        self.datatypes += (datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]


class Cell:
    """
    Collection of elements, both geometric objects and references to
    other cells.
    
    Parameters
    ----------
    name : string
        The name of the cell.
    exclude_from_global : bool
        If ``True``, the cell will not be included in the global list of cells
        maintained by ``gdspy``.
    """
    
    cell_list = {}
    """
    Dictionary containing all cells created, indexed by name.  This
    dictionary is updated automatically whenever a new ``Cell`` object is
    created without the ``exclude_from_global`` flag.
    """
    
    def __init__(self, name, exclude_from_global=False):
        self.name = name
        self.elements = []
        if name in Cell.cell_list:
            raise ValueError("gdspy - A cell named %s has already been created." % (name,))
        if not exclude_from_global:
            Cell.cell_list[name] = self

    def to_gds(self, multiplier):
        """
        Convert this cell to a GDSII structure.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            structure.
        
        
        Returns
        -------
        out : string
            The GDSII binary string that represents this cell.
        """
        now = datetime.datetime.today()
        name = self.name
        if len(name)%2 != 0:
            name = name + '\0'
        data = struct.pack('>16h', 28, 0x0502, now.year, now.month, now.day, now.hour, now.minute, now.second, now.year, now.month, now.day, now.hour, now.minute, now.second, 4 + len(name), 0x0606) + name
        for element in self.elements:
            data += element.to_gds(multiplier)
        return data + struct.pack('>2h', 4, 0x0700)
        
    def copy(self, name, exclude_from_global=False):
        """
        Creates a copy of this cell.

        Parameters
        ----------
        name : string
            The name of the cell.
        
        
        Returns
        -------
        out : ``Cell``
            The new copy of this cell.
        """
        new_cell = Cell(name, exclude_from_global)
        new_cell.elements = list(self.elements)
        return new_cell

    def add(self, element):
        """
        Add a new element or list of elements to this cell.
        
        Parameters
        ----------
        element : object, list
            The element or list of elements to be inserted in this cell.
        
        Returns
        -------
        out : ``Cell``
            This cell.
        """
        if (element.__class__ == [].__class__):
            self.elements += element
        else:
            self.elements.append(element)
        return self
    
    def area(self, by_layer=False):
        """
        Calculate the total area of the elements on this cell, including
        cell references and arrays.
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary with the areas
            of each individual layer.
        
        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if by_layer:
            cell_area = {}
            for element in self.elements:
                element_area = element.area(True)
                for ll in element_area.keys():
                    if cell_area.has_key(ll):
                        cell_area[ll] += element_area[ll]
                    else:
                        cell_area[ll] = element_area[ll]
        else:
            cell_area = 0
            for element in self.elements:
                cell_area += element.area()
        return cell_area

    def get_polygons(self, by_layer=False):
        """
        Returns a list of polygons in this cell.
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary with the
            polygons of each individual layer.
        
        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons in each
            layer (if ``by_layer`` is ``True``).
        """
        if by_layer:
            polygons = {}
            for element in self.elements:
                if isinstance(element, Polygon):
                    if polygons.has_key(element.layer):
                        polygons[element.layer].append(numpy.array(element.points))
                    else:
                        polygons[element.layer] = [numpy.array(element.points)]
                elif isinstance(element, PolygonSet):
                    for ii in range(len(element.polygons)):
                        if polygons.has_key(element.layers[ii]):
                            polygons[element.layers[ii]].append(numpy.array(element.polygons[ii]))
                        else:
                            polygons[element.layers[ii]] = [numpy.array(element.polygons[ii])]
                else:
                    cell_polygons = element.get_polygons(True)
                    for kk in cell_polygons.keys():
                        if polygons.has_key(kk):
                            polygons[kk] += cell_polygons[kk]
                        else:
                            polygons[kk] = cell_polygons[kk]
        else:
            polygons = []
            for element in self.elements:
                if isinstance(element, Polygon):
                    polygons.append(numpy.array(element.points))
                elif isinstance(element, PolygonSet):
                    for points in element.polygons:
                        polygons.append(numpy.array(points))
                else:
                    polygons += element.get_polygons()
        return polygons

    def get_dependencies(self):
        """
        Returns a list of the cells included in this cell as references.
        
        Returns
        -------
        out : list of ``Cell``
            List of the cells referenced by this cell.
        """
        dependencies = []
        for element in self.elements:
            if isinstance(element, CellReference) or isinstance(element, CellArray):
                dependencies.append(element.ref_cell)
        return dependencies

    def flatten(self, single_layer=None):
        """
        Flatten all ``CellReference`` and ``CellArray`` elements in this
        cell into real polygons, instead of references.
        
        Parameters
        ----------
        single_layer : integer or None
            If not ``None``, all polygons will be transfered to the
            layer indicated by this number.

        Returns
        -------
        out : ``Cell``
            This cell.
        """
        if single_layer is None:
            polyDic = self.get_polygons(True)
            self.elements = []
            for layer in polyDic.iterkeys():
                self.add(PolygonSet(layer, polyDic[layer]))
        else:
            polygons = self.get_polygons()
            self.elements = []
            self.add(PolygonSet(single_layer, polygons))
        return self
    

class CellReference:
    """
    Simple reference to an existing cell.

    Parameters
    ----------
    ref_cell : ``Cell`` or string
        The referenced cell or its name.
    origin : array-like[2]
        Position where the reference is inserted.
    rotation : number
        Angle of rotation of the reference (in *degrees*).
    magnification : number
        Magnification factor for the reference.
    x_reflection : bool
        If ``True``, the reference is reflected parallel to the x
        direction before being rotated.
    """

    def __init__(self, ref_cell, origin=(0, 0), rotation=None, magnification=None, x_reflection=False):
        self.origin = origin
        self.ref_cell = Cell.cell_list.get(ref_cell, ref_cell)
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
    
    def to_gds(self, multiplier):
        """
        Convert this object to a GDSII element.
        
        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.
        
        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        name = self.ref_cell.name
        if len(name)%2 != 0:
            name = name + '\0'
        data = struct.pack('>4h', 4, 0x0A00, 4 + len(name), 0x1206) + name
        if not (self.rotation is None) or not (self.magnification is None) or self.x_reflection:
            word = 0
            values = ''
            if self.x_reflection:
                word += 0x8000
            if not (self.magnification is None):
                word += 0x0004
                values += struct.pack('>2h', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                word += 0x0002
                values += struct.pack('>2h', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>2hH', 6, 0x1A01, word) + values
        return data + struct.pack('>2h2l2h', 12, 0x1003, int(round(self.origin[0] * multiplier)), int(round(self.origin[1] * multiplier)), 4, 0x1100)
    
    def area(self, by_layer=False):
        """
        Calculate the total area of the referenced cell with the
        magnification factor included.
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary with the areas
            of each individual layer.
        
        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if self.magnification is None:
            return self.ref_cell.area(by_layer)
        else:
            if by_layer:
                factor = self.magnification * self.magnification
                cell_area = self.ref_cell.area(True)
                for kk in cell_area.keys():
                    cell_area[kk] *= factor
                return cell_area
            else:
                return self.ref_cell.area() * self.magnification * self.magnification

    def get_polygons(self, by_layer=False):
        """
        Returns a list of polygons created by this reference.
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary with the
            polygons of each individual layer.
        
        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons in each
            layer (if ``by_layer`` is ``True``).
        """
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
        if by_layer:
            polygons = self.ref_cell.get_polygons(True)
            for kk in polygons.keys():
                for ii in range(len(polygons[kk])):
                    if self.x_reflection:
                        polygons[kk][ii] *= numpy.array([1, -1], dtype=int)
                    if self.magnification is not None:
                        polygons[kk][ii] *= numpy.array([self.magnification, self.magnification])
                    if self.rotation is not None:
                        polygons[kk][ii] = polygons[kk][ii] * ct + polygons[kk][ii][:,::-1] * numpy.array([-st, st])
                    if self.origin is not None:
                        polygons[kk][ii] = polygons[kk][ii] + numpy.array(self.origin)
        else:
            polygons = self.ref_cell.get_polygons()
            for ii in range(len(polygons)):
                if self.x_reflection:
                    polygons[ii] *= numpy.array([1, -1], dtype=int)
                if self.magnification is not None:
                    polygons[ii] *= numpy.array([self.magnification, self.magnification])
                if self.rotation is not None:
                    polygons[ii] = polygons[ii] * ct + polygons[ii][:,::-1] * numpy.array([-st, st])
                if self.origin is not None:
                    polygons[ii] = polygons[ii] + numpy.array(self.origin)
        return polygons
    

class CellArray:
    """
    Multiple references to an existing cell in an array format.

    Parameters
    ----------
    ref_cell : ``Cell`` or string
        The referenced cell or its name.
    columns : positive integer
        Number of columns in the array.
    rows : positive integer
        Number of columns in the array.
    spacing : array-like[2]
        distances between adjacent columns and adjacent rows.
    origin : array-like[2]
        Position where the cell is inserted.
    rotation : number
        Angle of rotation of the reference (in *degrees*).
    magnification : number
        Magnification factor for the reference.
    x_reflection : bool
        If ``True``, the reference is reflected parallel to the x
        direction before being rotated.
    """

    def __init__(self, ref_cell, columns, rows, spacing, origin=(0, 0), rotation=None, magnification=None, x_reflection=False):
        self.columns = columns
        self.rows = rows
        self.spacing = spacing
        self.origin = origin
        self.ref_cell = Cell.cell_list.get(ref_cell, ref_cell)
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
    
    def to_gds(self, multiplier):
        """
        Convert this object to a GDSII element.
        
        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.
        
        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        name = self.ref_cell.name
        if len(name)%2 != 0:
            name = name + '\0'
        data = struct.pack('>4h', 4, 0x0B00, 4 + len(name), 0x1206) + name
        x2 = self.origin[0] + self.columns * self.spacing[0]
        y2 = self.origin[1]
        x3 = self.origin[0]
        y3 = self.origin[1] + self.rows * self.spacing[1]
        if not (self.rotation is None) or not (self.magnification is None) or self.x_reflection:
            word = 0
            values = ''
            if self.x_reflection:
                word += 0x8000
                y3 = 2 * self.origin[1] - y3
            if not (self.magnification is None):
                word += 0x0004
                values += struct.pack('>2h', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                word += 0x0002
                sa = numpy.sin(self.rotation * numpy.pi / 180.0)
                ca = numpy.cos(self.rotation * numpy.pi / 180.0)
                tmp = (x2 - self.origin[0]) * ca - (y2 - self.origin[1]) * sa + self.origin[0]
                y2 = (x2 - self.origin[0]) * sa + (y2 - self.origin[1]) * ca + self.origin[1]
                x2 = tmp
                tmp = (x3 - self.origin[0]) * ca - (y3 - self.origin[1]) * sa + self.origin[0]
                y3 = (x3 - self.origin[0]) * sa + (y3 - self.origin[1]) * ca + self.origin[1]
                x3 = tmp
                values += struct.pack('>2h', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>2hH', 6, 0x1A01, word) + values
        return data + struct.pack('>6h6l2h', 8, 0x1302, self.columns, self.rows, 28, 0x1003, int(round(self.origin[0] * multiplier)), int(round(self.origin[1] * multiplier)), int(round(x2 * multiplier)), int(round(y2 * multiplier)), int(round(x3 * multiplier)), int(round(y3 * multiplier)), 4, 0x1100)

    def area(self, by_layer=False):
        """
        Calculate the total area of the cell array with the
        magnification factor included.
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary with the areas
            of each individual layer.
        
        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if self.magnification is None:
            factor = self.columns * self.rows
        else:
            factor = self.columns * self.rows * self.magnification * self.magnification
        if by_layer:
            cell_area = self.ref_cell.area(True)
            for kk in cell_area.keys():
                cell_area[kk] *= factor
            return cell_area
        else:
            return self.ref_cell.area() * factor

    def get_polygons(self, by_layer=False):
        """
        Returns a list of polygons created by this reference.
        
        Parameters
        ----------
        by_layer : bool
            If ``True``, the return value is a dictionary with the
            polygons of each individual layer.
        
        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons in each
            layer (if ``by_layer`` is ``True``).
        """
        if self.magnification is not None:
            mag = numpy.array([self.magnification, self.magnification])
        else:
            mag = numpy.ones(2)
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
        if by_layer:
            cell_polygons = self.ref_cell.get_polygons(True)
            polygons = {}
            for kk in cell_polygons.keys():
                polygons[kk] = []
                for ii in range(self.columns):
                    for jj in range(self.rows):
                        for points in cell_polygons[kk]:
                            polygons[kk].append(points * mag + numpy.array([self.spacing[0] * ii, self.spacing[1] * jj]))
                            if self.x_reflection:
                                polygons[kk][-1] *= numpy.array([1, -1], dtype=int)
                            if self.rotation is not None:
                                polygons[kk][-1] = polygons[kk][-1] * ct + polygons[kk][-1][:,::-1] * numpy.array([-st, st])
                            if self.origin is not None:
                                polygons[kk][-1] = polygons[kk][-1] + numpy.array(self.origin)
        else:
            cell_polygons = self.ref_cell.get_polygons()
            polygons = []
            for ii in range(self.columns):
                for jj in range(self.rows):
                    for points in cell_polygons:
                        polygons.append(points * mag + numpy.array([self.spacing[0] * ii, self.spacing[1] * jj]))
                        if self.x_reflection:
                            polygons[-1] *= numpy.array([1, -1], dtype=int)
                        if self.rotation is not None:
                            polygons[-1] = polygons[-1] * ct + polygons[-1][:,::-1] * numpy.array([-st, st])
                        if self.origin is not None:
                            polygons[-1] = polygons[-1] + numpy.array(self.origin)
        return polygons


class GdsImport:
    """
    Object used to import structures from a GDSII stream file.

    Parameters
    ----------
    infile : file or string
        GDSII stream file (or path) to be imported. It must be opened for
        reading in binary format.
    unit : number
        Unit (in *meters*) to use for the imported structures. If ``None``,
        the units used to create the GDSII file will be used.
    rename : dictionary
        Dictionary used to rename the imported cells. Keys and values must
        be strings.

    Attributes
    ----------
    cell_list : dictionary
        Dictionary will all imported cells, indexed by name.

    Notes
    -----
    Not all features from the GDSII specification are currently supported.
    A warning will be produced if any unsuported features are found in the
    imported file.

    Examples
    --------
    >>> gdsii = gdspy.GdsImport('gdspy-sample.gds')
    >>> for cell_name in gdsii.cell_list:
    ...     gdsii.extract(cell_name)
    """
    def __init__(self, infile, unit=None, rename={}):
        self.cell_list = {}
        self._incomplete = []

        _record_name = ('HEADER', 'BGNLIB', 'LIBNAME', 'UNITS', 'ENDLIB', 'BGNSTR', 'STRNAME', 'ENDSTR', 'BOUNDARY', 'PATH', 'SREF', 'AREF', 'TEXT', 'LAYER', 'DATATYPE', 'WIDTH', 'XY', 'ENDEL', 'SNAME', 'COLROW', 'TEXTNODE', 'NODE', 'TEXTTYPE', 'PRESENTATION', 'SPACING', 'STRING', 'STRANS', 'MAG', 'ANGLE', 'UINTEGER', 'USTRING', 'REFLIBS', 'FONTS', 'PATHTYPE', 'GENERATIONS', 'ATTRTABLE', 'STYPTABLE', 'STRTYPE', 'ELFLAGS', 'ELKEY', 'LINKTYPE', 'LINKKEYS', 'NODETYPE', 'PROPATTR', 'PROPVALUE', 'BOX', 'BOXTYPE', 'PLEX', 'BGNEXTN', 'ENDTEXTN', 'TAPENUM', 'TAPECODE', 'STRCLASS', 'RESERVED', 'FORMAT', 'MASK', 'ENDMASKS', 'LIBDIRSIZE', 'SRFNAME', 'LIBSECUR')
        _unused_records = (0x05, 0x00, 0x01, 0x02, 0x034, 0x38)

        if infile.__class__ == ''.__class__:
            infile = open(infile, 'rb')
            close = True
        else:
            close = False

        record = self._read_record(infile)
        args = []
        while record is not None:
            ## LAYER, DATATYPE
            if record[0] in (0x0d, 0x0e):
                args.append(record[1][0])
            ## XY
            elif record[0] == 0x10:
                args.append(factor * record[1])
            ## ENDEL
            elif record[0] == 0x11:
                if create_element is not None:
                    cell.add(create_element(*args))
                    create_element = None
                args = []
            ## BOUNDARY
            elif record[0] == 0x08:
                create_element = self._create_polygon
            ## SNAME
            elif record[0] == 0x12:
                args.append(rename.get(record[1], record[1]))
            ## COLROW
            elif record[0] == 0x13:
                args.append(record[1])
            ## STRANS
            elif record[0] == 0x1a:
                args += [(record[1][0] & 0x8000) > 0, 1, 0]
            ## MAG
            elif record[0] == 0x1b:
                args[-2] = record[1][0]
            ## ANGLE
            elif record[0] == 0x1c:
                args[-1] = record[1][0]
            ## SREF
            elif record[0] == 0x0a:
                create_element = self._create_reference
            ## AREF
            elif record[0] == 0x0b:
                create_element = self._create_array
            ## STRNAME
            elif record[0] == 0x06:
                name = rename.get(record[1], record[1])
                cell = Cell(name, exclude_from_global=True)
                self.cell_list[name] = cell
            ## ENDSTR
            elif record[0] == 0x07:
                cell = None
            ## UNITS
            elif record[0] == 0x03:
                if unit is None:
                    factor = record[1][0]
                else:
                    factor = record[1][1] / unit
            ## ENDLIB
            elif record[0] == 0x04:
                for ref in self._incomplete:
                    if ref.ref_cell in self.cell_list:
                        ref.ref_cell = self.cell_list[ref.ref_cell]
                    else:
                        ref.ref_cell = Cell.cell_list.get(ref.ref_cell, ref.ref_cell)
            ## Not supported
            elif record[0] not in _unused_records:
                print "gdspy - WARNING: Record type {0} not supported by gds_import.".format(_record_name[record[0]])
            record = self._read_record(infile)

        if close:
            infile.close()

    def _read_record(self, stream):
        """
        Read a complete record from a GDSII stream file.

        Parameters
        ----------
        stream : file
            GDSII stream file to be imported.

        Returns
        -------
        out : 2-tuple
            Record type and data (as a numpy.array)
        """
        header = stream.read(4)
        if len(header) < 4:
            return None
        size, rec_type = struct.unpack('>HH', header)
        data_type = (rec_type & 0x00ff)
        rec_type = rec_type // 256
        data = None
        if size > 4:
            if data_type == 0x01:
                data = numpy.array(struct.unpack('>{0}H'.format((size - 4) // 2), stream.read(size - 4)), dtype='uint')
            elif data_type == 0x02:
                data = numpy.array(struct.unpack('>{0}h'.format((size - 4) // 2), stream.read(size - 4)), dtype=int)
            elif data_type == 0x03:
                data = numpy.array(struct.unpack('>{0}l'.format((size - 4) // 4), stream.read(size - 4)), dtype=int)
            elif data_type == 0x05:
                data = numpy.array([_eight_byte_real_to_float(stream.read(8)) for _ in range((size - 4) // 8)])
            else:
                data = stream.read(size - 4)
                if data[-1] == '\0':
                    data = data[:-1]
        return (rec_type, data)

    def _create_polygon(self, layer, datatype, xy):
        return Polygon(layer, xy[:-2].reshape((xy.size // 2 - 1, 2)), datatype)

    def _create_reference(self, name, *args):
        if len(args) > 1:
            ref = CellReference(name, args[3], args[2], args[1], args[0])
        else:
            ref = CellReference(name, args[0])
        if not isinstance(ref.ref_cell, Cell):
            self._incomplete.append(ref)
        return ref

    def _create_array(self, name, *args):
        xy = args[-1]
        columns = args[-2][0]
        rows = args[-2][1]
        if len(args) > 3:
            sa = -numpy.sin(args[2] * numpy.pi / 180.0)
            ca = numpy.cos(args[2] * numpy.pi / 180.0)
            x2 = (xy[2] - xy[0]) * ca - (xy[3] - xy[1]) * sa + xy[0]
            y3 = (xy[4] - xy[0]) * sa + (xy[5] - xy[1]) * ca + xy[1]
            if args[0]:
                y3 = 2 * xy[1] - y3
            ref = CellArray(name, columns, rows, ((x2 - xy[0]) / columns, (y3 - xy[1]) / rows), xy[0:2], args[2], args[1], args[0])
        else:
            ref = CellArray(name, columns, rows, ((xy[2] - xy[0]) / columns, (xy[5] - xy[1]) / rows), xy[0:2])
        if not isinstance(ref.ref_cell, Cell):
            self._incomplete.append(ref)
        return ref

    def extract(self, cell):
        """
        Extract a cell from the imported GDSII file and include it in the
        current scope, including referenced dependencies.

        Parameters
        ----------
        cell : ``Cell`` or string
            Cell or name of the cell to be extracted from the imported
            file. Referenced cells will be automatically extracted as well.

        Returns
        -------
        out : ``Cell``
            The extracted cell.
        """
        cells = Cell.cell_list.itervalues()
        cell = self.cell_list.get(cell, cell)
        for c in cell.get_dependencies():
            if c not in cells:
                self.extract(c)
        if cell not in cells:
            Cell.cell_list[cell.name] = cell
        return cell
    

def slice(layer, objects, position, axis, datatype=0):
    """
    Slice polygons and polygon sets at given positions along an axis.

    Parameters
    ----------
    layer : integer, list
        The GDSII layer numbers for the elements between each division.
        If the number of layers in the list is less than the number of
        divided regions, the list is repeated.
    objects : ``Polygon``, ``PolygonSet``, or list
        Operand of the slice operation.  If this is a list, each element
        must be a ``Polygon``, ``PolygonSet``, ``CellReference``,
        ``CellArray``, or an array-like[N][2] of vertices of a polygon.
    position : number or list of numbers
        Positions to perform the slicing operation along the specified
        axis.
    axis : 0 or 1
        Axis along which the polygon will be sliced.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and 255).

    Returns
    -------
    out : list[N] of PolygonSet
        Result of the slicing operation, with N = len(positions) + 1.
        Each PolygonSet comprises all polygons between 2 adjacent slicing
        positions, in crescent order.

    Examples
    --------
    >>> ring = gdspy.Round(1, (0, 0), 10, inner_radius = 5)
    >>> result = gdspy.slice(1, ring, [-7, 7], 0)
    >>> cell.add(result[1])
    """
    if (layer.__class__ != [].__class__):
        layer = [layer]
    if (objects.__class__ != [].__class__):
        objects = [objects]
    if (position.__class__ != [].__class__):
        position = [position]
    position.sort()
    result = [[] for i in range(len(position) + 1)]
    polygons = []
    for obj in objects:
        if isinstance(obj, Polygon):
            polygons.append(obj.points)
        elif isinstance(obj, PolygonSet):
            polygons += obj.polygons
        elif isinstance(obj, CellReference) or isinstance(obj, CellArray):
            polygons += obj.get_polygons()
        else:
            polygons.append(obj)
    for i, p in enumerate(position):
        nxt_polygons = []
        for pol in polygons:
            (pol1, pol2) = chop(pol, p, axis)
            result[i] += pol1
            nxt_polygons += pol2
        polygons = nxt_polygons
    result[-1] = polygons
    for i in range(len(result)):
        result[i] = PolygonSet(layer[i % len(layer)], result[i], datatype)
    return result


def boolean(layer, objects, operation, precision=1.0e-9, max_points=199, datatype=0):
    """
    Execute any boolean operation on polygons and polygon sets.

    Parameters
    ----------
    layer : integer
        The GDSII layer number for the resulting element.
    objects : array-like
        Operands of the boolean operation. Each element of this array
        must be a ``Polygon``, ``PolygonSet``, ``CellReference``,
        ``CellArray``, or an array-like[N][2] of vertices of a polygon.
    operation : function
        Function that accepts as input an array with ``len(objects)``
        elements. Each element is an integer representing the incidence
        of the corresponding element in ``objects``. The function must
        return a bool or integer (interpreted as bool).
    precision : number
        Precision of the grid where the operation if performed. All
        vertices in the geometry will be rounded to this precision.
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure
        they have at most ``max_points`` vertices. This is not a
        tessellating function, so this number should be as high as
        possible. For example, it should be set to 199 for polygons
        being drawn in GDSII files.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and 255).

    Returns
    -------
    out : PolygonSet
        Result of the boolean operation.

    Notes
    -----
    Since ``operation`` receives a list of integers as input, it can be
    somewhat more general than boolean operations only. See the examples
    below.

    Examples
    --------
    >>> circle = gdspy.Round(0, 10, (0, 0))
    >>> triangle = gdspy.Round(0, 12, (0, 0), number_of_points=3)
    >>> bad_poly = gdspy.L1Path(1, (0, 0), '+y', 2,
            [6, 4, 4, 8, 4, 5, 10], [-1, -1, -1, 1, 1, 1])
    >>> union = gdspy.boolean(1, [circle, triangle],
            lambda p: p[0] or p[1], precision=0.001)
    >>> intersection = gdspy.boolean(1, [circle, triangle],
            lambda p: p[0] and p[1], precision=0.001)
    >>> subtraction = gdspy.boolean(1, [circle, triangle],
            lambda p: p[0] - p[1], precision=0.001)
    >>> multi_xor = gdspy.boolean(1, [badPath], lambda p: p[0] % 2))
    """
    polygons = []     
    indices = []
    special_function = False
    for obj in objects:
        if isinstance(obj, Polygon):
            polygons.append(obj.points)
            indices.append(1)
        elif isinstance(obj, PolygonSet):
            special_function = True
            polygons += obj.polygons
            indices.append(len(obj.polygons))
        elif isinstance(obj, CellReference) or isinstance(obj, CellArray):
            special_function = True
            a = obj.get_polygons()
            polygons += a
            indices.append(len(a))
        else:
            polygons.append(obj)
            indices.append(1)
    if special_function:
        def function(p):
            a = numpy.zeros(len(indices))
            ip = 0
            ia = 0
            for ii in indices:
                a[ia] = numpy.sum(p[ip:ip+ii])
                ia += 1
                ip += ii
            return operation(a)
        return PolygonSet(layer, clip(polygons, function, precision, max_points), datatype)
    else:
        return PolygonSet(layer, clip(polygons, operation, precision, max_points), datatype)


def gds_print(outfile, cells=None, name='library', unit=1.0e-6, precision=1.0e-9):
    """
    Output a list of cells as a GDSII stream library.

    The dimensions actually written on the GDSII file will be the
    dimensions of the objects created times the ratio ``unit/precision``.
    For example, if a circle with radius 1.5 is created and we set
    ``unit=1.0e-6`` (1 um) and ``precision=1.0e-9`` (1 nm), the radius of
    the circle will be 1.5 um and the GDSII file will contain the
    dimension 1500 nm.
    
    Parameters
    ----------
    outfile : file or string
        The file (or path) where the GDSII stream will be written. It must
        be opened for writing operations in binary format.
    cells : array-like
        The list of cells or cell names to be included in the library. If
        ``None``, all cells listed in ``Cell.cell_list`` are used.
    name : string
        Name of the GDSII library (file).
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).

    Examples
    --------
    >>> gdspy.gds_print('out-file.gds', unit=1.0e-6, precision=1.0e-9)
    """
    if outfile.__class__ == ''.__class__:
        outfile = open(outfile, 'wb')
        close = True
    else:
        close = False
    if cells == None:
        cells = Cell.cell_list.itervalues()
    else:
        cells = [Cell.cell_list.get(c, c) for c in cells]
        i = 0
        while i < len(cells):
            for cell in cells[i].get_dependencies():
                if cell not in cells:
                    cells.append(cell)
            i += 1
    now = datetime.datetime.today()
    if len(name)%2 != 0:
        name = name + '\0'
    outfile.write(struct.pack('>19h', 6, 0x0002, 0x0258, 28, 0x0102, now.year, now.month, now.day, now.hour, now.minute, now.second, now.year, now.month, now.day, now.hour, now.minute, now.second, 4+len(name), 0x0206) + name + struct.pack('>2h', 20, 0x0305) + _eight_byte_real(precision / unit) + _eight_byte_real(precision))
    for cell in cells:
        outfile.write(cell.to_gds(unit / precision))
    outfile.write(struct.pack('>2h', 4, 0x0400))
