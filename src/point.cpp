#include "point.h"

bool sort_by_gx_descending(Point a, Point b) {
    if (a.get_gx() > b.get_gx()) return true;
    if (a.get_gx() == b.get_gx() && a.get_gy() > b.get_gy()) return true;
    return false;
}

bool sort_by_gy_descending(Point a, Point b) {
    if (a.get_gy() > b.get_gy()) return true;
    if (a.get_gy() == b.get_gy() && a.get_gx() > b.get_gx()) return true;
    return false;
}

void Points::sort_by_y() {
    sort(_points.begin(), _points.end(), sort_by_gy_descending);
}

void Points::sort_by_x() {
    sort(_points.begin(), _points.end(), sort_by_gx_descending);
}

void Points::create_center_from(Points points) {
    int size = points.size();
    for (int i = 0; i < size; i++) {
        Point current = points.get_point(i);

        for (int j = i + 1; j < size; j++) {
            Point context = points.get_point(j);

            Point center;
            center.get_center_of(current, context);
            
            Points::append(center);
        }
    }
}

void Points::parse_1d_points_from(vector<double> points_1d) {
    for (unsigned int i = 0; i < 0.5 * points_1d.size(); i++) {
        Point point;
        unsigned int index = i * 2;
        point.set_point(points_1d[index], points_1d[index+1]);
        Points::append(point);
    }
}

int Points::size() {
    return _points.size();
}

Point Points::get_point(int index) {
    return _points[index];
}

void Points::append(Point point) {
    _points.push_back(point);
}

void Point::get_center_of(Point point1 , Point point2) {
    _gx = 0.5 * (point1._gx + point2._gx);
    _gy = 0.5 * (point1._gy + point2._gy);
}

void Point::set_point(double gx, double gy) {
    _gx = gx;
    _gy = gy;
}

double Point::get_gx() {
    return _gx;
}

double Point::get_gy() {
    return _gy;
}

void Points::dump() {
    for (unsigned int i = 0; i < _points.size(); i++) {
        printf("(%d): ", i);
        _points[i].dump();
    }
}

void Point::dump() {
    printf("point(x, y) = %f, %f\n", _gx, _gy);
}
