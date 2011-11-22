#include <iostream>
#include "point.h"

using namespace std;

void PointArray::append(Point point) {
}

void Point::get_center_point(Point point1 , Point point2) {
    _gx = 0.5 * (point1._gx + point2._gx);
    _gy = 0.5 * (point1._gy + point2._gy);
}

void Point::set_point(double gx, double gy) {
    _gx = gx;
    _gy = gy;
}

double Point::get_point_x() {
    return _gx;
}

double Point::get_point_y() {
    return _gy;
}

int main() {
    Point point1, point2, point3;
    point1.set_point(10, 20);
    point2.set_point(10, 40);

    point3.get_center_point(point1, point2);

    cout << "center x, y = " << point3.get_point_x() << ", "
        << point3.get_point_y() << endl;
    
    return 0;
}
