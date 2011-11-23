#include "point.h"

int PointArray::size() {
    return _array.size();
}

Point PointArray::get_point(int index) {
    return _array[index];
}

void PointArray::append(Point point) {
    _array.push_back(point);
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

    PointArray array;
    array.append(point1);
    array.append(point2);
    array.append(point3);

    cout << "array = " << endl;
    int i = 0;
    for (i = 0; i < array.size(); i++) {
        cout << array.get_point(i).get_point_x() << "\t"
            << array.get_point(i).get_point_y() << endl;
    }

    return 0;
}
