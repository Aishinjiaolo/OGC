#include "point.h"

void Points::parse_1d_array_from(vector<double> point_1d_array) {
    for (unsigned int i = 0; i < 0.5 * point_1d_array.size(); i++) {
        Point point;
        unsigned int index = i * 2;
        point.set_point(point_1d_array[index], point_1d_array[index+1]);
        Points::append(point);
    }
}

int Points::size() {
    return _array.size();
}

Point Points::get_point(int index) {
    return _array[index];
}

void Points::append(Point point) {
    _array.push_back(point);
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

int main() {
    double main_1d_array[] = {5, 5, 0, 5, 10, 15, 20, 25, 25, 5, 10, 0};
    vector<double> main_1d_points(main_1d_array, main_1d_array + 12);

    Points main_points;
    main_points.parse_1d_array_from(main_1d_points);

    for (int i = 0; i < main_points.size(); i++) {
        cout << main_points.get_point(i).get_gx() << "\t"
            << main_points.get_point(i).get_gy() << endl;
    }

    return 0;
}
