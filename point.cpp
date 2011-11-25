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

void Points::parse_1d_array_from(vector<double> point_1d_array) {
    for (unsigned int i = 0; i < 0.5 * point_1d_array.size(); i++) {
        Point point;
        unsigned int index = i * 2;
        point.set_point(point_1d_array[index], point_1d_array[index+1]);
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

int main() {
    struct timeval start, end;
    gettimeofday(&start, NULL);
   
    const int max = 100;
    double main_1d_array[max];
    srand(time(NULL));

    for (int i = 0; i < max; i++) {
        main_1d_array[i] = rand() % 100;
    };

    vector<double> main_1d_points(main_1d_array, main_1d_array + max);

    Points main_points;
    main_points.parse_1d_array_from(main_1d_points);
    main_points.sort_by_y();
    
    for (int i = 0; i < main_points.size(); i++) {
        cout << main_points.get_point(i).get_gx() << "\t"
            << main_points.get_point(i).get_gy() << endl;
    }
    
    gettimeofday(&end, NULL);
    cout << "Time elapsed " <<
        end.tv_usec - start.tv_usec << " usec" << endl;
    
    return 0;
}
