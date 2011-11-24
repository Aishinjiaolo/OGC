#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

class Point {
    public:
        double get_gx();
        double get_gy();
        
        void set_point(double gx, double gy);
        void get_center_of(Point point1, Point point2);

    private:
        double _gx;
        double _gy;
};

class Points {
    public:
        void append(Point point);
        void parse_1d_array_from(vector<double> point_1d_array);
        void sort_by_x();
        void sort_by_y();

        Point get_point(int index);
        int size();

    private:
        vector<Point> _points;
};

bool sort_by_gx_descending(Point a, Point b) {return a.get_gx() > b.get_gx();}
bool sort_by_gy_descending(Point a, Point b) {return a.get_gy() > b.get_gy();}
