#include <iostream>
#include <string>
#include <vector>

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
        Point get_point(int index);
        int size();

    private:
        vector<Point> _array;
};
