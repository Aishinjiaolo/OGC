#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/time.h>

using namespace std;

class Profiler {
    public:
        void start();
        void end();

    private:
        struct timeval _start, _end;
};

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
        void parse_1d_points_from(vector<double> points_1d);
        void sort_by_x();
        void sort_by_y();
        void create_center_from(Points potins);
        
        Point get_point(int index);
        int size();

    private:
        vector<Point> _points;
};


