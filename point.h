#include <string>
using namespace std;

class Point {
    public:
        void set_point(double gx, double gy);
        
        double get_point_x();
        double get_point_y();
        
        void get_center_point(Point point1, Point point2);

    private:
        double _gx;
        double _gy;
};

class PointArray {
    public:
        void append(Point point);

    private:
        Point *_array;
};
