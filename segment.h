#include <cmath>
#include "point.h"

class Segment {
    public:
        void set_segment(Point head, Point tail);
        
        Point  get_head();
        Point  get_tail();
        int    get_dir();
        double get_length();
        double get_angle();

    private:
        Point  _head;
        Point  _tail;
        int    _dir;
        double _length;
        double _angle;
};

class Segments {
    public:
        void append(Segment segment);
        void sort();
        
        Segment get_segment(int index);
        int size();

    private:
        vector<Segment> _segments;
};
