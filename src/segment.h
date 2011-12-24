#include "point.h"

class Segment {
    public:
        void set_segment(Point *head, Point *tail);
        void dump();
        void copy(Segment context);
        void free();

        Point  *get_head();
        Point  *get_tail();

        int    get_dir();

        double get_length();
        double get_angle();
        double get_segment_property(int kt_spt_index);

    private:
        Point  *_head;
        Point  *_tail;

        int    _dir;

        double _length;
        double _angle;

        double _property[KT_MAX_SPT_SIZE];
};

class Segments {
    public:
        void append(Segment *segment);
        //void sort_by_x();
        //void sort_by_y();
        void dump();

        Segment *get_segment(int index);

        unsigned int size();

    private:
        vector<Segment*> _segments;
};
