#include "point.h"

class Segment {
    public:
        //Segment() {printf("Segment()\n");};
        //~Segment() {printf("~Segment()\n");};

        void set_segment(Point *head, Point *tail);
        void set_segment_property(int kt_spt_index, double value);
        void set_segment_property_int(int kt_spt_index, int value);
        void dump();
        void copy(Segment context);
        void free();

        Point  *get_head();
        Point  *get_tail();

        int get_dir();
        int get_segment_property_int(int kt_spt_index);

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
        //Segments() {printf("Segments()\n");};
        //~Segments() {printf("~Segments()\n");};

        void append(Segment *segment);
        //void sort_by_x();
        //void sort_by_y();
        void dump();

        Segment *get_segment(int index);

        unsigned int size();

    private:
        vector<Segment*> _segments;
};
