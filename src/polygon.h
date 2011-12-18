#include "segment.h"

class Polygon {
    public:
        void set_polygon(Segments *segments);
        
        Point  get_center();
        double get_area();
        unsigned int get_segment_number();

        Segment *get_segment(int index);

        void dump();
        void copy(Polygon context);

    private:
        Segments *_segments;
        double _area;
        unsigned int _segment_numbers;
};

class Polygons {
    public:
        void append(Polygon *polygon);
        void sort_by_x();
        void sort_by_y();

        Polygon *get_polygon(int index);
        unsigned int size();

        void dump();

    private:
        vector<Polygon*> _polygons;
};


