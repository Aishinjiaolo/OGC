#include "segment.h"

class Polygon {
    public:
        void set_polygon(Segments *segments);
        
        Point  get_center();
        double get_area();
        int    get_segment_number();

    private:
        Segments *_segments;
        double _area;
        int    _segment_numbers;
};

class Polygons {
    public:
        void append(Polygon *polygon);
        void sort_by_x();
        void sort_by_y();

        Polygon get_polygon(int index);
        int get_polygon_number();

    private:
        vector<Polygon*> _polygons;
};
