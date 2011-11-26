#include "segment.h"

class Polygon {
    public:
        void set_polygon(Segments segments);
        
        Point  get_center();
        double get_area();
        int    size();

    private:
        double _area;
        int    _size;
};

class Polygons {
    public:
        void append(Polygon polygon);
        void sort();

        Polygon get_polygon(int index);
        int size();

    private:
        // each polygon is uncertain in size?
        vector<Polygon> _polygons;
};
