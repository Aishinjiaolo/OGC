#include "segment.h"

class Polygon {
    public:
        //Polygon() {printf("Polygon()\n");};
        //~Polygon() {printf("~Polygon()\n");};

        void set_polygon(Segments *segments);
        void dump();
        void copy(Polygon context);
        void free();

        //double get_area();

        unsigned int get_segment_number();

        Point  get_center();

        Segment *get_segment(int index);

    private:
        Segments *_segments;

        double _area;

        unsigned int _segment_numbers;
};

class Polygons {
    public:
        void append(Polygon *polygon);
        //void sort_by_x();
        //void sort_by_y();
        void dump();

        Polygon *get_polygon(int index);

        unsigned int size();

    private:
        vector<Polygon*> _polygons;
};


