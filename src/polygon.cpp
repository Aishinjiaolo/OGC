#include "polygon.h"

void Polygons::dump() {
    printf("\n In these polygons dump:\n");
    for (unsigned int i = 0; i < _polygons.size(); i++) {
        printf("\npolygon(%d):\n", i);
        _polygons[i]->dump();
    }
}

void Polygon::dump() {
    printf("\n In this polygon dump:\n");
    for (unsigned int i = 0; i < _segment_numbers; i++) {
        printf("segment(%d):\n", i);
        _segments->get_segment(i)->dump();
    }
}

Polygon *Polygons::get_polygon(int index) {
    return _polygons[index];
}

unsigned int Polygons::size() {
    return _polygons.size();
}

void Polygons::append(Polygon *polygon) {
    _polygons.push_back(polygon);
}

Segment *Polygon::get_segment(int index) {
    return _segments->get_segment(index);
}

void Polygon::set_polygon(Segments *segments) {
    _segments = segments;
    _segment_numbers = segments->size();
}

Point Polygon::get_center() {
    double x_total = 0, y_total = 0;
    for (unsigned int i = 0; i < _segment_numbers; i++) {
        x_total += _segments->get_segment(i)->get_tail()->get_gx();
        y_total += _segments->get_segment(i)->get_tail()->get_gy();
    }

    Point center;
    center.set_point(x_total/_segment_numbers, y_total/_segment_numbers);
    return center;
}

unsigned int Polygon::get_segment_number() {
    return _segment_numbers;
}
