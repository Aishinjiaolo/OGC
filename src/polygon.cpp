#include "polygon.h"

Segment *Polygon::get_segment(int index) {
    return _segments->get_segment(index);
}

void Polygon::set_polygon(Segments *segments) {
    _segments = segments;
    _segment_numbers = segments->size();
}

Point Polygon::get_center() {
    double x_total = 0, y_total = 0;
    for (int i = 0; i < _segment_numbers; i++) {
        x_total += _segments->get_segment(i)->get_tail()->get_gx();
        y_total += _segments->get_segment(i)->get_tail()->get_gy();
    }

    Point center;
    center.set_point(x_total/_segment_numbers, y_total/_segment_numbers);
    return center;
}

int Polygon::get_segment_number() {
    return _segment_numbers;
}
