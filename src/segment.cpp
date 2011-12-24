#include "segment.h"

double Segment::get_segment_property(int kt_spt_index) {
    return _property[kt_spt_index];
}

void Segment::free() {
    delete _head;
    delete _tail;
}

void Segment::copy(Segment context) {
    double head_x = context.get_head()->get_gx();
    double head_y = context.get_head()->get_gy();
    double tail_x = context.get_tail()->get_gx();
    double tail_y = context.get_tail()->get_gy();

    Point *new_head = new Point;
    Point *new_tail = new Point;

    new_head->set_point(head_x, head_y);
    new_tail->set_point(tail_x, tail_y);

    this->set_segment(new_head, new_tail);
}

void Segment::dump() {
    printf("head: ");
    _head->dump();
    printf("tail: ");
    _tail->dump();
}

void Segments::dump() {
    for (unsigned int i = 0; i < _segments.size(); i++) {
        printf("segment(%d):\n", i);
        _segments[i]->dump();
    }
}

unsigned int Segments::size() {
    return _segments.size();
}

void Segments::append(Segment *segment) {
    _segments.push_back(segment);
}

Segment *Segments::get_segment(int index) {
    return _segments[index];
}

Point *Segment::get_head() {
    return _head; 
}

Point *Segment::get_tail() {
    return _tail;
}

int Segment::get_dir() {
    return _property[KT_SPT_DIRECTION];
}

double Segment::get_length() {
    return _property[KT_SPT_LENGTH];
}

double Segment::get_angle() {
    return _property[KT_SPT_ANGLE];
}

void Segment::set_segment(Point *head, Point *tail) {
    _head = head;
    _tail = tail;
    
    double bias_x = tail->get_gx() - head->get_gx();
    double bias_y = tail->get_gy() - head->get_gy();
    bool is_same_point =
        abs(bias_x) <= EPS && abs(bias_y) <= EPS ? true : false;

    double angle, length;
    length = is_same_point ? 0 : sqrt(pow(bias_x, 2) + pow(bias_y, 2));
    angle  = is_same_point ? 0 : atan(bias_y / bias_x) * 180 / M_PI;

    if (bias_x < 0) angle += HALF_RADIUS;
    if (angle  < 0) angle += RADIUS;
    
    int angle_snap = static_cast<int>(angle);
    int dir;
    switch(angle_snap) {
    case 0:
        dir = 0;
        break;
    case 45:
        dir = 1;
        break;
    case 90:
        dir = 2;
        break;
    case 135:
        dir = 3;
        break;
    case 180:
        dir = 4;
        break;
    case 225:
        dir = 5;
        break;
    case 270:
        dir = 6;
        break;
    case 315:
        dir = 7;
        break;
    default:
        dir = 0;
    }

    _property[KT_SPT_GX]         = tail->get_gx();
    _property[KT_SPT_GY]         = tail->get_gy();
    _property[KT_SPT_LENGTH]     = length;
    _property[KT_SPT_ANGLE]      = angle;
    _property[KT_SPT_DIRECTION]  = dir;
}
