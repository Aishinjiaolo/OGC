#include "segment.h"

void Segment::dump() {
    printf("head: ");
    _head->dump();
    printf("tail: ");
    _tail->dump();
}

void Segments::dump() {
    for (int i = 0; i < _segments.size(); i++) {
        printf("segment(%d):\n", i);
        _segments[i]->dump();
    }
}

int Segments::size() {
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
    return _dir;
}

double Segment::get_length() {
    return _length;
}

double Segment::get_angle() {
    return _angle;
}

void Segment::set_segment(Point *head, Point *tail) {
    _head = head;
    _tail = tail;
    
    double bias_x = tail->get_gx() - head->get_gx();
    double bias_y = tail->get_gy() - head->get_gy();
    bool is_same_point =
        abs(bias_x) <= EPS && abs(bias_y) <= EPS ? true : false;

    _length = is_same_point ? 0 : sqrt(pow(bias_x, 2) + pow(bias_y, 2));
    _angle  = is_same_point ? 0 : atan(bias_y / bias_x) * 180 / M_PI;
    if (bias_x < 0) _angle += 180;
    if (_angle < 0) _angle += 360;

    int angle = static_cast<int>(_angle);
    switch(angle) {
    case 0:
        _dir = 0;
        break;
    case 45:
        _dir = 1;
        break;
    case 90:
        _dir = 2;
        break;
    case 135:
        _dir = 3;
        break;
    case 180:
        _dir = 4;
        break;
    case 225:
        _dir = 5;
        break;
    case 270:
        _dir = 6;
        break;
    case 315:
        _dir = 7;
        break;
    default:
        _dir = 0;
    }
}
