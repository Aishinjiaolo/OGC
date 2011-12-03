#include "api.h"

void loop_segment(Polygon *polygon) {
    if (segment_function == NULL) return;

    for (int i = 0; i < polygon->get_segment_number(); i++) {
        ktInterp kt;
        kt.seg_index = i;
        kt.polygon   = polygon;
        segment_function(&kt);
    }
}

void Profiler::end() {
    gettimeofday(&_end, NULL);
    cout << "Time elapsed " <<
        _end.tv_usec - _start.tv_usec << " usec" << endl;
}

void Profiler::start() {
    gettimeofday(&_start, NULL);
}
