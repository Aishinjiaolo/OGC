#include "api.h"

static void (*segment_function)(ktInterp *kt) = NULL;

void set_segment_function(
        void (*seg_function)(ktInterp *kt)) {
    segment_function = seg_function;
}

void loop_segment(Polygon *polygon) {
    if (segment_function == NULL) {
        printf("no segment function set!!\n");
        return;
    }

    for (unsigned int i = 0; i < polygon->get_segment_number(); i++) {
        ktInterp kt;
        kt.seg_index = i;
        kt.polygon   = polygon;
        segment_function(&kt);
    }
    
    segment_function = NULL;
}

void Profiler::end() {
    gettimeofday(&_end, NULL);
    cout << "Time elapsed " <<
        _end.tv_usec - _start.tv_usec << " usec" << endl;
}

void Profiler::start() {
    gettimeofday(&_start, NULL);
}
