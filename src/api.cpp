#include "api.h"

static void (*segment_function)(ktInterp *kt) = NULL;

void set_segment_function(
        void (*seg_function)(ktInterp *kt)) {
    segment_function = seg_function;
}

void loop_segment(Polygons *polygons) {
    if (segment_function == NULL) {
        printf("no segment function set!!\n");
        return;
    }

    for (unsigned int i = 0; i < polygons->size(); i++) {
        Polygon *polygon = polygons->get_polygon(i);
        for (unsigned int j = 0; j < polygon->get_segment_number(); j++) {
            ktInterp kt;
            kt.seg_index = j;
            kt.polygon   = polygon;
            segment_function(&kt);
        }
    }
    
    segment_function = NULL;
}

void Profiler::end() {
    gettimeofday(&_end, NULL);
    printf("Time elapsed %d usec.\n", _end.tv_usec - _start.tv_usec);
}

void Profiler::start() {
    gettimeofday(&_start, NULL);
}
