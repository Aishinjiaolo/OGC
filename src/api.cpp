#include "api.h"

double ktGetSegmentProperty(
        ktInterp *kt, int figure_type, int offset, int kt_spt_index) {
    return kt->polygon->get_segment(kt->seg_index+offset)->
        get_segment_property(kt_spt_index);
}

static void (*segment_function)(ktInterp *kt) = NULL;

void ktSetSegmentFunction(
        void (*seg_function)(ktInterp *kt)) {
    segment_function = seg_function;
}

void ktLoopSegment(Polygons *polygons) {
    if (segment_function == NULL) {
        printf("no segment function set!!\n");
        return;
    }

    ktInterp kt;
    unsigned int total_polygon_number = polygons->size();
    for (unsigned int i = 0; i < total_polygon_number; i++) {
        kt.polygon = polygons->get_polygon(i);
        unsigned int total_segment_number = kt.polygon->get_segment_number();
        for (unsigned int j = 0; j < total_segment_number; j++) {
            kt.seg_index = j;
            segment_function(&kt);
        }
    }

    segment_function = NULL;
}

void Profiler::end() {
    gettimeofday(&_end, NULL);
    printf("Function '%s' 1 call, time elapsed %d usec.\n",
            _function_name, _end.tv_usec - _start.tv_usec);
}

void Profiler::start(char *function_name) {
    _function_name = function_name;
    gettimeofday(&_start, NULL);
}
