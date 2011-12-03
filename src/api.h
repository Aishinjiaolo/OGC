#include "polygon.h"

typedef struct ktInterp_struct {
    int seg_index;
    int poly_index;
    Polygon *polygon;
} ktInterp;

static void (*segment_function)(ktInterp *kt) = NULL;

void set_segment_function(
        void (*seg_function)(ktInterp *kt)) {
    segment_function = seg_function;
}

void loop_segment(Polygon *polygon);

class Profiler {
    public:
        void start();
        void end();

    private:
        struct timeval _start, _end;
};

