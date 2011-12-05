#include "polygon.h"
#include <string>

typedef struct ktInterp_struct {
    int seg_index;
    int poly_index;
    Polygon *polygon;
} ktInterp;

void set_segment_function(void (*seg_function)(ktInterp *kt));

void loop_segment(Polygons *polygons);

class Profiler {
    public:
        void start(char *function_name);
        void end();

    private:
        struct timeval _start, _end;
        char *_function_name;
};

