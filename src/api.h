#include "polygon.h"
#include <string>

typedef struct ktInterp_struct {
    int seg_index;
    int poly_index;
    Polygon *polygon;
} ktInterp;

double ktGetSegmentProperty(
        ktInterp *kt, int figure_type, int offset, int kt_spt_index);

void ktSetSegmentFunction(void (*seg_function)(ktInterp *kt));
void ktLoopSegment(Polygons *polygons);

class Profiler {
    public:
        void start(char *function_name);
        void end();

    private:
        struct timeval _start, _end;
        char *_function_name;
};

