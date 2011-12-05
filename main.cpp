#include "api.h"

// callback segment function
void grow_45_degree(ktInterp *kt) {
    int i = kt->seg_index;
    double tail_x = kt->polygon->get_segment(i)->get_tail()->get_gx();
    double tail_y = kt->polygon->get_segment(i)->get_tail()->get_gy();

    if (tail_x == tail_y) {
        double new_tail_x = tail_x >= 0 ? tail_x + 1 : tail_x - 1;
        double new_tail_y = new_tail_x;
        kt->polygon->get_segment(i)->get_tail()
            ->set_point(new_tail_x, new_tail_y);
    }
}

int main() {
    // profiler
    Profiler profiler;
    profiler.start();

    // random generate 1d vector
    // assign maximum size =>
    const int max = 100;
    vector<double> main_1d_points(max);

    srand(time(NULL));
    for (int i = 0; i < max; i++) {
        main_1d_points[i] = rand() % 100;
    };

    // transform to Points, then sorting
    Points main_points;
    main_points.parse_1d_points_from(main_1d_points);
    main_points.sort_by_y();

    main_points.dump();

    Points center;
    center.create_center_from(main_points);

    center.dump();

    // create a set of points
    Point a, b, c, d, e, f, g, h;
    a.set_point( 0,  0);
    b.set_point( 5,  0);
    c.set_point( 5,  5);
    d.set_point( 0,  5);
    e.set_point(-5,  5);
    f.set_point(-5,  0);
    g.set_point(-5, -5);
    h.set_point( 0, -5);

    // form each segment
    Segment ab, bc, cd, de, ef, fg, gh, ha;
    ab.set_segment(&a, &b);
    bc.set_segment(&b, &c);
    cd.set_segment(&c, &d);
    de.set_segment(&d, &e);
    ef.set_segment(&e, &f);
    fg.set_segment(&f, &g);
    gh.set_segment(&g, &h);
    ha.set_segment(&h, &a);

    // form a set of segments
    Segments all;
    all.append(&ab);
    all.append(&bc);
    all.append(&cd);
    all.append(&de);
    all.append(&ef);
    all.append(&fg);
    all.append(&gh);
    all.append(&ha);

    all.dump();

    // set the set of segments to a polygon
    // and test polygon function
    Polygon polygon1;
    polygon1.set_polygon(&all);
    polygon1.dump();
    printf("(%d) segments in this polygon\n", polygon1.get_segment_number());
    printf("center point of this polygon is (%f, %f)\n",
            polygon1.get_center().get_gx(),
            polygon1.get_center().get_gy());

    // set the set of segments to another polygon
    // yet, they reference the same segments set
    Polygon polygon2;
    polygon2.set_polygon(&all);
    polygon2.dump();

    // append the polygons into a set of polygons
    Polygons polygons;
    polygons.append(&polygon1);
    polygons.append(&polygon2);

    // set segment function and loop all polygons set
    set_segment_function(grow_45_degree);
    loop_segment(&polygons);

    polygons.dump();

    // try loop segment without function set
    loop_segment(&polygons);

    profiler.end();

    return 0;
}
