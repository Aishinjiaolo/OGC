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
    profiler.start((char*)"main");

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

    // try segment copy
    Segment copy_test;
    copy_test.copy(ab);
    printf("copy segment:\n");
    copy_test.dump();
    copy_test.free();

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
    printf("\nthis is polygon2\n");
    Polygon polygon2;
    polygon2.set_polygon(&all);
    polygon2.dump();

    // try polygon copy
    printf("\nthis is polygon3\n");
    Polygon polygon3;
    polygon3.copy(polygon2);
    polygon3.dump();

    // append the polygons into a set of polygons
    Polygons polygons;
    polygons.append(&polygon1);
    polygons.append(&polygon3);

    // set segment function and loop all polygons set
    ktSetSegmentFunction(grow_45_degree);
    ktLoopSegment(&polygons);

    polygons.dump();

    // try loop segment without function set
    ktLoopSegment(&polygons);

    polygon3.free();
    profiler.end();

    return 0;
}
