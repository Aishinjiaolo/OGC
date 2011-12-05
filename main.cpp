#include "api.h"

// callback segment function
void print_vertex(ktInterp *kt) {
    int i = kt->seg_index;
    printf("(%d)\n", i);
    kt->polygon->get_segment(i)->dump();
}

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

    Point a, b, c, d, e, f, g, h;
    a.set_point( 0,  0);
    b.set_point( 5,  0);
    c.set_point( 5,  5);
    d.set_point( 0,  5);
    e.set_point(-5,  5);
    f.set_point(-5,  0);
    g.set_point(-5, -5);
    h.set_point( 0, -5);
    
    Segment ab, bc, cd, de, ef, fg, gh, ha;
    ab.set_segment(&a, &b);
    bc.set_segment(&b, &c);
    cd.set_segment(&c, &d);
    de.set_segment(&d, &e);
    ef.set_segment(&e, &f);
    fg.set_segment(&f, &g);
    gh.set_segment(&g, &h);
    ha.set_segment(&h, &a);

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

    Polygon polygon1;
    polygon1.set_polygon(&all);

    set_segment_function(grow_45_degree);
    loop_segment(&polygon1);

    all.dump();

    // test polygon function
    Polygon polygon2;
    polygon2.set_polygon(&all);
    cout << polygon2.get_segment_number()
        << " segments in this polygon" << endl;
    cout << "center point of this polygon is: "
        << polygon2.get_center().get_gx() << ", "
        << polygon2.get_center().get_gy() << endl;

    cout << "try segment loop on a polygon :" << endl;

    set_segment_function(print_vertex);
    loop_segment(&polygon2);

    // try loop segment without function set
    loop_segment(&polygon2);
    profiler.end();

    printf("\ntry to dump polygon\n");
    polygon1.dump();
    polygon2.dump();

    Polygons polygons;
    polygons.append(&polygon1);
    polygons.append(&polygon2);
    printf("\ntry to dump polygons\n");
    polygons.dump();

    return 0;
}
