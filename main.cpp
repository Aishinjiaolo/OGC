#include "polygon.h"

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
    
    for (int i = 0; i < main_points.size(); i++) {
        cout << main_points.get_point(i).get_gx() << "\t"
            << main_points.get_point(i).get_gy() << endl;
    }

    Points center;
    center.create_center_from(main_points);

    for (int i = 0; i < center.size(); i++) {
        cout << center.get_point(i).get_gx() << "\t"
            << center.get_point(i).get_gy() << endl;
    }

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

    for (int i = 0; i < all.size(); i++) {
        cout << "segment " << i << ": "
            << all.get_segment(i)->get_head()->get_gx() << ", "
            << all.get_segment(i)->get_head()->get_gy() << endl;
        cout << "           "
            << all.get_segment(i)->get_tail()->get_gx() << ", "
            << all.get_segment(i)->get_tail()->get_gy() << endl;
    }

    for (int i = 0; i < all.size(); i++) {
        double tail_x = all.get_segment(i)->get_tail()->get_gx();
        double tail_y = all.get_segment(i)->get_tail()->get_gy();
        if (tail_x == tail_y) {
            double new_tail_x = tail_x >= 0 ? tail_x + 1 : tail_x - 1;
            double new_tail_y = new_tail_x;

            Point new_tail;
            new_tail.set_point(new_tail_x, new_tail_y);

            // doest not work as expected
            all.get_segment(i)->set_segment(all.get_segment(i)->get_head(),
                    &new_tail);
        }
    }

    for (int i = 0; i < all.size(); i++) {
        cout << "segment " << i << ": "
            << all.get_segment(i)->get_head()->get_gx() << ", "
            << all.get_segment(i)->get_head()->get_gy() << endl;
        cout << "           "
            << all.get_segment(i)->get_tail()->get_gx() << ", "
            << all.get_segment(i)->get_tail()->get_gy() << endl;
    }

    profiler.end();

    return 0;
}
