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

    Point a, b;
    a.set_point(0, 0);
    b.set_point(0, 0);
    
    Segment c;
    c.set_segment(a, b);
    
    cout << "head = " << c.get_head().get_gx() << ", "
        << c.get_head().get_gy() << endl;
    cout << "tail = " << c.get_tail().get_gx() << ", "
        << c.get_tail().get_gy() << endl;
    
    cout << "length = " << c.get_length() << endl;
    cout << "angle = " << c.get_angle() << endl;
    cout << "dir = " << c.get_dir() << endl;
    
    profiler.end();
    
    return 0;
}
