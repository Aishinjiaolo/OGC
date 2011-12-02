#include "polygon.h"
#include "gtest/gtest.h"

TEST(PolygonTest, SetPolygon) {
    srand(time(NULL));
    int size = 1000000;
    for (int i = 0; i < size; i++) {
        double digit = (rand() % 100) * GRID;
        double x1 = (rand() % 100) * digit;
        double y1 = (rand() % 100) * digit;
        double x2 = (rand() % 100) * digit;
        double y2 = (rand() % 100) * digit;
        double x3 = (rand() % 100) * digit;
        double y3 = (rand() % 100) * digit;
        double x4 = (rand() % 100) * digit;
        double y4 = (rand() % 100) * digit;

        Point v1, v2, v3, v4;
        v1.set_point(x1, y1);
        v2.set_point(x2, y2);
        v3.set_point(x3, y3);
        v4.set_point(x4, y4);

        Segment s1, s2, s3, s4;
        s1.set_segment(&v1, &v2);
        s2.set_segment(&v2, &v3);
        s3.set_segment(&v3, &v4);
        s4.set_segment(&v4, &v1);

        Segments ss;
        ss.append(&s1);
        ss.append(&s2);
        ss.append(&s3);
        ss.append(&s4);

        Polygon p;
        p.set_polygon(&ss);

        EXPECT_EQ(4, p.get_segment_number());
    }
}
