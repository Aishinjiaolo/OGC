#include "point.h"
#include "gtest/gtest.h"

TEST(PointTest, SetPoint) {
    srand(time(NULL));
    int max = 1000000;
    for (int i = 0 ; i < max; i++) {
        Point test;
        double x = rand() % 100;
        double y = rand() % 100;
        test.set_point(x, y);

        EXPECT_EQ(x, test.get_gx());
        EXPECT_EQ(y, test.get_gy());
    }
}

TEST(PointsTest, EmptySize) {
    Points points;
    EXPECT_EQ(0, points.size());
}

TEST(PointsTest, PointSize) {
    srand(time(NULL));
    int size = rand() % 1000000;
    Points points;
    for (int i = 0 ; i < size; i++) {
        Point test;
        points.append(test);
    }

    EXPECT_EQ(size, points.size());
}


