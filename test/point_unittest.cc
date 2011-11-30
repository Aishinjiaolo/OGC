#include "point.h"
#include "gtest/gtest.h"

TEST(PointTest, SetPoint) {
    Point test;
    double x =  5;
    double y = 10;
    test.set_point(x, y);

    EXPECT_EQ(x, test.get_gx());
    EXPECT_EQ(y, test.get_gy());
}

