#include "point.h"
#include "gtest/gtest.h"

TEST(PointTest, SetPoint) {
    srand(time(NULL));
    for (int i = 0 ; i < 100; i++) {
        Point test;
        double x = rand() % 100;
        double y = rand() % 100;
        test.set_point(x, y);

        EXPECT_EQ(x, test.get_gx());
        EXPECT_EQ(y, test.get_gy());
    }
}


