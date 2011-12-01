#include "segment.h"
#include "gtest/gtest.h"

TEST(SegmentTest, SetSegment) {
    srand(time(NULL));
    int size = 1000000;
    for (int i = 0; i < size; i++) {
        double digit = (rand() % 100) * GRID;
        double x1 = (rand() % 100) * digit;
        double y1 = (rand() % 100) * digit;
        double x2 = (rand() % 100) * digit;
        double y2 = (rand() % 100) * digit;

        Point head, tail;
        head.set_point(x1, y1);
        tail.set_point(x2, y2);

        Segment segment;
        segment.set_segment(&head, &tail);

        EXPECT_EQ(x1, segment.get_head()->get_gx());
        EXPECT_EQ(y1, segment.get_head()->get_gy());
        EXPECT_EQ(x2, segment.get_tail()->get_gx());
        EXPECT_EQ(y2, segment.get_tail()->get_gy());
      
        EXPECT_TRUE(segment.get_dir() >= 0);
        EXPECT_TRUE(segment.get_dir() <= 7);

        EXPECT_TRUE(segment.get_angle() >=   0);
        EXPECT_TRUE(segment.get_angle() <  360);
    }
}

TEST(SegmentTest, GetProperty) {
    srand(time(NULL));
    int size = 1000000;
    for (int i = 0; i < size; i++) {
        double digit = (rand() % 100) * GRID;
        double length = (rand() % 100000) * digit;
        double x1 = (rand() % 100) * digit;
        double y1 = (rand() % 100) * digit;
        double x2 = x1;
        double y2 = y1 + length;

        Point head, tail;
        head.set_point(x1, y1);
        tail.set_point(x2, y2);

        Segment segment;
        segment.set_segment(&head, &tail);

        EXPECT_EQ(length, segment.get_length());
        EXPECT_TRUE(segment.get_length() >= 0);
        if (length > EPS) EXPECT_EQ(2, segment.get_dir());
    }

    for (int i = 0; i < size; i++) {
        int digit = (rand() % 100) * GRID;
        double length = (rand() % 100000) * digit;
        double x1 = (rand() % 100) * digit;
        double y1 = (rand() % 100) * digit;
        double x2 = x1 + length;
        double y2 = y1;

        Point head, tail;
        head.set_point(x1, y1);
        tail.set_point(x2, y2);

        Segment segment;
        segment.set_segment(&head, &tail);

        EXPECT_EQ(length, segment.get_length());
        EXPECT_TRUE(segment.get_length() >= 0);
        if (length > EPS) EXPECT_EQ(0, segment.get_dir());
    }
}
