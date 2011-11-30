#include "segment.h"
#include "gtest/gtest.h"

TEST(SegmentTest, SetSegment) {
    srand(time(NULL));
    int size = 1000000;
    for (int i = 0; i < size; i++) {
        double x1 = rand() % 100;
        double y1 = rand() % 100;
        double x2 = rand() % 100;
        double y2 = rand() % 100;

        Point head, tail;
        head.set_point(x1, y1);
        tail.set_point(x2, y2);

        Segment segment;
        segment.set_segment(&head, &tail);

        EXPECT_EQ(x1, segment.get_head()->get_gx());
        EXPECT_EQ(y1, segment.get_head()->get_gy());
        EXPECT_EQ(x2, segment.get_tail()->get_gx());
        EXPECT_EQ(y2, segment.get_tail()->get_gy());
    }
}

TEST(SegmentTest, GetLength) {
    srand(time(NULL));
    int size = 1000000;
    for (int i = 0; i < size; i++) {
        double length = rand() % 100000;
        double x1 = rand() % 100;
        double y1 = rand() % 100;
        double x2 = x1;
        double y2 = y1 + length;

        Point head, tail;
        head.set_point(x1, y1);
        tail.set_point(x2, y2);

        Segment segment;
        segment.set_segment(&head, &tail);

        EXPECT_EQ(length, segment.get_length());
        EXPECT_TRUE(segment.get_length() >= 0);
    }

    for (int i = 0; i < size; i++) {
        double length = rand() % 100000;
        double x1 = rand() % 100;
        double y1 = rand() % 100;
        double x2 = x1 + length;
        double y2 = y1;

        Point head, tail;
        head.set_point(x1, y1);
        tail.set_point(x2, y2);

        Segment segment;
        segment.set_segment(&head, &tail);

        EXPECT_EQ(length, segment.get_length());
        EXPECT_TRUE(segment.get_length() >= 0);
    }
}
