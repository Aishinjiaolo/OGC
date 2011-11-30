#include "segment.h"
#include "gtest/gtest.h"

TEST(SegmentTest, SetSegment) {
    srand(time(NULL));
    for (int i = 0; i < 100; i++) {
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


