#include <cstdio>
#include"RNA.h"
#include "gtest/gtest.h"


//test constructors


TEST(RNA, TestRNAConstructor) {
    Nucleotid n = C;
    RNA rn(1000, n);
    RNA rn2;
    bool tmp = true;
    for (size_t i = 0;i<rn.capacity();i++) {
        if (Nucleotid (rn[i]) != n) {
        tmp = false;
        break;
        }
    }

    EXPECT_EQ(1000, rn.capacity());
    EXPECT_EQ(0, rn2.getlength());
    EXPECT_TRUE(tmp);
}
//copy and "=="
TEST(TestCopy, TestCopyRNA) {
    Nucleotid n = C;
    RNA rn1(1000, n);                    //chenged to const operator[] - need to check
    RNA rn2(rn1);
    EXPECT_EQ(rn1.getlength(), rn2.getlength());
    EXPECT_EQ(rn1.capacity(), rn2.capacity());
    bool tmp = true;
    bool IsEqual = (rn1 == rn2);
    for (size_t i = 0;i<rn1.getlength();i++) {
        if (Nucleotid (rn1[i]) != Nucleotid (rn2[i])) {
            tmp = false;
            break;
        }
    }
    EXPECT_TRUE(tmp);
    EXPECT_TRUE(IsEqual);
}
//operator=
TEST(RNA, TestOperatorAssigningRNK) {
    Nucleotid n = C;
    RNA rn1(1000, n);
    RNA rn2;
    rn2 = rn1;
    EXPECT_EQ(rn1.getlength(), rn2.getlength());
    EXPECT_EQ(rn1.capacity(), rn2.capacity());
    bool tmp = true;
    bool IsEqual = (rn1 == rn2);
        for (size_t i = 0;i<rn1.getlength();i++) {
            if (Nucleotid (rn1[i]) != Nucleotid (rn2[i])) {
                tmp = false;
                break;
            }
        }
    EXPECT_TRUE(tmp);
    EXPECT_TRUE(IsEqual);
}
TEST(TestAssigning, TestAssigningThisToThis) {
    Nucleotid n = C;
    RNA rn1(1000, n);
    RNA rn2(rn1);
    rn1 = rn1;
    bool tmp = true;
    for (size_t i = 0;i<rn1.getlength();i++) {
        if (Nucleotid (rn1[i]) != Nucleotid (rn2[i])) {
            tmp = false;
            break;
        }
    }
    EXPECT_TRUE(rn1.getlength()== rn2.getlength());
    EXPECT_TRUE(tmp);
}



//void trim();
TEST(TrimTest, TestTrimRNK) {
    RNA rn2(1000, G);
    rn2=rn2.trim(0);
    EXPECT_EQ(rn2.getlength(),1);
    EXPECT_EQ(rn2.capacity(),1);
}

//run time
TEST(TestTime, RunTimeSmall) {
    RNA rn;
    int size = 10000;
//	size_t start = clock();
    for (int i = 0;i<size;i++) {
        rn[i] =G;//(Nucleotide)((i*i) % 4);
    }
//	size = 1000000000;
//	size_t end = clock();
//	std::cout << "runTime: " << end - start << "\n";
    EXPECT_TRUE(true);
}
//run time with 100000000
TEST(TestTime, RunTimeBig) {
    RNA rn;
    int size = 100000000;
//	size_t start = clock();
    for (int i = 0;i<size;i++) {
//		rn[i] = G;//(Nucleotide)((i*i) % 4);
}
    EXPECT_TRUE(true);
}

TEST(TestOperatorSum, TestCorrect) {
    size_t size1 = 1000;
    size_t size2 = 1000;
    Nucleotid n1 = C;
    Nucleotid n2 = G;

    RNA rna1(1000, C);

    RNA rna2(1000, n2);

    RNA rna3 = rna1 + rna2;

    bool isEqual = true;
    for (size_t i = 0;i<size1;i++) {
        if (Nucleotid (rna3[i]) != n1) {
            isEqual = false;
            break;
        }
    }
    if (isEqual) {
        for (size_t i = size1;i<size1+size2;i++) {
            if (Nucleotid (rna3[i]) != n2) {
                isEqual = false;
                break;
            }
         }
    }
    EXPECT_TRUE(isEqual);
}
TEST(TestOperatorSum, TestBoundaryCauseSizeOf) {
    size_t size1 = sizeof(size_t);
    size_t size2 = 3 * sizeof(size_t);
    Nucleotid n1 = C;
    Nucleotid n2 = G;
    RNA rn1(size1, n1);
    RNA rn2(size2, n2);
    RNA rn3 = rn1 + rn2;
    bool isEqual = true;
    for (size_t i = 0;i<size1;i++) {
        if (Nucleotid (rn3[i]) != n1) {
            isEqual = false;
            break;
        }
    }
    if (isEqual == true) {
        for (size_t i = size1;i<rn3.getlength();i++) {
            if (Nucleotid (rn3[i]) != n2) {
                isEqual = false;
                break;
            }
        }
    }
    EXPECT_TRUE(isEqual);
}
TEST(TestOperatorSum, TestBoundaryCauseSizeOfMinusOne) {
    size_t size1 = (2 * sizeof(size_t)) - 1;
    size_t size2 = 3 * sizeof(size_t);
    Nucleotid n1 = C;
    Nucleotid n2 = G;
    RNA rn1(size1, n1);
    RNA rn2(size2, n2);
    RNA rn3 = rn1 + rn2;
    bool isEqual = true;
    for (size_t i = 0;i<size1;i++) {
        if (Nucleotid (rn3[i]) != n1) {
            isEqual = false;
            break;
        }
    }
    if (isEqual) {
        for (size_t i = size1;i<rn3.getlength();i++) {
            if (Nucleotid (rn3[i]) != n2) {
                isEqual = false;
                break;
            }
        }
    }
EXPECT_TRUE(isEqual);
}
//operator!=

TEST(TestIsCompl, TestIsComplRnk) {
    size_t size = 10000;
    Nucleotid n = C;
    RNA rn1(size, n);
    RNA rn2 = rn1;
    EXPECT_FALSE(rn1.isComplementary(rn2));
}

TEST(TestSplit, TestSplitRunTime) {
    size_t size = 99;// 000;
    size_t splitNum = 2;//00;
    Nucleotid n = C;
    RNA rn1(size, n);

    RNA rn2 = rn1.split(splitNum);

    size_t tmpSize = rn1.getlength();
    bool tmp = true;
    for (size_t i = 0;i<splitNum;i++) {
        if (((Nucleotid) rn1[i]) != n) {
            tmp = false;
            break;
        }
    }
    for (size_t j = 0; j < size-splitNum; ++j) {
        if (Nucleotid (rn2[j]) != n) {
            tmp = false;
            break;
        }
    }


    EXPECT_TRUE(tmp);
}
//
// Created by Юлия Затолоцкая on 17.11.2019.
////
// Created by juliyazatolockaya on 18.11.2019.
//

