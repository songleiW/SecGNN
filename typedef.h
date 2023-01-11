#ifndef PLAINGCN_TYPEDEF_H
#define PLAINGCN_TYPEDEF_H

#include <stdio.h>
#include <string>
#include "vector"
#include "map"
#include "cstring"
#include "cstdlib"
#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include "limits"
#include "mpcUtil.h"
#include <unistd.h>
#include<bitset>
#include "gmp.h"
#include "gmpxx.h"

#define HIDDENNODE 16
#define LEARNINGRATE 0.2   // 学习速率（注意：越高虽然越快 也容易误差较大）
#define MAXEPOC 200
#define WINDOWSIZE 5
#define t 32768
#define DATALENGTH 64
#define dataRate 100
#define NODENUM 2708//19717//3312//
#define OUTNODE 7//3//6//
#define FEATURENUM 1433//501//3703//
using namespace std;


typedef struct inputNode {

    int neibNum;
    std::vector<long long> neibors, edgeWeight;
    std::vector<long long > coef;
    std::vector<long long> feature;
    std::vector<long long> embeddingFirst, embeddingSecond; //未激活的
    std::vector<long long> embeddingFirstLayer, embeddingSecondLayer; //加权激活过的
    std::vector<long long> trueValue;
    std::vector<bool> msb;
};

typedef struct hiddenNode {
    std::vector<long long> weights, wDelta;
};

typedef struct outputNode {
    std::vector<long long> weights, wDelta;
};


#endif //PLAINGCN_TYPEDEF_H
