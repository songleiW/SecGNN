#ifndef PLAINGCN_UTIL_H
#define PLAINGCN_UTIL_H

#include "typedef.h"
#include<climits>

using namespace std;

long long getRandom(int);

void Init(bool);

void readData(int);

void selectTrainSet(int);

void getNodeData(const string, vector<inputNode> &);

void getEdgeData(const string, vector<inputNode> &);

void getEdgeWeightData(const string, vector<inputNode> &);

void getNeibNum(const string, vector<inputNode> &);

void selectTrainSet(int);

void train();

void forwardFirstLayer();

void forwardSecondLayer(int);

void backward(int);

void testValidaSet();

void writeModel(int);

long long crossEntropy(std::vector<long long>, std::vector<long long>,
                       std::vector<long long>, std::vector<long long>);

void writeData(vector<inputNode>);

void sofMax(int);
long long Online();
long long Offline();
#endif //PLAINGCN_UTIL_H
