
#include <vector>
#include "iostream"
#include "typedef.h"
using namespace std;
void initMPC();

void mpcMulti(vector<long long>, vector<long long>,
              vector<long long>, vector<long long>, long long *);

void mpcMulti(long long, long long,
              long long, long long, long long *);

void mpcMsb(long long, long long, bool *);

void mpcMultiBool(bool *, bool *, bool *, bool *, int, bool *);

bool* mpcReLU(long long, long long, long long *);

void mpcExp(vector<long long>, vector<long long>, long long *);

void mpcExp(long long, long long, long long *result);

void mpcMax(long long, long long, long long *);

void mpcReci(long long, long long, long long *);

void mpcSquareRoot(long long, long long, long long *);

void mpcArrayAccessInput(vector<vector<long long> >, vector<vector<long long> >,
                    vector<vector<long long> > &, vector<vector<long long> > &,
                    vector<vector<long long> > &,vector<vector<long long> > &);
void mpcArrayAccessHidd(vector<long long > , vector<long long> ,
                        vector<long long > &, vector<long long> &);
void mpcLog(vector<long long > , vector<long long> ,long long *);
#ifndef ENCRYPTGCN_MPCUTIL_H
#define ENCRYPTGCN_MPCUTIL_H
#endif //ENCRYPTGCN_MPCUTIL_H
