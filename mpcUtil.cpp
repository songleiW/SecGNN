#include "mpcUtil.h"
#include "typedef.h"
std::vector<long long> multiTriplesParty1, multiTriplesParty2;
std::vector<bool> multiBoolTriplesParty1, multiBoolTriplesParty2;
bool gParty1[DATALENGTH], pParty1[DATALENGTH], gParty2[DATALENGTH], pParty2[DATALENGTH],
        x1Bits[DATALENGTH], x2Bits[DATALENGTH], bitx[DATALENGTH - 1], res[2 * DATALENGTH];
int key12 = 45, key13 = 64, key23 = 65;
extern vector<inputNode> inputLayerParty1, inputLayerParty2;
extern vector<hiddenNode> hiddenLayerParty1, hiddenLayerParty2;
long long z1, z2, e, f, a1, b1, c1, a2, b2, c2,
        alpha1[NODENUM][1 + FEATURENUM], alpha2[NODENUM][1 + FEATURENUM], alpha3[NODENUM][1 + FEATURENUM];
long long numMutil=0,numMutilBool=0, lengthData=0;
void initMPC() {
    long long temp;
    string str;
    ifstream infile("../cora/triplesParty1.txt", ios::in);
    while (getline(infile, str)) {
        stringstream input(str);
        input >> temp;
        multiTriplesParty1.push_back(temp);
        input >> temp;
        multiTriplesParty1.push_back(temp);
        input >> temp;
        multiTriplesParty1.push_back(temp);


    }
    infile.close();
    ifstream infile1("../cora/triplesParty2.txt", ios::in);
    while (getline(infile1, str)) {
        stringstream input(str);
        input >> temp;
        multiTriplesParty2.push_back(temp);
        input >> temp;
        multiTriplesParty2.push_back(temp);
        input >> temp;
        multiTriplesParty2.push_back(temp);
    }
    infile1.close();


    bool tempBool;
    ifstream infile2("../cora/triplesBoolParty1.txt", ios::in);
    while (getline(infile2, str)) {
        stringstream input(str);
        input >> tempBool;
        multiBoolTriplesParty1.push_back(tempBool);
        input >> tempBool;
        multiBoolTriplesParty1.push_back(tempBool);
        input >> tempBool;
        multiBoolTriplesParty1.push_back(tempBool);


    }
    infile2.close();

    ifstream infile3("../cora/triplesBoolParty2.txt", ios::in);
    while (getline(infile3, str)) {
        stringstream input(str);
        input >> tempBool;
        multiBoolTriplesParty2.push_back(tempBool);
        input >> tempBool;
        multiBoolTriplesParty2.push_back(tempBool);
        input >> tempBool;
        multiBoolTriplesParty2.push_back(tempBool);
    }
    infile3.close();
}

void mpcMulti(vector<long long> x1, vector<long long> y1,
              vector<long long> x2, vector<long long> y2, long long *result) {

    numMutil+=x1.size();
    for (int i = 0; i < x1.size(); i++) {

        c1 = multiTriplesParty1[multiTriplesParty1.size() - 1];
        // multiTriplesParty1.pop_back();


        b1 = multiTriplesParty1[multiTriplesParty1.size() - 2];
        // multiTriplesParty1.pop_back();

        a1 = multiTriplesParty1[multiTriplesParty1.size() - 3];
        // multiTriplesParty1.pop_back();

        c2 = multiTriplesParty2[multiTriplesParty2.size() - 1];
        //multiTriplesParty2.pop_back();

        b2 = multiTriplesParty2[multiTriplesParty2.size() - 2];
        //multiTriplesParty2.pop_back();

        a2 = multiTriplesParty2[multiTriplesParty2.size() - 3];
        //multiTriplesParty2.pop_back();
        e = x1[i] - a1 + x2[i] - a2;
        f = y1[i] - b1 + y2[i] - b2;
        result[i] = (e * f + f * a1 + e * b1 + c1) / t;
        result[x1.size() + i] = (f * a2 + e * b2 + c2) / t;
    }
}

void mpcMulti(long long x1, long long y1,
              long long x2, long long y2, long long *result) {
    numMutil++;
    c1 = multiTriplesParty1[multiTriplesParty1.size() - 1];
    // multiTriplesParty1.pop_back();
    b1 = multiTriplesParty1[multiTriplesParty1.size() - 2];
    // multiTriplesParty1.pop_back();
    a1 = multiTriplesParty1[multiTriplesParty1.size() - 3];
    // multiTriplesParty1.pop_back();
    c2 = multiTriplesParty2[multiTriplesParty2.size() - 1];
    //multiTriplesParty2.pop_back();
    b2 = multiTriplesParty2[multiTriplesParty2.size() - 2];
    //multiTriplesParty2.pop_back();
    a2 = multiTriplesParty2[multiTriplesParty2.size() - 3];
    //multiTriplesParty2.pop_back();

    e = x1 - a1 + x2 - a2;
    f = y1 - b1 + y2 - b2;
    result[0] = (e * f + f * a1 + e * b1 + c1) / t;
    result[1] = (f * a2 + e * b2 + c2) / t;
}

void mpcMultiBool(bool *x1Begin, bool *y1Begin,
                  bool *x2Begin, bool *y2Begin, int length, bool *result) {

    bool a1, b1, c1, a2, b2, c2, e, f;
    numMutilBool+=length;

    for (int i = 0; i < length; i++) {

        c1 = multiBoolTriplesParty1[multiBoolTriplesParty1.size() - 1];
        // multiBoolTriplesParty1.pop_back();

        b1 = multiBoolTriplesParty1[multiBoolTriplesParty1.size() - 2];
        // multiBoolTriplesParty1.pop_back();

        a1 = multiBoolTriplesParty1[multiBoolTriplesParty1.size() - 3];
        // multiBoolTriplesParty1.pop_back();

        c2 = multiBoolTriplesParty2[multiBoolTriplesParty2.size() - 1];
        //multiBoolTriplesParty2.pop_back();

        b2 = multiBoolTriplesParty2[multiBoolTriplesParty2.size() - 2];
        //multiBoolTriplesParty2.pop_back();

        a2 = multiBoolTriplesParty2[multiBoolTriplesParty2.size() - 3];
        //multiBoolTriplesParty2.pop_back();

        e = *(x1Begin + i) ^ a1 ^ *(x2Begin + i) ^ a2;
        f = *(y1Begin + i) ^ b1 ^ *(y2Begin + i) ^ b2;

        *(result + i) = e & f ^ f & a1 ^ e & b1 ^ c1;
        *(result + length + i) = f & a2 ^ e & b2 ^ c2;
    }
}

void mpcMsb(long long x1, long long x2, bool *msb) {
    bitset<DATALENGTH> Bitx1(x1), Bitx2(x2);

    for (int i = 0; i < 64; ++i) {
        x1Bits[i] = Bitx1[i];
        x2Bits[i] = Bitx2[i];
    }


    mpcMultiBool(x1Bits, bitx, bitx, x2Bits, DATALENGTH - 1, res);

    gParty1[0] = false;
    gParty2[0] = false;

    for (int i = 0; i < DATALENGTH - 1; ++i) {
        gParty1[i + 1] = res[i];
        gParty2[i + 1] = res[DATALENGTH - 1 + i];

        pParty1[i + 1] = x1Bits[i];
        pParty2[i + 1] = x2Bits[i];
    }

    int length = DATALENGTH;

    for (int i = 0; i < 6; ++i) {

        mpcMultiBool(pParty1 + 1, gParty1, pParty2 + 1, gParty2, 1, res);
        gParty1[0] = res[0] ^ gParty1[1];
        gParty2[0] = res[1] ^ gParty2[1];

        int num = 1;
        for (int j = 2; j < length; j +=2) {
            mpcMultiBool(pParty1 + j + 1, gParty1 + j, pParty2 + j + 1, gParty2 + j, 1, res);

            gParty1[num] = res[0] ^ gParty1[j + 1];
            gParty2[num] = res[1] ^ gParty2[j + 1];

            mpcMultiBool(pParty1 + j + 1, pParty1 + j, pParty2 + j + 1, pParty2 + j, 1, res);

            pParty1[num] = res[0];
            pParty2[num++] = res[1];
        }
        length /= 2;
    }

    *msb = x1Bits[DATALENGTH - 1] ^ gParty1[0];
    *(msb + 1) = x2Bits[DATALENGTH - 1] ^ gParty2[0];
}

bool* mpcReLU(long long x1, long long x2, long long *relu) {

    bool msb[2];
    long long random;
    mpcMsb(x1, x2, msb);

    msb[0] = !msb[0];

    *relu = sqrt(rand());
    *(relu + 1) = (msb[0] ^ msb[1]) * x1 - *relu;

    random = sqrt(rand());

    *(relu + 1) += random;
    *relu += (msb[0] ^ msb[1]) * x2 - random;

    return msb;
}

void mpcMax(long long x1, long long x2, long long *max) {

    mpcReLU(max[0] - x1, max[1] - x2, max);
    *max += x1;
    *(max + 1) += x2;
}

void mpcExp(vector<long long> x1, vector<long long> x2, long long *result) {
    for (int i = 0; i < x1.size(); ++i) {
        x1[i] = t / 2 + x1[i] / pow(2, 9);
        x2[i] = t / 2 + x2[i] / pow(2, 9);
    }

    for (int i = 0; i < 8; ++i) {
        mpcMulti(x1, x1, x2, x2, result);
        for (int j = 0; j < x1.size(); ++j) {
            x1[j] = result[j];
            x2[j] = result[j + x1.size()];
        }
    }
    mpcMulti(x1, x1, x2, x2, result);
}

void mpcExp(long long x1, long long x2, long long *result) {

    x1 = t / 2 + x1 / pow(2, 9);
    x2 = t / 2 + x2 / pow(2, 9);

    for (int i = 0; i < 8; ++i) {
        mpcMulti(x1, x1, x2, x2, result);
        x1 = result[0];
        x2 = result[1];
    }
    mpcMulti(x1, x1, x2, x2, result);
}

void mpcReci(long long  x1, long long x2, long long *result){

    long long y1,y2;

    mpcExp(0.25*t-x1,0.25*t-x2,result);

    y1=3*result[0]+0.0015*t;
    y2=3*result[1]+0.0015*t;

    for (int i = 0; i < 13; ++i) {
        mpcMulti(y1,y1,y2,y2,result);
        mpcMulti(x1,result[0],x2,result[1],result);
        y1=2*y1-result[0];
        y2=2*y2-result[1];
    }
    result[0]=y1;
    result[1]=y2;
}

void mpcSquareRoot(long long x1, long long x2, long long *result) {
    long long y1, y2;

    mpcExp(0.5 - x1, 0.5 - x2, result);

    y1 = sqrt(3 * result[0] + 0.003);
    y2 = sqrt(3 * result[1] + 0.003);

    for (int i = 0; i < 17; ++i) {
        mpcMulti(y1, y1, y2, y2, result);
        mpcMulti(y1, result[0], y2, result[1], result);
        mpcMulti(x1, result[0], x2, result[1], result);

        y1 = 1.5 * y1 - result[0] / 2;
        y2 = 1.5 * y2 - result[1] / 2;
    }
    result[0] = y1;
    result[1] = y2;
}

void mpcArrayAccessInput(vector<vector<long long> > index1, vector<vector<long long> > index2,
                    vector<vector<long long> > &neibNumParty1, vector<vector<long long> > &neibNumParty2,
                    vector<vector<long long> > &neibFeatureParty1, vector<vector<long long> > &neibFeatureParty2) {
    lengthData+=2*NODENUM*(1+FEATURENUM);
    vector<inputNode> nodeParty1(inputLayerParty1.size()), nodeParty2(inputLayerParty1.size()),
            nodeParty11(inputLayerParty1.size()), nodeParty22(inputLayerParty1.size());
    long long r1, r2, random1, random2, random3, index;

    r1 = rand();
    r2 = rand();
    for (int i = 0; i < NODENUM; ++i) {
        srand(key12 + i);
        random1 = rand();
        srand(key13 + i);
        random2 = rand();
        srand(key23 + i);
        random3 = rand();

        alpha1[i][0] = random1 - random2;
        alpha2[i][0] = random2 - random3;
        alpha3[i][0] = random3 - random1;

        for (int j = 1; j < FEATURENUM + 1; ++j) {
            srand(key12 + i + j);
            random1 = rand();
            srand(key13 + i + j);
            random2 = rand();
            srand(key23 + i + j);
            random3 = rand();

            alpha1[i][j] = random1 - random2;
            alpha2[i][j] = random2 - random3;
            alpha3[i][j] = random3 - random1;
        }
    }

    //party1
    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r1) % NODENUM;

        nodeParty1[index] = inputLayerParty1[j];
        nodeParty1[index].neibNum += alpha1[index][0];

        for (int i = 0; i < FEATURENUM; ++i) {
            nodeParty1[index].feature[i] += alpha1[index][i + 1];
        }
    }

    //party3
    for (int j = 0; j < NODENUM; ++j) {
        nodeParty1[j].neibNum += alpha3[j][0];
        for (int i = 0; i < FEATURENUM; ++i) {
            nodeParty1[j].feature[i] += alpha3[j][i + 1];
        }
    }

    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r2) % NODENUM;

        nodeParty11[index] = nodeParty1[j];
    }

    //party2
    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r1) % NODENUM;

        nodeParty2[index] = inputLayerParty2[j];
        nodeParty2[index].neibNum += alpha2[index][0];

        for (int i = 0; i < FEATURENUM; ++i) {
            nodeParty2[index].feature[i] += alpha2[index][i + 1];
        }
    }

    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r2) % NODENUM;
        nodeParty22[index] = nodeParty2[j];
    }

    for (int i = 0; i < NODENUM; ++i) {

        for (int j = 0; j < index1[0].size(); ++j) {

            index=(index1[i][j] + index2[i][j] + r1 + r2) % NODENUM;

            neibNumParty1[i][j] = nodeParty11[index].neibNum;
            neibNumParty2[i][j] = nodeParty22[index].neibNum;

            for (int k = 0; k < FEATURENUM; ++k) {
                neibFeatureParty1[i] [j*FEATURENUM+ k] = nodeParty11[index].feature[k];
                neibFeatureParty2[i] [j*FEATURENUM+ k] = nodeParty22[index].feature[k];
            }
        }
    }
}


void mpcArrayAccessHidd(vector<long long > index1, vector<long long> index2,
                         vector<long long > &neibFeatureParty1, vector<long long> &neibFeatureParty2) {
    lengthData+=2*NODENUM*HIDDENNODE;
    vector<inputNode> nodeParty1(inputLayerParty1.size()), nodeParty2(inputLayerParty1.size()),
            nodeParty11(inputLayerParty1.size()), nodeParty22(inputLayerParty1.size());
    long long r1, r2, random1, random2, random3, index;

    r1 = rand();
    r2 = rand();
    for (int i = 0; i < NODENUM; ++i) {
        srand(key12 + i);
        random1 = rand();
        srand(key13 + i);
        random2 = rand();
        srand(key23 + i);
        random3 = rand();

        for (int j = 0; j < HIDDENNODE; ++j) {
            srand(key12 + i + j);
            random1 = rand();
            srand(key13 + i + j);
            random2 = rand();
            srand(key23 + i + j);
            random3 = rand();

            alpha1[i][j] = random1 - random2;
            alpha2[i][j] = random2 - random3;
            alpha3[i][j] = random3 - random1;
        }
    }

    //party1
    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r1) % NODENUM;
        nodeParty1[index] = inputLayerParty1[j];
        for (int i = 0; i < HIDDENNODE; ++i) {
            nodeParty1[index].embeddingFirstLayer[i] += alpha1[index][i];
        }
    }

    //party3
    for (int j = 0; j < NODENUM; ++j) {
        for (int i = 0; i < HIDDENNODE; ++i) {
            nodeParty1[j].embeddingFirstLayer[i] += alpha3[j][i];
        }
    }

    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r2) % NODENUM;

        nodeParty11[index] = nodeParty1[j];
    }

    //party2
    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r1) % NODENUM;

        nodeParty2[index] = inputLayerParty2[j];

        for (int i = 0; i < HIDDENNODE; ++i) {
            nodeParty2[index].embeddingFirstLayer[i] += alpha2[index][i];
        }
    }

    for (int j = 0; j < NODENUM; ++j) {
        index = (j + r2) % NODENUM;
        nodeParty22[index] = nodeParty2[j];
    }



    for (int j = 0; j < index1.size(); ++j) {

        index=(index1[j] + index2[j] + r1 + r2) % NODENUM;

        for (int k = 0; k < HIDDENNODE; ++k) {
            neibFeatureParty1[j*HIDDENNODE+ k] = nodeParty11[index].embeddingFirstLayer[k];
            neibFeatureParty2 [j*HIDDENNODE+ k] = nodeParty22[index].embeddingFirstLayer[k];
        }
    }
}

void mpcLog(vector<long long > x1, vector<long long> x2,long long *result){
    long long res[2*x1.size()];
    vector<long long > y0, y1,yy0(x1.size()),yy1(x1.size()),
    h1(x1.size()),h2(x1.size()),hh1(x1.size()),hh2(x1.size());

    for (int i = 0; i < x1.size(); ++i) {
        y0.push_back(-2*x1[i]-t/2.0);
        y1.push_back(-2*x2[i]-t/2.0);
    }

    mpcExp(y0,y1,res);

    for (int i = 0; i < x1.size(); ++i) {
        y0[i]=x1[i]/120.0-20*res[i]+1.5*t;
        y1[i]=x2[i]/120.0-20*res[i+x1.size()]+1.5*t;
    }

    copy( y0.begin(), y0.end(), yy0.begin() );
    copy( y1.begin(), y1.end(), yy1.begin() );

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < x1.size(); ++j) {
            y0[j]*=-1;
            y1[j]*=-1;
        }
        mpcExp(y0,y1,res);
        for (int j = 0; j < x1.size(); ++j) {
            y0[j]=res[j];
            y1[j]=res[j+x1.size()];
        }

        mpcMulti(y0,x1,y1,x2,res);

        for (int j = 0; j < x1.size(); ++j) {

            h1[j]=t/2.0-res[j];
            h2[j]=t/2.0-res[x1.size()+j];
            hh1[j]=t/2.0-res[j];
            hh2[j]=t/2.0-res[x1.size()+j];

            yy0[j]-=h1[j];
            yy1[j]-=h2[j];
        }

        for (int j = 0; j < 7; ++j) {

            mpcMulti(h1,hh1,h2,hh2,res);

            for (int k = 0; k < x1.size(); ++k) {
                hh1[k]=res[k];
                hh2[k]=res[k+x1.size()];

                yy0[k]-=hh1[k]/(j+2);
                yy1[k]-=hh2[k]/(j+2);
            }
        }

        copy( yy0.begin(), yy0.end(), y0.begin() );
        copy( yy1.begin(), yy1.end(), y1.begin() );
    }

    for (int i = 0; i < x1.size(); ++i) {
        result[i]=y0[i];
        result[i+x1.size()]=y1[i];
    }
}