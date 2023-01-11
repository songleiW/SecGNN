#include "test.h"
#include "typedef.h"
#include "util.h"

extern std::vector<hiddenNode> hiddenLayerParty1, hiddenLayerParty2;
extern std::vector<outputNode> outputLayerParty1, outputLayerParty2;
extern std::vector<inputNode> inputLayerParty1, inputLayerParty2;
extern vector<int> indexTrainSet;

void readModel(int flag) {
    string file1,file2;
    if (flag==0) {
        file1="../cora/modelParty1.txt";
        file2="../cora/modelParty2.txt";
    }else if(flag==1){
        file1="../citeseer/modelParty1.txt";
        file2="../citeseer/modelParty2.txt";
    }else{
        file1="../pubmed/modelParty1.txt";
        file2="../pubmed/modelParty2.txt";
    }

    ifstream infile(file1, ios::in);
    if (!infile) {  // 判断文件是否存在
        cerr << "open error." << endl;
        exit(1); // 退出程序
    }
    string str; // 定义字符数组用来接受读取一行的数据
    int weight;
    for (int i = 0; i < HIDDENNODE; ++i) {
        getline(infile, str);
        stringstream input(str);
        hiddenNode newNode;
        for (int j = 0; j < FEATURENUM; ++j) {
            input >> weight;
            newNode.weights.push_back(weight);
        }
        hiddenLayerParty1.push_back(newNode);
    }

    for (int i = 0; i < OUTNODE; ++i) {
        getline(infile, str);
        stringstream input(str);
        outputNode newNode;
        for (int j = 0; j < HIDDENNODE; ++j) {
            input >> weight;
            newNode.weights.push_back(weight);
        }
        outputLayerParty1.push_back(newNode);
    }
    infile.close();

    ifstream infile2(file2, ios::in);
    if (!infile) {  // 判断文件是否存在
        cerr << "open error." << endl;
        exit(1); // 退出程序
    }
    for (int i = 0; i < HIDDENNODE; ++i) {
        getline(infile2, str);
        stringstream input(str);
        hiddenNode newNode;
        for (int j = 0; j < FEATURENUM; ++j) {
            input >> weight;
            newNode.weights.push_back(weight);
        }
        hiddenLayerParty2.push_back(newNode);
    }

    for (int i = 0; i < OUTNODE; ++i) {
        getline(infile2, str);
        stringstream input(str);
        outputNode newNode;
        for (int j = 0; j < HIDDENNODE; ++j) {
            input >> weight;
            newNode.weights.push_back(weight);
        }
        outputLayerParty2.push_back(newNode);
    }
    infile2.close();


}

bool test(int index) {
    long long degree[2], res[2], degreeReci[2], degreeNeib[2];

    mpcReci(inputLayerParty1[index].neibNum, inputLayerParty2[index].neibNum, degree);

    for (int j = 0; j < HIDDENNODE; ++j) {

        mpcMulti(degree[0], inputLayerParty1[index].embeddingFirstLayer[j],
                 degree[1], inputLayerParty2[index].embeddingFirstLayer[j], res);

        inputLayerParty1[index].embeddingSecond[j] = res[0];
        inputLayerParty2[index].embeddingSecond[j] = res[1];
    }

    mpcSquareRoot(inputLayerParty1[index].neibNum, inputLayerParty2[index].neibNum, degree);

    for (int j = 0; j < inputLayerParty1[index].neibors.size(); ++j) {

        int neibId = inputLayerParty1[index].neibors[j] + inputLayerParty2[index].neibors[j];

        mpcSquareRoot(inputLayerParty1[neibId].neibNum,
                      inputLayerParty2[neibId].neibNum, degreeNeib);

        mpcMulti(degree[0], degreeNeib[0], degree[1], degreeNeib[1], degreeReci);

        for (int k = 0; k < HIDDENNODE; ++k) {

            mpcMulti(degreeReci[0], inputLayerParty1[neibId].embeddingFirstLayer[k],
                     degreeReci[1], inputLayerParty2[neibId].embeddingFirstLayer[k], res);

            mpcMulti(res[0], inputLayerParty1[index].edgeWeight[j],
                     res[1], inputLayerParty2[index].edgeWeight[j], res);

            inputLayerParty1[index].embeddingSecond[k] += res[0];
            inputLayerParty2[index].embeddingSecond[k] += res[1];
        }
    }

    for (int j = 0; j < OUTNODE; ++j) {
        long long res[2][HIDDENNODE];
        mpcMulti(inputLayerParty1[index].embeddingSecond, outputLayerParty1[j].weights,
                 inputLayerParty2[index].embeddingSecond, outputLayerParty2[j].weights, *res);

        inputLayerParty1[index].embeddingSecondLayer[j] = 0;
        inputLayerParty2[index].embeddingSecondLayer[j] = 0;

        for (int k = 0; k < HIDDENNODE; ++k) {
            inputLayerParty1[index].embeddingSecondLayer[j] += res[0][k];
            inputLayerParty2[index].embeddingSecondLayer[j] += res[1][k];
        }
    }

    sofMax(index);

    int loc = 0;
    long long temp = 0;
    for (int i = 0; i < OUTNODE; ++i) {
        if (inputLayerParty1[index].embeddingSecondLayer[i] + inputLayerParty2[index].embeddingSecondLayer[i] > temp) {
            loc = i;
            temp = inputLayerParty1[index].embeddingSecondLayer[i] + inputLayerParty2[index].embeddingSecondLayer[i];
        }
    }
    return inputLayerParty1[index].trueValue[loc] + inputLayerParty2[index].trueValue[loc] == 1;
}


double testData() {

    forwardFirstLayer();

    double rightNum = 0,num=0;
    for (int i = 0; i < NODENUM; ++i) {
        if (num>=1000){
            break;
        }
        if (std::find(indexTrainSet.begin(), indexTrainSet.end(), i) != indexTrainSet.end()) {
            continue;
        }
        if (test(i)) {
            rightNum++;
        }
        num++;
    }

    return rightNum * 1.0 / 1000;
}