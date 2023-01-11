#include "util.h"
#include "test.h"
vector<inputNode> inputLayerParty1, inputLayerParty2;
vector<hiddenNode> hiddenLayerParty1, hiddenLayerParty2;
vector<outputNode> outputLayerParty1, outputLayerParty2;
static long long lastError = 100 * t;
static long long error = 100 * t;
vector<int> indexTrainSet,indexValidSet;
extern long long lengthData,numMutilBool,numMutil;
vector<vector<long long> > indexParty1(NODENUM), indexParty2(NODENUM);
int dummyNeibNum;
double accuracy;
void readData(int flag) {
    if (flag==0){
        getNodeData("../cora/contentParty1.txt", inputLayerParty1);
        getEdgeData("../cora/citesParty1.txt", inputLayerParty1);
        getEdgeWeightData("../cora/edgeWeightParty1.txt", inputLayerParty1);
        getNeibNum("../cora/NeibNumParty1.txt", inputLayerParty1);

        getNodeData("../cora/contentParty2.txt", inputLayerParty2);
        getEdgeData("../cora/citesParty2.txt", inputLayerParty2);
        getEdgeWeightData("../cora/edgeWeightParty2.txt", inputLayerParty2);
        getNeibNum("../cora/NeibNumParty2.txt", inputLayerParty2);

    }else if(flag==1){
        getNodeData("../citeseer/contentParty1.txt", inputLayerParty1);
        getEdgeData("../citeseer/citesParty1.txt", inputLayerParty1);
        getEdgeWeightData("../citeseer/edgeWeightParty1.txt", inputLayerParty1);
        getNeibNum("../citeseer/NeibNumParty1.txt", inputLayerParty1);

        getNodeData("../citeseer/contentParty2.txt", inputLayerParty2);
        getEdgeData("../citeseer/citesParty2.txt", inputLayerParty2);
        getEdgeWeightData("../citeseer/edgeWeightParty2.txt", inputLayerParty2);
        getNeibNum("../citeseer/NeibNumParty2.txt", inputLayerParty2);
    }else{
        getNodeData("../pubmed/contentParty1.txt", inputLayerParty1);
        getEdgeData("../pubmed/citesParty1.txt", inputLayerParty1);
        getEdgeWeightData("../pubmed/edgeWeightParty1.txt", inputLayerParty1);
        getNeibNum("../pubmed/NeibNumParty1.txt", inputLayerParty1);

        getNodeData("../pubmed/contentParty2.txt", inputLayerParty2);
        getEdgeData("../pubmed/citesParty2.txt", inputLayerParty2);
        getEdgeWeightData("../pubmed/edgeWeightParty2.txt", inputLayerParty2);
        getNeibNum("../pubmed/NeibNumParty2.txt", inputLayerParty2);
    }
    dummyNeibNum=inputLayerParty1[0].neibors.size();
}

void getNodeData(const string fileName, vector<inputNode> &inputLayer) {
    ifstream infile(fileName, ios::in);

    if (!infile) {  
        cerr << "open error." << endl;
        exit(1); 
    }

    double data;
    string str; 
    while (getline(infile, str)) {
        inputNode newNode;
        stringstream input(str);
        int num=0;
        for (int i = 0; i < FEATURENUM; ++i) { //features
            input >> data;
            newNode.feature.push_back(data*t);
        }

        for (int i = 0; i < OUTNODE; ++i) {
            input >> data;
            newNode.trueValue.push_back(data);
        }
        inputLayer.push_back(newNode);
    }
    infile.close();
}

void getEdgeData(const string fileName, vector<inputNode> &inputLayer) {
    ifstream input(fileName);
    if (!input) {
        cerr << "open error." << endl;
        exit(1); // 退出程序
    }
    string buff;
    int num = 0;
    while (getline(input, buff)) {
        char *datas = (char *) buff.c_str();
        const char *spilt = " ";
        // strtok字符串拆分函数
        char *data = strtok(datas, spilt);
        while (data != NULL) {
            // atof是stdlib头文件下转化字符串为数字的函数
            inputLayer[num].neibors.push_back(atof(data));
            // NULL代表从上次没拆分完地方继续拆
            data = strtok(NULL, spilt);
        }
        num++;
    }
    input.close();
}

void getEdgeWeightData(const string fileName, vector<inputNode> &inputLayer) {
    ifstream input(fileName);
    if (!input) {
        cerr << "open error." << endl;
        exit(1); // 退出程序
    }
    string buff;
    int num = 0;
    while (getline(input, buff)) {
        char *datas = (char *) buff.c_str();
        const char *spilt = " ";
        // strtok字符串拆分函数
        char *data = strtok(datas, spilt);
        while (data != NULL) {
            // atof是stdlib头文件下转化字符串为数字的函数
            inputLayer[num].edgeWeight.push_back(atof(data) * t);
            // NULL代表从上次没拆分完地方继续拆
            data = strtok(NULL, spilt);
        }
        num++;
    }
    input.close();
}

void getNeibNum(const string fileName, vector<inputNode> &inputLayer) {
    ifstream input(fileName);
    if (!input) {
        cerr << "open error." << endl;
        exit(1); // 退出程序
    }
    string buff;
    int num = 0;
    while (getline(input, buff)) {
        char *datas = (char *) buff.c_str();
        const char *spilt = " ";
        // strtok字符串拆分函数
        char *data = strtok(datas, spilt);
        inputLayer[num].neibNum = atof(data) * t;
        num++;
    }
    input.close();
}

void Init(bool flag) {

    long long randomWeight, random,res[2],reci[2];
    for (int i = 0; i < NODENUM; ++i) {
        long long numFeartParty1=0,numFeartParty2=0;
        for (int j = 0; j < FEATURENUM; ++j) {
            numFeartParty1+=inputLayerParty1[i].feature[j];
            numFeartParty2+=inputLayerParty2[i].feature[j];
        }
        mpcReci(numFeartParty1,numFeartParty2,reci);

        for (int j = 0; j < FEATURENUM; ++j) {
            mpcMulti(inputLayerParty1[i].feature[j],reci[0],
                     inputLayerParty2[i].feature[j],reci[1],res);
            inputLayerParty1[i].feature[j]=res[0];
            inputLayerParty2[i].feature[j]=res[1];
        }

    }

    for (int i = 0; i < HIDDENNODE; ++i) {
        hiddenNode newNode1, newNode2;
        for (int j = 0; j < FEATURENUM; ++j) {
            newNode1.wDelta.push_back(0);
            newNode2.wDelta.push_back(0);

            if (flag) {
                randomWeight = getRandom(HIDDENNODE);
                random = sqrt(rand());

                newNode1.weights.push_back(random);
                newNode2.weights.push_back(randomWeight - random);
                lengthData++;
            }
        }
        hiddenLayerParty1.push_back(newNode1);
        hiddenLayerParty2.push_back(newNode2);
    }

    for (int i = 0; i < OUTNODE; ++i) {
        outputNode newNode1, newNode2;
        for (int j = 0; j < HIDDENNODE; ++j) {
            newNode1.wDelta.push_back(0);
            newNode2.wDelta.push_back(0);

            if (flag) {
                randomWeight = getRandom(OUTNODE);
                random = sqrt(rand());
                newNode1.weights.push_back(random);
                newNode2.weights.push_back(randomWeight - random);
                lengthData++;
            }
        }

        outputLayerParty1.push_back(newNode1);
        outputLayerParty2.push_back(newNode2);
    }

    for (int i = 0; i < NODENUM; ++i) {
        for (int j = 0; j < HIDDENNODE; ++j) {
            inputLayerParty1[i].embeddingFirstLayer.push_back(0);
            inputLayerParty2[i].embeddingFirstLayer.push_back(0);
            inputLayerParty1[i].msb.push_back(false);
            inputLayerParty2[i].msb.push_back(false);
        }

        for (int j = 0; j < OUTNODE; ++j) {
            inputLayerParty1[i].embeddingSecondLayer.push_back(0);
            inputLayerParty2[i].embeddingSecondLayer.push_back(0);

        }
    }

    vector<vector<long long> >
            neibNumParty1(NODENUM,vector<long long > (dummyNeibNum)),
    neibNumParty2(NODENUM,vector<long long > (dummyNeibNum)),
    neibFeatureParty1(NODENUM,vector<long long > (dummyNeibNum*FEATURENUM)),
    neibFeatureParty2(NODENUM,vector<long long > (dummyNeibNum*FEATURENUM));

    for (int i = 0; i < NODENUM; ++i) {
        for (int j = 0; j < dummyNeibNum; ++j) {
            indexParty1[i].push_back(inputLayerParty1[i].neibors[j]);
            indexParty2[i].push_back(inputLayerParty2[i].neibors[j]);
        }
    }

    mpcArrayAccessInput(indexParty1, indexParty2,neibNumParty1, neibNumParty2,neibFeatureParty1, neibFeatureParty2);

    for (int i = 0; i < NODENUM; ++i) {
        long long degree[2], res[2], degreeReci[2], degreeNeib[2];

        mpcReci(inputLayerParty1[i].neibNum, inputLayerParty2[i].neibNum, degree);

        inputLayerParty1[i].coef.push_back(degree[0]);
        inputLayerParty2[i].coef.push_back(degree[1]);

        for (int j = 0; j < FEATURENUM; ++j) {
            mpcMulti(degree[0],inputLayerParty1[i].feature[j],degree[1],inputLayerParty2[i].feature[j],res);
            inputLayerParty1[i].embeddingFirst.push_back(res[0]);
            inputLayerParty2[i].embeddingFirst.push_back(res[1]);
        }

        mpcSquareRoot(inputLayerParty1[i].neibNum, inputLayerParty2[i].neibNum, degree);

        for (int j = 0; j < dummyNeibNum; ++j) {

            mpcSquareRoot(neibNumParty1[i][j], neibNumParty2[i][j], degreeNeib);

            mpcMulti(degree[0], degreeNeib[0], degree[1], degreeNeib[1], degreeReci);

            mpcMulti(degreeReci[0], inputLayerParty1[i].edgeWeight[j],
                     degreeReci[1], inputLayerParty2[i].edgeWeight[j], res);

            inputLayerParty1[i].coef.push_back(res[0]);
            inputLayerParty2[i].coef.push_back(res[1]);

            for (int k = 0; k < FEATURENUM; ++k) {
                long long result[2];
                mpcMulti(res[0], neibFeatureParty1[i][j * FEATURENUM + k],
                         res[1], neibFeatureParty2[i][j * FEATURENUM + k], result);

                inputLayerParty1[i].embeddingFirst[k] += result[0];
                inputLayerParty2[i].embeddingFirst[k] += result[1];
            }
        }
    }

    for (int i = 0; i < NODENUM; ++i) {
        for (int j = 0; j < HIDDENNODE; ++j) {
            inputLayerParty1[i].embeddingSecond.push_back(0);
            inputLayerParty2[i].embeddingSecond.push_back(0);
        }
    }
}

void selectTrainSet(int numEachClass){
    int max=0;
    for (int i = 0; i < OUTNODE; ++i) {
        for (int j = 0,num=0;num<numEachClass; ++j) {
            if (inputLayerParty1[j].trueValue[i]+inputLayerParty2[j].trueValue[i]==1){
                if (num<numEachClass){
                    indexTrainSet.push_back(j);
                }
                num++;
            }
        }
    }

    for (int i = 0; i < NODENUM; ++i) {
        int nCount = std::count(indexTrainSet.begin(), indexTrainSet.end(), i);
        if (nCount <= 0)
        {
            indexValidSet.push_back(i);
        }
        if (indexValidSet.size()>=500){
            break;
        }
    }

}

void train() {
    int epoch = 0;
    int flag = 0;

    while (epoch <= MAXEPOC) {

        if (abs(lastError-error)<0.02*t*indexTrainSet.size()){
            flag++;
        }
        else{
            flag=0;
        }

        if (flag>=WINDOWSIZE){
            break;
        }
        lastError=error;
        error = 0;
        //初始化隐藏层的变化量
        for (int i = 0; i < HIDDENNODE; ++i) {
            hiddenLayerParty1[i].wDelta.assign(FEATURENUM, 0);
            hiddenLayerParty2[i].wDelta.assign(FEATURENUM, 0);
        }

        //初始化输出层的变化量
        for (int i = 0; i < OUTNODE; ++i) {
            outputLayerParty1[i].wDelta.assign(HIDDENNODE, 0);
            outputLayerParty2[i].wDelta.assign(HIDDENNODE, 0);
        }

        forwardFirstLayer();

        for (int iter = 0; iter < indexTrainSet.size(); iter++) {

            forwardSecondLayer(indexTrainSet[iter]);
            backward(indexTrainSet[iter]);
        }

        // 修改加权
        for (int i = 0; i < HIDDENNODE; ++i) {
            for (int j = 0; j < FEATURENUM; ++j) {
                hiddenLayerParty1[i].weights[j] -=
                        LEARNINGRATE  * hiddenLayerParty1[i].wDelta[j];
                hiddenLayerParty2[i].weights[j] -=
                        LEARNINGRATE * hiddenLayerParty2[i].wDelta[j];
            }
        }


        for (int i = 0; i < OUTNODE; ++i) {
            for (int j = 0; j < HIDDENNODE; ++j) {
                outputLayerParty1[i].weights[j] -=
                        LEARNINGRATE  * outputLayerParty1[i].wDelta[j];
                outputLayerParty2[i].weights[j] -=
                        LEARNINGRATE  * outputLayerParty2[i].wDelta[j];
            }
        }

        cout << "epoch:" << ++epoch << "\t\t" << "error: " << 1.0 * error / (t*indexTrainSet.size()) << endl;

        testValidaSet();
    }
}

void forwardFirstLayer() {

    //初始化第一层的输出
    for (int i = 0; i < NODENUM; ++i) {
        inputLayerParty1[i].embeddingFirstLayer.assign(HIDDENNODE, 0);
        inputLayerParty2[i].embeddingFirstLayer.assign(HIDDENNODE, 0);
    }

    //计算节点第一层的embedding
    for (int i = 0; i < NODENUM; ++i) {
        for (int j = 0; j < HIDDENNODE; ++j) {

            long long res[2][FEATURENUM];

            mpcMulti(inputLayerParty1[i].embeddingFirst, hiddenLayerParty1[j].weights,
                     inputLayerParty2[i].embeddingFirst, hiddenLayerParty2[j].weights, *res);

            for (int k = 0; k < FEATURENUM; ++k) {
                inputLayerParty1[i].embeddingFirstLayer[j] += res[0][k];
                inputLayerParty2[i].embeddingFirstLayer[j] += res[1][k];
            }

            long long result[2];
            bool msb[2];
            *msb=mpcReLU(inputLayerParty1[i].embeddingFirstLayer[j], inputLayerParty2[i].embeddingFirstLayer[j], result);
            inputLayerParty1[i].embeddingFirstLayer[j] = result[0];
            inputLayerParty2[i].embeddingFirstLayer[j] = result[1];

            inputLayerParty1[i].msb[j]=*msb;
            inputLayerParty2[i].msb[j]=*(msb+1);
        }
    }
}

void forwardSecondLayer(int index){

    long long res[2];

    for (int j = 0; j < HIDDENNODE; ++j) {

        mpcMulti(inputLayerParty1[index].coef[0],inputLayerParty1[index].embeddingFirstLayer[j],
                 inputLayerParty2[index].coef[0],inputLayerParty2[index].embeddingFirstLayer[j],res);

        inputLayerParty1[index].embeddingSecond[j]=res[0];
        inputLayerParty2[index].embeddingSecond[j]=res[1];
    }

    vector<long long > neibFeatureParty1(dummyNeibNum*HIDDENNODE), neibFeatureParty2(dummyNeibNum*HIDDENNODE);

    mpcArrayAccessHidd(indexParty1[index], indexParty2[index],neibFeatureParty1, neibFeatureParty2);

    for (int j = 0; j < dummyNeibNum; ++j) {

        for (int k = 0; k < HIDDENNODE; ++k) {

            mpcMulti(neibFeatureParty1[j*HIDDENNODE+k],inputLayerParty1[index].coef[j+1],
                     neibFeatureParty2[j*HIDDENNODE+k],inputLayerParty2[index].coef[j+1],res);

            inputLayerParty1[index].embeddingSecond[k]+=res[0];
            inputLayerParty2[index].embeddingSecond[k]+=res[1];
        }
    }

    for (int j = 0; j < OUTNODE; ++j) {
        long long res[2][HIDDENNODE] ;
        mpcMulti(inputLayerParty1[index].embeddingSecond,outputLayerParty1[j].weights,
                 inputLayerParty2[index].embeddingSecond,outputLayerParty2[j].weights,*res);

        inputLayerParty1[index].embeddingSecondLayer[j]=0;
        inputLayerParty2[index].embeddingSecondLayer[j]=0;

        for (int k = 0; k < HIDDENNODE; ++k) {
            inputLayerParty1[index].embeddingSecondLayer[j]+=res[0][k];
            inputLayerParty2[index].embeddingSecondLayer[j]+=res[1][k];
        }
    }

    sofMax(index);
}

void backward(int index) {

    error += crossEntropy(inputLayerParty1[index].trueValue, inputLayerParty1[index].embeddingSecondLayer,
                          inputLayerParty2[index].trueValue, inputLayerParty2[index].embeddingSecondLayer);

    for (int i = 0; i < HIDDENNODE; ++i) {
        for (int j = 0; j < OUTNODE; ++j) {

            vector<long long> tempVectorx1, tempVectorx2, tempVectory1, tempVectory2;

            tempVectorx1.push_back(
                    inputLayerParty1[index].embeddingSecondLayer[j] - inputLayerParty1[index].trueValue[j] * t);
            tempVectorx2.push_back(
                    inputLayerParty2[index].embeddingSecondLayer[j] - inputLayerParty2[index].trueValue[j] * t);
            tempVectory1.push_back(inputLayerParty1[index].embeddingSecond[i]);
            tempVectory2.push_back(inputLayerParty2[index].embeddingSecond[i]);

            long long res[2][1];

            mpcMulti(tempVectorx1, tempVectory1, tempVectorx2, tempVectory2, *res);

            outputLayerParty1[j].wDelta[i] += res[0][0];
            outputLayerParty2[j].wDelta[i] += res[1][0];
        }
    }

    for (int i = 0; i < FEATURENUM; ++i) {
        for (int j = 0; j < HIDDENNODE; ++j) {
            if (inputLayerParty1[index].msb[j]^inputLayerParty2[index].msb[j]==false){
               // continue;
            }
            for (int k = 0; k < OUTNODE; ++k) {
                vector<long long> tempVectorx1, tempVectorx2, tempVectory1, tempVectory2;

                tempVectorx1.push_back(inputLayerParty1[index].embeddingFirst[i]);
                tempVectory1.push_back(
                        inputLayerParty1[index].embeddingSecondLayer[k] - inputLayerParty1[index].trueValue[k] * t);
                tempVectorx2.push_back(inputLayerParty2[index].embeddingFirst[i]);
                tempVectory2.push_back(
                        inputLayerParty2[index].embeddingSecondLayer[k] - inputLayerParty2[index].trueValue[k] * t);

                long long res[2][1];

                mpcMulti(tempVectorx1, tempVectory1, tempVectorx2, tempVectory2, *res);

                tempVectory1[0] = outputLayerParty1[k].weights[j];
                tempVectory2[0] = outputLayerParty2[k].weights[j];
                tempVectorx1[0] = res[0][0];
                tempVectorx2[0] = res[1][0];

                mpcMulti(tempVectorx1, tempVectory1, tempVectorx2, tempVectory2, *res);

                hiddenLayerParty1[j].wDelta[i] += res[0][0];
                hiddenLayerParty2[j].wDelta[i] += res[1][0];

            }
        }
    }
}

void testValidaSet(){

    forwardFirstLayer();

    int rightNum=0;
    for (int i = 0; i < indexValidSet.size(); ++i) {
        if (test(indexValidSet[i])) {
            rightNum++;
        }
    }
    accuracy=1.0*rightNum/indexValidSet.size();
    cout<<"Accuracy of validation set: "<<accuracy<<endl;
}

void writeModel(int flag) {
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
    ofstream outfileParty1(file1, ios::trunc);
    ofstream outfileParty2(file2, ios::trunc);

    for (int i = 0; i < HIDDENNODE; ++i) {
        for (int j = 0; j < FEATURENUM; ++j) {
            outfileParty1 << hiddenLayerParty1[i].weights[j];
            outfileParty1 << " ";

            outfileParty2 << hiddenLayerParty2[i].weights[j];
            outfileParty2 << " ";

        }
        outfileParty1 << endl;
        outfileParty2 << endl;
    }

    for (int i = 0; i < OUTNODE; ++i) {
        for (int j = 0; j < HIDDENNODE; ++j) {
            outfileParty1 << outputLayerParty1[i].weights[j];
            outfileParty1 << " ";

            outfileParty2 << outputLayerParty2[i].weights[j];
            outfileParty2 << " ";
        }
        outfileParty1 << endl;
        outfileParty2 << endl;
    }
}

long long crossEntropy(std::vector<long long> trueValue, std::vector<long long> predValue,
                       std::vector<long long> trueValue1, std::vector<long long> predValue1) {

    long long res[OUTNODE*2], result=0;

    mpcLog(predValue,predValue1,res);

    for (int i = 0; i < OUTNODE; ++i) {
        predValue[i]=res[i];
        predValue1[i]=res[i+OUTNODE];
        trueValue[i]*=t;
        trueValue1[i]*=t;
    }

    mpcMulti(predValue,trueValue, predValue1,trueValue1,res);

    for (int i = 0; i < OUTNODE*2; ++i) {
        result-=res[i];
    }

    if (isnan(result)) {
        exit(2);
    }
    return result;
}

void sofMax(int index){

    long long MAX[2]={LONG_MIN/3,LONG_MIN/3},res[2*OUTNODE],reci[2],sum[2];
    sum[0]=0;
    sum[1]=0;
    for (int i = 0; i < OUTNODE; ++i) {
        mpcMax(inputLayerParty1[index].embeddingSecondLayer[i],
               inputLayerParty2[index].embeddingSecondLayer[i],MAX);
    }
    for (int i = 0; i < OUTNODE; ++i) {

        inputLayerParty1[index].embeddingSecondLayer[i]-=MAX[0];
        inputLayerParty2[index].embeddingSecondLayer[i]-=MAX[1];
    }

    mpcExp(inputLayerParty1[index].embeddingSecondLayer,
           inputLayerParty2[index].embeddingSecondLayer,res);

    for (int i = 0; i < OUTNODE; ++i) {
        sum[0]+=res[i];
        sum[1]+=res[i+OUTNODE];
    }

    mpcReci(sum[0],sum[1],reci);

    for (int i = 0; i < OUTNODE; ++i) {

        mpcMulti(res[i],reci[0],res[i+OUTNODE],reci[1],sum);
        inputLayerParty1[index].embeddingSecondLayer[i]=sum[0];
        inputLayerParty2[index].embeddingSecondLayer[i]=sum[1];
    }
}

long long getRandom(int num) {
    return t * ((2.0 * (double) rand() / RAND_MAX) / sqrt(num) - 1 / sqrt(num));
}

long long Online(){
    long long length;
    length=lengthData;
    length+=numMutil*2*2;
    length+=numMutilBool*2*2/64;
    return length/(1024*1024*8);
}

long long Offline(){
    long long length;
    length=numMutil*6;
    length+=numMutilBool*6/64;
    return length/(1024*1024*8);
}
