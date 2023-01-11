#include "typedef.h"
#include "util.h"
#include "test.h"
#include<ctime>
using namespace std;
extern long long lengthData,numMutilBool,numMutil;
int main() {
    int flag=0;
    clock_t startClock, endClock;
    startClock = clock();        //程序开始计时
    initMPC();
    cout << "---随机数传输完成---" << endl;
    readData(flag);
    cout << "---密文数据读取完成---" << endl;

    selectTrainSet(40);

    if (1) {

        Init(1);

        cout << "---权重初始化完成---" << endl;

        train();

        endClock = clock();        //程序结束用时
        double endtime = (double) (endClock - startClock) / CLOCKS_PER_SEC;
        long long onCom=Online();
        long long offCom=Offline();

        cout <<"Online comm: " <<onCom<<" MB"<<endl;
        cout <<"Offline comm: " <<offCom<<" MB"<<endl;

        cout << "Training time: " << endtime / 120.0+ onCom/(dataRate*60.0)<< " minutes" << endl;
        cout << "Offline Time: " << offCom/(dataRate*60.0)<< " minutes" << endl;

        writeModel(flag);
        cout << "---写入模型成功---" << endl;

    }
    else{
        readModel(flag);
        cout << "---密文模型读取完成---" << endl;
        Init(0);
        cout << "---初始化完成---" << endl;

        double accuracy = testData();

        endClock = clock();        //程序结束用时

        cout << "---准确率：" << accuracy << "---" << endl;

        double endtime = (double) (endClock - startClock) / CLOCKS_PER_SEC;
        long long onCom=Online();
        long long offCom=Offline();

        cout <<"Online comm: " <<onCom<<" MB"<<endl;
        cout <<"Offline comm: " <<offCom<<" MB"<<endl;

        cout << "Inference time: " << endtime / 120.0+ onCom/(dataRate*60.0)<< " minutes" << endl;
        cout << "Offline Time: " << offCom/(dataRate*60.0)<< " minutes" << endl;

    }
    return 0;
}

