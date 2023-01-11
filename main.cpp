#include "typedef.h"
#include "util.h"
#include "test.h"
#include<ctime>
using namespace std;
extern long long lengthData,numMutilBool,numMutil;
int main() {
    int flag=0;
    clock_t startClock, endClock;
    startClock = clock();        
    initMPC();
    readData(flag);

    selectTrainSet(40);

    if (1) {

        Init(1);


        train();

        endClock = clock();        
        double endtime = (double) (endClock - startClock) / CLOCKS_PER_SEC;
        long long onCom=Online();
        long long offCom=Offline();

        cout <<"Online comm: " <<onCom<<" MB"<<endl;
        cout <<"Offline comm: " <<offCom<<" MB"<<endl;

        cout << "Training time: " << endtime / 120.0+ onCom/(dataRate*60.0)<< " minutes" << endl;
        cout << "Offline Time: " << offCom/(dataRate*60.0)<< " minutes" << endl;

        writeModel(flag);

    }
    else{
        readModel(flag);
        Init(0);

        double accuracy = testData();

        endClock = clock();        

        cout << "---Accuracy：" << accuracy << "---" << endl;

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

