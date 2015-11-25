#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

vector<double> xCoordinates;

double einstainSummation(int n)
{
    double sum=0;

    for(int i=0; i<n;i++) {
        sum+=xCoordinates[i];
    }

    return sum;
}

double einstainMultiplication(int n)
{
    double result=1;

    for (int i=1; i<=n; i++){
        result = result * cos(xCoordinates.at(i-1)/i);
    }

    return result;
}
int main()
{
    double myValues[] = {1,2,3,4};
    xCoordinates.assign(myValues, myValues+4);
    double result = einstainMultiplication(2);
    cout << "result multiplication: " << result << endl;
    result = einstainSummation(2);
    cout << "result summation: " << result << endl;


    return 0;
}
