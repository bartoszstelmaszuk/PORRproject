#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

int n = 2;

struct Point
{
    double *coordinates = (double *)malloc(n*sizeof(double));
}Point;

struct weightPoint
{
    vector<double> dimensionValues;
    vector<struct Point> blockCorners;
}weightPoint;

double einstainSummation(struct weightPoint value)
{
    double sum=0;

    for(int i=0; i<value.dimensionValues.size();i++) {
        sum+=value.dimensionValues[i];
    }

    return sum;
}

double einstainMultiplication(struct weightPoint value)
{
    double result=1;

    for (int i=1; i<=value.dimensionValues.size(); i++){
        result = result * cos(value.dimensionValues.at(i-1)/i);
    }

    return result;
}

double evaluationFunction(struct weightPoint value)
{
    return (1/40)*einstainSummation(value)+1-einstainMultiplication(value);
}

/*double approximationFunction(vector<double> value)
{
    return evaluationFunction(value) - (L/2)*sqrt();

}*/

struct weightPoint initFunctionDomain()
{
    struct weightPoint initialPoint;

    double *initValues = (double *)malloc(n*sizeof(double));
    for(int i=0; i<n; i++){
        initValues[i] = 0;
    }
    initialPoint.dimensionValues.assign(initValues, initValues+n);

    return initialPoint;

}
int main()
{
    struct weightPoint examplePoint;
    double myValues[] = {0,0};
    examplePoint.dimensionValues.assign(myValues, myValues+n);

    double result = einstainMultiplication(examplePoint);
    cout << "result multiplication: " << result << endl;

    result = einstainSummation(examplePoint);
    cout << "result summation: " << result << endl;

    result = evaluationFunction(examplePoint);
    cout << "result evaluationFunction: " << result << endl;


    return 0;
}
