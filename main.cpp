#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

int n = 2;
int devisionConstant = 2;

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

void devideBlock(struct weightPoint point)
{
    double longestDistance = 0;
    struct Point beginFinalPoint;
    struct Point endFinalPoint;

    for (int i=0; i<2*n; i++){
        struct Point beginPoint = point.blockCorners[i];
        struct Point endPoint = point.blockCorners[i+1];

        for (int j=0; j<2; j++)
        {
            double newDistance = beginPoint.coordinates[j] = endPoint.coordinates[j];
            if(newDistance>longestDistance) {
                longestDistance = newDistance;
                beginFinalPoint = beginPoint;
                endFinalPoint = endPoint;
            }
        }

    }

    struct weightPoint *newBlocks = (struct weightPoint*)malloc(devisionConstant*sizeof(struct weightPoint));
    struct weightPoint newPoint;
    newPoint.blockCorners = point.blockCorners;

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
