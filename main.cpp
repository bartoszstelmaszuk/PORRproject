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
    double approxFunctionValue;
    vector<double>dimensionsValue;
    vector<struct Point> blockCorners;
}weightPoint;

double einstainSummation(struct weightPoint value)
{
    double sum=0;

    for(int i=0; i<value.dimensionsValue.size();i++) {
        sum+=value.dimensionsValue[i];
    }

    return sum;
}

double einstainMultiplication(struct weightPoint value)
{
    double result=1;

    for (int i=1; i<=value.dimensionsValue.size(); i++){
        result = result * cos(value.dimensionsValue.at(i-1)/i);
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
    initialPoint.dimensionsValue.assign(initValues, initValues+n);

    return initialPoint;

}

struct weightPoint *devideBlock(struct weightPoint point)
{
    double longestDistance = 0;
    struct Point beginFinalPoint;
    struct Point endFinalPoint;

    for (int i=0; i<2*n-1; i++){
        struct Point beginPoint = point.blockCorners[i];
        struct Point endPoint = point.blockCorners[i+1];

        for (int j=0; j<n; j++)
        {
            double newDistance = abs(beginPoint.coordinates[j]) + abs(endPoint.coordinates[j]);

            if(newDistance>longestDistance) {
                longestDistance = newDistance;
                beginFinalPoint.coordinates = beginPoint.coordinates;
                endFinalPoint.coordinates = endPoint.coordinates;
                if(beginFinalPoint.coordinates[j] >=0) {
                    beginFinalPoint.coordinates[j] = beginPoint.coordinates[j]-longestDistance/2;
                } else {
                    beginFinalPoint.coordinates[j] = beginPoint.coordinates[j]+longestDistance/2;
                }
                if(endFinalPoint.coordinates[j] >=0) {
                    endFinalPoint.coordinates[j] = endPoint.coordinates[j]-longestDistance/2;
                } else {
                    endFinalPoint.coordinates[j] = endPoint.coordinates[j]+longestDistance/2;
                }
            }
        }

    }

    struct weightPoint *newBlocks = (struct weightPoint*)malloc(devisionConstant*sizeof(struct weightPoint));
    struct weightPoint newPoint;
    newPoint.blockCorners = point.blockCorners;
    newPoint.blockCorners[0] = beginFinalPoint;
    newPoint.blockCorners[1] = endFinalPoint;

    struct weightPoint newPoint2;
    newPoint2.blockCorners = point.blockCorners;
    newPoint2.blockCorners[2] = beginFinalPoint;
    newPoint2.blockCorners[3] = endFinalPoint;

    newBlocks[0] = newPoint;
    newBlocks[1] = newPoint2;

    return newBlocks;

}
int main()
{
    struct weightPoint examplePoint;
    vector<struct Point> points;

    struct Point point1;
    point1.coordinates[0] = -40;
    point1.coordinates[1] = -40;
    points.push_back(point1);
    struct Point point2;
    point2.coordinates[0] = 40;
    point2.coordinates[1] = -40;
    points.push_back(point2);
    struct Point point3;
    point3.coordinates[0] = 40;
    point3.coordinates[1] = 40;
    points.push_back(point3);
    struct Point point4;
    point4.coordinates[0] = -40;
    point4.coordinates[1] = 40;
    points.push_back(point4);

    examplePoint.blockCorners = points;

    double myValues[] = {0,0};
    examplePoint.dimensionsValue.assign(myValues, myValues+n);

    double result = einstainMultiplication(examplePoint);
    cout << "result multiplication: " << result << endl;

    result = einstainSummation(examplePoint);
    cout << "result summation: " << result << endl;

    result = evaluationFunction(examplePoint);
    cout << "result evaluationFunction: " << result << endl;

    struct weightPoint *devidedBlocks;
    devidedBlocks = devideBlock(examplePoint);
    struct Point newPoint;

    for(int i=0;i<4;i++){
        newPoint = devidedBlocks[0].blockCorners[i];
        cout << "newPoint" << i<< " x: " << newPoint.coordinates[0] << endl;
        cout << "newPoint" << i << " y: " << newPoint.coordinates[1] << endl;
    }

    return 0;
}
