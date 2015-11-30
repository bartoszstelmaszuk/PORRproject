#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <string.h>

using namespace std;

int n = 2;
int devisionConstant = 2;

typedef struct Point
{
    double *coordinates;
}Point;

typedef struct weightPoint
{
    double approxFunctionValue;
    double longestSide;
    vector<double>dimensionsValue;
    vector<Point> blockCorners;
}weightPoint;

Point *newPoint ()
{
    Point *retPoint = (Point *)malloc(sizeof(Point));
    if (retPoint == NULL)
        return NULL;

    retPoint->coordinates = (double *)malloc(n*sizeof(double));
    if (retPoint->coordinates == NULL) {
        free (retPoint);
        return NULL;
    }

    return retPoint;
}

struct weightPoint *newWeightPoint ()
{
    weightPoint *retPoint = (weightPoint*)malloc(sizeof(weightPoint));
    if (retPoint == NULL)
        return NULL;

    retPoint->longestSide = 0;
    return retPoint;
}

weightPoint initFunctionDomain()
{
    weightPoint initialPoint;

    double *initValues = (double *)malloc(n*sizeof(double));
    for(int i=0; i<n; i++){
        initValues[i] = 0;
    }
    initialPoint.dimensionsValue.assign(initValues, initValues+n);

    return initialPoint;

}

void showBlockCorners(weightPoint weightPoint)
{
    struct Point newPoint;

    for(int i=0;i<weightPoint.blockCorners.size();i++){
        newPoint = weightPoint.blockCorners[i];
        cout << "weightPoint corner newPoint" << i<< " x: " << newPoint.coordinates[0] << endl;
        cout << "weightPoint corner newPoint" << i << " y: " << newPoint.coordinates[1] << endl;
    }
}

void showPoint(Point selectedPoint)
{
    cout << "Point x: " << selectedPoint.coordinates[0] << endl;
    cout << "Point y: " << selectedPoint.coordinates[1] << endl;
}

void copy_Point(Point *p1, Point *p2)
{
    for(int i=0; i<n; i++) {
        p1->coordinates[i] = p2->coordinates[i];
    }
}

void copy_weightPoint(weightPoint *p1, weightPoint *p2)
{
    for(int i=0; i < p2->blockCorners.size() ;i++) {
            Point *point = newPoint();
            copy_Point(point, &(p2->blockCorners[i]));
            p1->blockCorners.push_back(*point);
    }
    p1->longestSide = p2->longestSide;
}

double einstainSummationSquered(weightPoint value)
{
    double sum=0;

    for(int i=0; i<value.dimensionsValue.size();i++) {
        sum+= (value.dimensionsValue[i]*value.dimensionsValue[i]);
    }

    return sum;
}

double einstainSummation(weightPoint value)
{
    double sum=1;

    for(int i=0; i<n;i++) {
        sum+= ((i+1)*abs(value.longestSide)*(value.longestSide));
    }

    return sum;
}

double einstainMultiplication(weightPoint value)
{
    double result=1;

    for (int i=1; i<=value.dimensionsValue.size(); i++){
        result = result * cos(value.dimensionsValue.at(i-1)/i);
    }

    return result;
}

double evaluationFunction(weightPoint value)
{
    return (1/40)*einstainSummationSquered(value)+1-einstainMultiplication(value);
}

void showElementsOfApproxFunc(weightPoint point)
{
    double result = einstainMultiplication(point);
    cout << "result multiplication: " << result << endl;

    result = einstainSummationSquered(point);
    cout << "result summation: " << result << endl;

    result = evaluationFunction(point);
    cout << "result evaluationFunction: " << result << endl;
}

void showWeightPoint(weightPoint point)
{
    showBlockCorners(point);
    for(int i=0; i<n; i++) {
        cout << "Dimensional Value" << i <<": " << point.dimensionsValue[i] << endl;
    }
    cout << "Longest side: " << point.longestSide << endl;
    cout << "Approximation function: " << point.approxFunctionValue << endl;
    showElementsOfApproxFunc(point);
}

double calculateLongestSide(weightPoint *point)
{
    double longestDistance = 0;
    Point *beginPoint = newPoint();
    Point *endPoint = newPoint();

    for (int i=0; i<4; i++){

        copy_Point(beginPoint, &point->blockCorners[i]);
        copy_Point(endPoint, &point->blockCorners[(i+1)%4]);

        for (int j=0; j<2; j++)
        {
            double distance1 = beginPoint->coordinates[j];
            double distance2 = endPoint->coordinates[j];
            double newDistance;
            if(distance1 == distance2) {
                newDistance = 0;
            }else if ((distance1 <= 0 && distance2 >= 0) || (distance1 >= 0 || distance2 <= 0)){
                newDistance = fabs(distance1) + fabs(distance2);
            } else {
                newDistance = fabs(distance1 - distance2);
            }

            if(newDistance>longestDistance) {

                longestDistance = newDistance;
            }
        }
    }
    point->longestSide = longestDistance;

    return longestDistance;
}
void calculateDimensionValues(weightPoint *point)
{
    Point *beginPoint = newPoint();
    Point *endPoint = newPoint();
    double middlePoint;

    for (int i=0; i<n; i++){
        beginPoint = &point->blockCorners[i];
        endPoint = &point->blockCorners[(i+1)];

        for (int j=0; j<2; j++){
            double distance1 = beginPoint->coordinates[j];
            double distance2 = endPoint->coordinates[j];
            double newDistance;
            if(distance1 == distance2) {
                newDistance = 0;
            } else {
                newDistance = distance2 - distance1;
            }

            middlePoint = beginPoint->coordinates[j] + newDistance/2;

            if((i==j)){
                point->dimensionsValue.push_back(middlePoint);
            }
        }
    }
    free(beginPoint);
    free(endPoint);
}

void calculateApproximationFunction(weightPoint *point)
{
    point->approxFunctionValue = (evaluationFunction(*point) - (1.5)*sqrt(einstainSummation(*point)));
}

vector<weightPoint *> devideBlock(weightPoint *point)
{
    double longestDistance = 0;
    Point *firstDevisionPoint = newPoint();
    Point *secondDevisionPoint = newPoint();
    Point *beginPoint = newPoint();
    Point *endPoint = newPoint();
    Point *secondPoint = newPoint();
    int longestCornerIndex = 0;

    for (int i=0; i<4; i++){

        copy_Point(beginPoint, &point->blockCorners[i]);
        copy_Point(endPoint, &point->blockCorners[(i+1)%4]);
        copy_Point(secondPoint, &point->blockCorners[(i+2)%4]);

        for (int j=0; j<2; j++)
        {
            double distance1 = beginPoint->coordinates[j];
            double distance2 = endPoint->coordinates[j];
            double newDistance;
            if(distance1 == distance2) {
                newDistance = 0;
            }else if ((distance1 <= 0 && distance2 >= 0) || (distance1 >= 0 || distance2 <= 0)){
                newDistance = fabs(distance1) + fabs(distance2);
            } else {
                newDistance = fabs(distance1 - distance2);
            }

            if(newDistance > longestDistance) {

                longestDistance = newDistance;
                longestCornerIndex = i;
                copy_Point(firstDevisionPoint,endPoint);

                if(firstDevisionPoint->coordinates[j] >= 0) {
                    firstDevisionPoint->coordinates[j] = firstDevisionPoint->coordinates[j]-longestDistance/2;
                } else {
                    firstDevisionPoint->coordinates[j] = firstDevisionPoint->coordinates[j]+longestDistance/2;
                }

                copy_Point(secondDevisionPoint, secondPoint);

                if(secondDevisionPoint->coordinates[j] >=0) {
                    secondDevisionPoint->coordinates[j] = secondDevisionPoint->coordinates[j]-longestDistance/2;
                } else {
                    secondDevisionPoint->coordinates[j] = secondDevisionPoint->coordinates[j]+longestDistance/2;
                }
            }
        }

    }



    cout << "firstDevisionPoint " << endl;
    showPoint(*firstDevisionPoint);
    cout << "secondDevisionPoint " << endl;
    showPoint(*secondDevisionPoint);

    weightPoint *newPoint = newWeightPoint();
    weightPoint *newPoint2 = newWeightPoint();

    copy_weightPoint(newPoint, point);
    copy_weightPoint(newPoint2, point);

    if(longestCornerIndex ==0) {
        copy_Point(&(newPoint2->blockCorners[1]), firstDevisionPoint);
        copy_Point(&(newPoint2->blockCorners[2]), secondDevisionPoint);

        copy_Point(&(newPoint->blockCorners[0]), firstDevisionPoint);
        copy_Point(&(newPoint->blockCorners[3]), secondDevisionPoint);
    } else if (longestCornerIndex == 1) {
        copy_Point(&(newPoint->blockCorners[2]), firstDevisionPoint);
        copy_Point(&(newPoint->blockCorners[3]), secondDevisionPoint);

        copy_Point(&(newPoint2->blockCorners[1]), firstDevisionPoint);
        copy_Point(&(newPoint2->blockCorners[0]), secondDevisionPoint);
    }

    calculateDimensionValues(newPoint);
    calculateDimensionValues(newPoint2);
    calculateLongestSide(newPoint);
    calculateLongestSide(newPoint2);
    calculateApproximationFunction(newPoint);
    calculateApproximationFunction(newPoint2);

    vector<weightPoint *> newBlocks;
    newBlocks.push_back(newPoint);
    newBlocks.push_back(newPoint2);

    free(firstDevisionPoint);
    free(secondDevisionPoint);
    free(beginPoint);
    free(endPoint);
    free(secondPoint);

    return newBlocks;
}

int main()
{
    vector<struct Point> points;

    Point *point1 = newPoint();
    point1->coordinates[0] = -20;
    point1->coordinates[1] = -40;
    points.push_back(*point1);
    Point *point2 = newPoint();
    point2->coordinates[0] = 20;
    point2->coordinates[1] = -40;
    points.push_back(*point2);
    Point *point3 = newPoint();
    point3->coordinates[0] = 20;
    point3->coordinates[1] = 40;
    points.push_back(*point3);
    Point *point4 = newPoint();
    point4->coordinates[0] = -20;
    point4->coordinates[1] = 40;
    points.push_back(*point4);

    weightPoint *examplePoint = newWeightPoint();
    memcpy(&examplePoint->blockCorners, &points, sizeof(points));


    vector<weightPoint *> devidedBlocks;
    devidedBlocks = devideBlock(examplePoint);

    /*showWeightPoint(*(devidedBlocks[0]));
    cout << "Second block: " << endl;
    showWeightPoint(*(devidedBlocks[1]));

    devidedBlocks = devideBlock((devidedBlocks[0]));
    showWeightPoint(*(devidedBlocks[0]));

    showWeightPoint(*(devidedBlocks[1]));*/

    weightPoint *chosenPoint = newWeightPoint();
    double d = 0.0001; //punkt stopu

    while (chosenPoint.longestSide > d;) {

    }

    return 0;
}
