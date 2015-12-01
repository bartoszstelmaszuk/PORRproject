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
    vector<double> coordinates;
}Point;

typedef struct weightPoint
{
    double approxFunctionValue;
    double longestSide;
    vector<double>dimensionsValue;
    vector<Point> blockCorners;
}weightPoint;

/*Point *newPoint ()
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
*/
/*weightPoint *newWeightPoint ()
{
    weightPoint *retPoint = (weightPoint*)malloc(10000*sizeof(weightPoint));
    if (retPoint == NULL)
        return NULL;

    retPoint->longestSide = 0;
    return retPoint;
}
*/
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
    Point newPoint;

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
        p1->coordinates.push_back(p2->coordinates[i]);
    }
}

void substitute_Point(Point *p1, Point *p2)
{
    for(int i=0; i<n; i++) {
        p1->coordinates[i] = p2->coordinates[i];
    }
}

void copy_weightPoint(weightPoint *p1, weightPoint *p2)
{
    for(int i=0; i < p2->blockCorners.size() ;i++) {
            Point *point = new Point;
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

void calculateApproximationFunction(weightPoint *point)
{
    point->approxFunctionValue = (evaluationFunction(*point) - (1.5)*sqrt(einstainSummation(*point)));
}

void calculateLongestSide(weightPoint *point)
{
    double longestDistance = 0;
    Point *beginPoint = new Point;
    Point *endPoint = new Point;

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

    delete(beginPoint);
    delete(endPoint);

    calculateApproximationFunction(point);
}
void calculateDimensionValues(weightPoint *point)
{
    Point *beginPoint = new Point;
    Point *endPoint = new Point;
    double middlePoint;

    for (int i=0; i<n; i++){
        copy_Point(beginPoint, &point->blockCorners[i]);
        copy_Point(endPoint, &point->blockCorners[i+1]);

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

    delete(beginPoint);
    delete(endPoint);

    calculateLongestSide(point);
}

void calculateAttributes(weightPoint *point)
{
    calculateDimensionValues(point);
}

vector<weightPoint *> devideBlock(weightPoint *point)
{
    double longestDistance = 0;
    Point *firstDevisionPoint = new Point;
    Point *secondDevisionPoint = new Point;
    Point *beginPoint = new Point;
    Point *endPoint = new Point;
    Point *secondPoint = new Point;
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

    weightPoint *newPoint = new weightPoint;
    weightPoint *newPoint2 = new weightPoint;

    copy_weightPoint(newPoint, point);
    copy_weightPoint(newPoint2, point);

    if(longestCornerIndex ==0) {
        substitute_Point(&(newPoint2->blockCorners[1]), firstDevisionPoint);
        substitute_Point(&(newPoint2->blockCorners[2]), secondDevisionPoint);

        substitute_Point(&(newPoint->blockCorners[0]), firstDevisionPoint);
        substitute_Point(&(newPoint->blockCorners[3]), secondDevisionPoint);
    } else if (longestCornerIndex == 1) {
        substitute_Point(&(newPoint->blockCorners[2]), firstDevisionPoint);
        substitute_Point(&(newPoint->blockCorners[3]), secondDevisionPoint);

        substitute_Point(&(newPoint2->blockCorners[1]), firstDevisionPoint);
        substitute_Point(&(newPoint2->blockCorners[0]), secondDevisionPoint);
    }

    calculateAttributes(newPoint);
    calculateAttributes(newPoint2);

    vector<weightPoint *> newBlocks;
    newBlocks.push_back(newPoint);
    newBlocks.push_back(newPoint2);

    delete(firstDevisionPoint);
    delete(secondDevisionPoint);
    delete(beginPoint);
    delete(endPoint);
    delete(secondPoint);

    return newBlocks;
}

int main()
{
    vector<struct Point> points;

    weightPoint *examplePoint = new weightPoint;
    Point *point1 = new Point;
    point1->coordinates.push_back(-20);
    point1->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point1);
    Point *point2 = new Point;
    point2->coordinates.push_back(20);
    point2->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point2);
    Point *point3 = new Point;
    point3->coordinates.push_back(20);
    point3->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point3);
    Point *point4 = new Point;
    point4->coordinates.push_back(-20);
    point4->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point4);

    calculateAttributes(examplePoint);

    showWeightPoint(*examplePoint);

    weightPoint *chosenPoint;
    vector<weightPoint *> approxFuncArray;

    chosenPoint = examplePoint;

    //showWeightPoint(*examplePoint);

    vector<weightPoint *> devidedBlocks;
    devidedBlocks = devideBlock(chosenPoint);

    calculateAttributes(devidedBlocks[0]);
    calculateAttributes(devidedBlocks[1]);

    //approxFuncArray.push_back(devidedBlocks[0]);
    //approxFuncArray.push_back(devidedBlocks[1]);

    showWeightPoint(*(devidedBlocks[0]));
    showWeightPoint(*(devidedBlocks[1]));

    //double minApproxFunction = chosenPoint->approxFunctionValue;

    double d = 1; //punkt stopu

    //devidedBlocks = devideBlock(devidedBlocks[0]);

    //showWeightPoint(*(devidedBlocks[0]));

   /* while (chosenPoint->longestSide > d) {
        int nextToDevideBlockIndex = 0;
        weightPoint *point;

        for(int i=0; i<approxFuncArray.size(); i++){
            point = approxFuncArray[i];
            if(point->approxFunctionValue < minApproxFunction) {
                nextToDevideBlockIndex = i;
            }
        }
        copy_weightPoint(chosenPoint, approxFuncArray[nextToDevideBlockIndex]);

        approxFuncArray.erase(approxFuncArray.begin()+nextToDevideBlockIndex);

        devidedBlocks = devideBlock(chosenPoint);

        calculateAttributes(devidedBlocks[0]);
        calculateAttributes(devidedBlocks[1]);

        approxFuncArray.push_back(devidedBlocks[0]);
        approxFuncArray.push_back(devidedBlocks[1]);
    }

    showWeightPoint(*chosenPoint);*/

    return 0;
}
