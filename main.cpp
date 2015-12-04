#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <string.h>

using namespace std;

int n = 2;
int devisionConstant = 2;
double L;

int lowerXBounder = -40;
int higherXBounder = 40;

int lowerYBounder = -40;
int higherYBounder = 40;

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

void showBlockCorners(weightPoint weightPoint)
{
    Point newPoint;

    for(int i=0;i<4;i++){
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
    for(int i=0; i<p2->coordinates.size(); i++) {
        p1->coordinates.push_back(p2->coordinates[i]);
    }
}

void substitute_Point(Point *p1, Point *p2)
{
    for(int i=0; i<p2->coordinates.size(); i++) {
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
    p1->approxFunctionValue = p2->approxFunctionValue;
    for(int i=0; i < p2->dimensionsValue.size() ;i++) {
            p1->dimensionsValue.push_back(p2->dimensionsValue[i]);
    }
}

void substitute_weightPoint(weightPoint *p1, weightPoint *p2)
{
    p1->blockCorners.clear();
    p1->dimensionsValue.clear();
    for(int i=0; i < p2->blockCorners.size() ;i++) {
            Point *point = new Point;
            copy_Point(point, &(p2->blockCorners[i]));
            p1->blockCorners.push_back(*point);
    }
    p1->longestSide = p2->longestSide;
    p1->approxFunctionValue = p2->approxFunctionValue;
    for(int i=0; i < p2->dimensionsValue.size() ;i++) {
            p1->dimensionsValue.push_back(p2->dimensionsValue[i]);
    }
}

double einstainSummationSquered(weightPoint value)
{
    double sum=0;

    for(int i=0; i<value.dimensionsValue.size();i++) {
        sum+= (value.dimensionsValue[i]*value.dimensionsValue[i]);
    }

    return sum;
}

double einstainSummationSquered(Point value)
{
    double sum=0;

    for(int i=0; i<value.coordinates.size();i++) {
        sum+= (value.coordinates[i]*value.coordinates[i]);
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

double einstainMultiplication(Point value)
{
    double result=1;

    for (int i=1; i<=value.coordinates.size(); i++){
        result = result * cos(value.coordinates.at(i-1)/i);
    }

    return result;
}

double evaluationFunction(weightPoint value)
{
    return (0.025)*einstainSummationSquered(value)+1-einstainMultiplication(value);
}

double evaluationFunction(Point value)
{
    return (0.025)*einstainSummationSquered(value)+1-einstainMultiplication(value);
}

void showElementsOfApproxFunc(weightPoint point)
{
    /*double result = einstainMultiplication(point);
    cout << "result multiplication: " << result << endl;

    result = einstainSummationSquered(point);
    cout << "result summation: " << result << endl;*/

    double result = evaluationFunction(point);
    cout << "result evaluationFunction: " << result << endl;
}

void showWeightPoint(weightPoint point)
{
    showBlockCorners(point);
    for(int i=0; i<point.dimensionsValue.size(); i++) {
        cout << "Dimensional Value" << i <<": " << point.dimensionsValue[i] << endl;
    }
    cout << "Longest side: " << point.longestSide << endl;
    cout << "Approximation function: " << point.approxFunctionValue << endl;
    showElementsOfApproxFunc(point);
}

void calculateApproximationFunction(weightPoint *point)
{
    point->approxFunctionValue = (evaluationFunction(*point) - (L)*sqrt(einstainSummation(*point)));
}

void calculateLongestSide(weightPoint *point)
{
    double longestDistance = 0;
    Point *beginPoint = new Point;
    Point *endPoint = new Point;

    copy_Point(beginPoint, &point->blockCorners[0]);
    copy_Point(endPoint, &point->blockCorners[1]);

    for (int i=0; i<2; i++){

        substitute_Point(beginPoint, &point->blockCorners[i]);
        substitute_Point(endPoint, &point->blockCorners[(i+1)%4]);

        for (int j=0; j<2; j++)
        {
            double distance1 = beginPoint->coordinates[j];
            double distance2 = endPoint->coordinates[j];
            double newDistance;
            if(distance1 == distance2) {
                newDistance = 0;
            }else if ((distance1 < 0 && distance2 > 0) || (distance1 > 0 && distance2 < 0)){
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

    point->dimensionsValue.clear();

    copy_Point(beginPoint, &point->blockCorners[0]);
    copy_Point(endPoint, &point->blockCorners[1]);

    for (int i=0; i<n; i++){
        substitute_Point(beginPoint, &point->blockCorners[i]);
        substitute_Point(endPoint, &point->blockCorners[i+1]);

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

double calculateDistanceBetweenPoints(Point p1, Point p2)
{
    double answer = 0;
    vector<double> result;

    for(int i=0; i<n; i++) {
        result.push_back(p1.coordinates[i]-p2.coordinates[i]);
    }

    for(int i=0; i<n; i++) {
        answer += result.at(i)*result.at(i);
    }

    answer = sqrt(answer);

    return answer;
}

double calculateLipschitzConstant()
{
    double answer;

    double r = ((double) rand() / (RAND_MAX));

    Point *point1 = new Point;
    point1->coordinates.push_back(lowerXBounder+(higherXBounder-lowerXBounder)*r);
    r = ((double) rand() / (RAND_MAX));
    point1->coordinates.push_back(lowerYBounder+(higherYBounder-lowerYBounder)*r);

    Point *point2 = new Point;
    r = ((double) rand() / (RAND_MAX));
    point2->coordinates.push_back(lowerXBounder+(higherXBounder-lowerXBounder)*r);
    r = ((double) rand() / (RAND_MAX));
    point2->coordinates.push_back(lowerYBounder+(higherYBounder-lowerYBounder)*r);

    double functionDiff = evaluationFunction(*point1)-evaluationFunction(*point2);
    double argDiff = calculateDistanceBetweenPoints(*point1,*point2);

    answer = fabs(functionDiff/argDiff);

    return answer;
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

    copy_Point(beginPoint, &point->blockCorners[0]); //copy point to declare sizes of vector
    copy_Point(endPoint, &point->blockCorners[1]);  // after that we use substitute_point to change values of points
    copy_Point(secondPoint, &point->blockCorners[2]);

    copy_Point(firstDevisionPoint,endPoint);
    copy_Point(secondDevisionPoint, secondPoint);

    for (int i=0; i<4; i++){

        substitute_Point(beginPoint, &point->blockCorners[i]);
        substitute_Point(endPoint, &point->blockCorners[(i+1)%4]);
        substitute_Point(secondPoint, &point->blockCorners[(i+2)%4]);

        for (int j=0; j<2; j++)
        {
            double distance1 = beginPoint->coordinates[j];
            double distance2 = endPoint->coordinates[j];
            double newDistance;
            if(distance1 == distance2) {
                newDistance = 0;
            }else if ((distance1 < 0 && distance2 > 0) || (distance1 > 0 && distance2 < 0)){
                newDistance = fabs(distance1) + fabs(distance2);
            } else {
                newDistance = fabs(distance1 - distance2);
            }

            if(newDistance > longestDistance) {

                longestDistance = newDistance;
                longestCornerIndex = i;
                substitute_Point(firstDevisionPoint,endPoint);

                if(firstDevisionPoint->coordinates[j] >= 0) {
                    firstDevisionPoint->coordinates[j] = firstDevisionPoint->coordinates[j]-longestDistance/2;
                } else {
                    firstDevisionPoint->coordinates[j] = firstDevisionPoint->coordinates[j]+longestDistance/2;
                }

                substitute_Point(secondDevisionPoint, secondPoint);

                if(secondDevisionPoint->coordinates[j] >=0) {
                    secondDevisionPoint->coordinates[j] = secondDevisionPoint->coordinates[j]-longestDistance/2;
                } else {
                    secondDevisionPoint->coordinates[j] = secondDevisionPoint->coordinates[j]+longestDistance/2;
                }
            }
        }

    }

    /*cout << "firstDevisionPoint " << endl;
    showPoint(*firstDevisionPoint);
    cout << "secondDevisionPoint " << endl;
    showPoint(*secondDevisionPoint);*/

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
    L = calculateLipschitzConstant()/2;

    cout << "L: " << L << endl;

    vector<struct Point> points;

    weightPoint *examplePoint = new weightPoint;
    /*
    //2D
    Point *point1 = new Point;
    point1->coordinates.push_back(-40);
    point1->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point1);
    Point *point2 = new Point;
    point2->coordinates.push_back(40);
    point2->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point2);
    Point *point3 = new Point;
    point3->coordinates.push_back(40);
    point3->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point3);
    Point *point4 = new Point;
    point4->coordinates.push_back(-40);
    point4->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point4);*/

    Point *point1 = new Point;
    point1->coordinates.push_back(-40);
    point1->coordinates.push_back(-40);
    point1->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point1);
    Point *point2 = new Point;
    point2->coordinates.push_back(40);
    point2->coordinates.push_back(-40);
    point2->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point2);
    Point *point3 = new Point;
    point3->coordinates.push_back(40);
    point3->coordinates.push_back(40);
    point3->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point3);
    Point *point4 = new Point;
    point4->coordinates.push_back(-40);
    point4->coordinates.push_back(40);
    point4->coordinates.push_back(40);
    examplePoint->blockCorners.push_back(*point4);
    Point *point5 = new Point;
    point5->coordinates.push_back(-40);
    point5->coordinates.push_back(-40);
    point5->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point5);
    Point *point6 = new Point;
    point6->coordinates.push_back(40);
    point6->coordinates.push_back(-40);
    point6->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point6);
    Point *point7 = new Point;
    point7->coordinates.push_back(40);
    point7->coordinates.push_back(40);
    point7->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point7);
    Point *point8 = new Point;
    point8->coordinates.push_back(-40);
    point8->coordinates.push_back(40);
    point8->coordinates.push_back(-40);
    examplePoint->blockCorners.push_back(*point8);

    calculateAttributes(examplePoint);

    showWeightPoint(*examplePoint);

    vector<weightPoint *> approxFuncArray;

    vector<weightPoint *> devidedBlocks;
    devidedBlocks = devideBlock(examplePoint);

    weightPoint *block1 = new weightPoint();
    copy_weightPoint(block1, devidedBlocks[0]);

    weightPoint *block2 = new weightPoint();
    copy_weightPoint(block2, devidedBlocks[1]);

    calculateAttributes(block1);
    calculateAttributes(block2);

    approxFuncArray.push_back(block1);
    approxFuncArray.push_back(block2);

    //showWeightPoint(*(approxFuncArray[0]));
    //showWeightPoint(*(approxFuncArray[1]));

    weightPoint *chosenPoint = new weightPoint();
    copy_weightPoint(chosenPoint, block2);

    double minApproxFunction;

    double d = 0.1; //punkt stopu
    int nextToDevideBlockIndex = 0;

    int k=0;

    while (chosenPoint->longestSide > d) {
        nextToDevideBlockIndex = 0;
        cout << k << " " << chosenPoint->longestSide << endl;
        k++;

        weightPoint *point = new weightPoint();
        copy_weightPoint(point, approxFuncArray[0]);
        minApproxFunction = point->approxFunctionValue;
        delete point;

        for(int i=0; i<approxFuncArray.size(); i++){
            weightPoint *point = new weightPoint();
            copy_weightPoint(point, approxFuncArray[i]);
            //cout << "counter: " << i << " " << point->approxFunctionValue << endl;
            if(point->approxFunctionValue < minApproxFunction) {
                minApproxFunction = point->approxFunctionValue;
                nextToDevideBlockIndex = i;
            }
            delete point;
        }

        substitute_weightPoint(chosenPoint, approxFuncArray[nextToDevideBlockIndex]);

        approxFuncArray.erase(approxFuncArray.begin()+nextToDevideBlockIndex);

        //showWeightPoint(*chosenPoint);
        devidedBlocks = devideBlock(chosenPoint);

        calculateAttributes(devidedBlocks[0]);
        calculateAttributes(devidedBlocks[1]);

        weightPoint *block1 = new weightPoint();
        copy_weightPoint(block1, devidedBlocks[0]);

        weightPoint *block2 = new weightPoint();
        copy_weightPoint(block2, devidedBlocks[1]);

        approxFuncArray.push_back(block1);
        approxFuncArray.push_back(block2);

        /*for (int l=0; l<approxFuncArray.size(); l++){
            showWeightPoint(*(approxFuncArray[l]));
        }
        showWeightPoint(*(devidedBlocks[0]));
        showWeightPoint(*(devidedBlocks[1]));*/
    }

    showWeightPoint(*chosenPoint);

    return 0;
}
