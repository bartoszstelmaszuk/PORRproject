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

    return retPoint;
}

double einstainSummation(weightPoint value)
{
    double sum=0;

    for(int i=0; i<value.dimensionsValue.size();i++) {
        sum+=value.dimensionsValue[i];
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
    return (1/40)*einstainSummation(value)+1-einstainMultiplication(value);
}

/*double approximationFunction(vector<double> value)
{
    return evaluationFunction(value) - (L/2)*sqrt();

}*/

weightPoint initFunctionDomain()
{
    weightPoint initialPoint;

    cout << "--------------" << endl;
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
}

vector<weightPoint *> devideBlock(weightPoint *point)
{
    double longestDistance = 0;
    struct Point *firstDevisionPoint = newPoint();
    struct Point *secondDevisionPoint = newPoint();
    struct Point *beginPoint = newPoint();
    struct Point *endPoint = newPoint();
    struct Point *secondPoint = newPoint();

    for (int i=0; i<3; i++){

        copy_Point(beginPoint, &point->blockCorners[i]);
        copy_Point(endPoint, &point->blockCorners[(i+1)%(2*n)]);
        copy_Point(secondPoint, &point->blockCorners[(i+2)%(2*n)]);

        for (int j=0; j<n; j++)
        {
            double distance1 = beginPoint->coordinates[j];
            double distance2 = endPoint->coordinates[j];
            double newDistance = abs(distance1) + abs(distance2);

            if(newDistance>longestDistance) {

                longestDistance = newDistance;
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

    copy_Point(&(newPoint->blockCorners[0]), firstDevisionPoint);
    copy_Point(&(newPoint->blockCorners[3]), secondDevisionPoint);

    copy_Point(&(newPoint2->blockCorners[1]), firstDevisionPoint);
    copy_Point(&(newPoint2->blockCorners[2]), secondDevisionPoint);

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
    weightPoint *examplePoint = newWeightPoint();
    vector<struct Point> points;

    Point *point1 = newPoint();
    point1->coordinates[0] = -40;
    point1->coordinates[1] = -40;
    points.push_back(*point1);
    Point *point2 = newPoint();
    point2->coordinates[0] = 40;
    point2->coordinates[1] = -40;
    points.push_back(*point2);
    Point *point3 = newPoint();
    point3->coordinates[0] = 40;
    point3->coordinates[1] = 40;
    points.push_back(*point3);
    Point *point4 = newPoint();
    point4->coordinates[0] = -40;
    point4->coordinates[1] = 40;
    points.push_back(*point4);

    memcpy(&examplePoint->blockCorners, &points, sizeof(points));

    /*double myValues[] = {0,0};
    examplePoint.dimensionsValue.assign(myValues, myValues+n);

    double result = einstainMultiplication(examplePoint);
    cout << "result multiplication: " << result << endl;

    result = einstainSummation(examplePoint);
    cout << "result summation: " << result << endl;

    result = evaluationFunction(examplePoint);
    cout << "result evaluationFunction: " << result << endl;*/

    vector<weightPoint *> devidedBlocks;
    devidedBlocks = devideBlock(examplePoint);

    showBlockCorners(*(devidedBlocks[0]));

    return 0;
}
