#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include "mpi.h"


using namespace std;

int n = 3;              // wymiar
double d = 0.020005; //punkt stopu
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
        sum+= ((i+1)*fabs(value.longestSide)*(value.longestSide));
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
    double middlePoint = 0;

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

    delete point1;
    delete point2;
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

void clearWeightPointVec(std::vector<weightPoint*>& vec){

    for(std::vector<weightPoint*>::iterator it = vec.begin(); it != vec.end(); ++it){
        delete (*it);
    }
    vec.clear();
}

void initWeightPoint(weightPoint* p, double (&tab)[8][3]){
    for(int i=0; i<8; ++i){
        Point *tmp = new Point;
        for(int j=0; j<3; ++j){
            tmp->coordinates.push_back(tab[i][j]);
        }
        p->blockCorners.push_back(*tmp);
    }
}
// for debug
void showTab(double (&tab)[8][3]){
    for(int i=0; i<8; ++i){
        cout << "[" << i << "]\t";
        for(int j=0; j<3; ++j){
            cout << tab[i][j] << "\t";
        }
        cout << endl;
    }
}

// DO WERYFIKACJI POPRAWNOSC OBLICZEN
// niech tab1 - tablica dobra dalej dzielona
void updateBounds(double (&tabToDevide)[8][3], double (&tab2)[8][3]){
    // wartosci min i max
    double min = tabToDevide[0][1], max = tabToDevide[2][1];
    double distance = fabs(min - max);
    // przepisz tablice dzielona do drugiej tablicy
    for(int i=0; i<8; ++i){
        for(int j=0; j<3; ++j){
            tab2[i][j] = tabToDevide[i][j];
        }
    }
    // zaktualizuj odpowiednie punkty obu tablic -> wspolrzedne y
    tabToDevide[2][1] -= distance/2;
    tabToDevide[3][1] -= distance/2;
    tabToDevide[6][1] -= distance/2;
    tabToDevide[7][1] -= distance/2;

    tab2[0][1] += distance/2;
    tab2[1][1] += distance/2;
    tab2[4][1] += distance/2;
    tab2[5][1] += distance/2;
}

int main(int argc, char* argv[])
{
    // zmienne dla MPI
    MPI_Status status;
    const int tagApprox = 4, tagSide = 5, tagGoOn = 6;
    const int RANK_MASTER = 0, RANK_SLAVE = 1;

    //clock_t startTime = clock();    // start time
    L = calculateLipschitzConstant()/2;
    cout << "L: " << L << endl;

    weightPoint *examplePoint1 = new weightPoint;
    weightPoint *examplePoint2 = new weightPoint;

    // tablice int przechowujace aktualne wspolrzedne blokow (wersja dla 3D) - 8 punktow po 3 wspolrzedne
    double tab1[8][3] = {
        {-40, -40, 40},
        {40, -40, 40},
        {40, 10, 40},
        {-40, 10, 40},
        {-40, -40, 10},
        {40, -40, 10},
        {40, 10, 10},
        {-40, 10, 10}
    };

    double tab2[8][3] = {
        {-40, 10, 10},
        {40, 10, 10},
        {40, 40, 10},
        {-40, 40, 10},
        {-40, 10, -40},
        {40, 10, -40},
        {40, 40, -40},
        {-40, 40, -40}
    };

    int rank = 0, numberProc = 0;
    double endCondition = 100.0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberProc);
    int kon = 0;
    while (kon < 3) {
        ++kon;
        //while (endCondition > d) {
        //showWeightPoint(*examplePoint1);
        //showWeightPoint(*examplePoint2);
        if(rank == 0){
            cout << "Process RANK = 0\n nic nie robieeee\n";
            // mieli dla pkt 1
            // czeka na wynik z drugiego procesu
            // porownuje ktory przedzial lepszy (gdzie approx jest mniejsze)
            // odpala metode update'ujaca przedzialy
            // i wszystko znowu sie mieli od nowa dla nowych tab1 i tab2
            //cout << "TAB1\n";
            //showTab(tab1);
            initWeightPoint(examplePoint1, tab1);
            calculateAttributes(examplePoint1);
            vector<weightPoint *> approxFuncArray;
            vector<weightPoint *> devidedBlocks;

            devidedBlocks = devideBlock(examplePoint1);

            weightPoint *block1 = new weightPoint();
            weightPoint *block2 = new weightPoint();
            weightPoint *chosenPoint = new weightPoint();
            copy_weightPoint(block1, devidedBlocks[0]);
            copy_weightPoint(block2, devidedBlocks[1]);
            calculateAttributes(block1);
            calculateAttributes(block2);
            approxFuncArray.push_back(block1);
            approxFuncArray.push_back(block2);
            copy_weightPoint(chosenPoint, block2);
            double minApproxFunction = 0.0;
            int nextToDevideBlockIndex = 0;
            double min = 5;
            int counter = 0;

            //while (chosenPoint->longestSide > d) {
            while(counter < 10){
                clearWeightPointVec(devidedBlocks);
                nextToDevideBlockIndex = 0;
                if(min > chosenPoint->longestSide){
                    min = chosenPoint->longestSide;
                }

                weightPoint *point = new weightPoint();
                copy_weightPoint(point, approxFuncArray[0]);
                minApproxFunction = point->approxFunctionValue;
                delete point;

                for(int i=0; i<approxFuncArray.size(); i++){
                    weightPoint *point = new weightPoint();
                    copy_weightPoint(point, approxFuncArray[i]);
                    // cout << "counter: " << i << " " << point->approxFunctionValue << endl;
                    if(point->approxFunctionValue < minApproxFunction) {
                        minApproxFunction = point->approxFunctionValue;
                        nextToDevideBlockIndex = i;
                    }
                    delete point;
                }
                substitute_weightPoint(chosenPoint, approxFuncArray[nextToDevideBlockIndex]);
                approxFuncArray.erase(approxFuncArray.begin()+nextToDevideBlockIndex);
                devidedBlocks = devideBlock(chosenPoint);
                calculateAttributes(devidedBlocks[0]);
                calculateAttributes(devidedBlocks[1]);
                weightPoint *block1 = new weightPoint();
                weightPoint *block2 = new weightPoint();
                copy_weightPoint(block1, devidedBlocks[0]);
                copy_weightPoint(block2, devidedBlocks[1]);
                approxFuncArray.push_back(block1);
                approxFuncArray.push_back(block2);
                counter++;
            }
            cout << "RANK = 0 DONE!\n ";
            double approxFunctionValueReceived = 0.0,
                    longestSideReceived = 0.0;
            // odbierz
            MPI_Recv(&approxFunctionValueReceived, 1, MPI_DOUBLE, RANK_SLAVE, tagApprox, MPI_COMM_WORLD, &status);
            MPI_Recv(&longestSideReceived, 1, MPI_DOUBLE, RANK_SLAVE, tagSide, MPI_COMM_WORLD, &status);
            cout << "ODEBRALEM: " << approxFunctionValueReceived << " " << longestSideReceived << endl;
            // i porownaj - najpierw sprawdz warunek zakonczenia - moze nie ma co dalej dzielic dziedziny

            if(approxFunctionValueReceived < chosenPoint->approxFunctionValue){   // czyli tab2 bylo lepsze
                updateBounds(tab2, tab1);
                endCondition = longestSideReceived;
                cout << "USTAWIAM END CONDITION1: " << endCondition;
            }
            else{
                updateBounds(tab1, tab2);
                endCondition = chosenPoint->longestSide;
                cout << "\tUSTAWIAM END CONDITION2: " << endCondition << endl;
            }
            // po updacie trzeba wyslac nowa tablice do drugiego procesu
            MPI_Send(&(tab2[0][0]), 24, MPI_DOUBLE, RANK_SLAVE, tagGoOn, MPI_COMM_WORLD );

            //showWeightPoint(*chosenPoint);
            delete block1;
            delete block2;

            // send signal that computations can be continued
            // MPI_Send(&CALC_ON, 1, MPI_INT, dest, tagGoOn, MPI_COMM_WORLD );

            /*delete chosenPoint;
            clearWeightPointVec(approxFuncArray);
            clearWeightPointVec(devidedBlocks);*/
        }

        else{
            cout << "TAB2\n";
            showTab(tab2);
            cout << "Process RANK = 1\n liczeeee\n";
            initWeightPoint(examplePoint2, tab2);
            calculateAttributes(examplePoint2);

            vector<weightPoint *> approxFuncArray;
            vector<weightPoint *> devidedBlocks;

            devidedBlocks = devideBlock(examplePoint2);

            weightPoint *block1 = new weightPoint();
            weightPoint *block2 = new weightPoint();
            weightPoint *chosenPoint = new weightPoint();
            copy_weightPoint(block1, devidedBlocks[0]);
            copy_weightPoint(block2, devidedBlocks[1]);
            calculateAttributes(block1);
            calculateAttributes(block2);
            approxFuncArray.push_back(block1);
            approxFuncArray.push_back(block2);

            copy_weightPoint(chosenPoint, block2);

            double minApproxFunction = 0.0;
            int nextToDevideBlockIndex = 0;
            double min = 5;
            int counter = 0;
            //while (chosenPoint->longestSide > d) {
            while(counter < 10){
                clearWeightPointVec(devidedBlocks);
                nextToDevideBlockIndex = 0;
                if(min > chosenPoint->longestSide){
                    min = chosenPoint->longestSide;
                }

                weightPoint *point = new weightPoint();
                copy_weightPoint(point, approxFuncArray[0]);
                minApproxFunction = point->approxFunctionValue;
                delete point;

                for(int i=0; i<approxFuncArray.size(); i++){
                    weightPoint *point = new weightPoint();
                    copy_weightPoint(point, approxFuncArray[i]);
                    // cout << "counter: " << i << " " << point->approxFunctionValue << endl;
                    if(point->approxFunctionValue < minApproxFunction) {
                        minApproxFunction = point->approxFunctionValue;
                        nextToDevideBlockIndex = i;
                    }
                    delete point;
                }

                substitute_weightPoint(chosenPoint, approxFuncArray[nextToDevideBlockIndex]);

                approxFuncArray.erase(approxFuncArray.begin()+nextToDevideBlockIndex);

                devidedBlocks = devideBlock(chosenPoint);

                calculateAttributes(devidedBlocks[0]);
                calculateAttributes(devidedBlocks[1]);

                weightPoint *block1 = new weightPoint();
                weightPoint *block2 = new weightPoint();
                copy_weightPoint(block1, devidedBlocks[0]);
                copy_weightPoint(block2, devidedBlocks[1]);
                approxFuncArray.push_back(block1);
                approxFuncArray.push_back(block2);

                counter++;
            }
            cout << "RANK = 1 DONE!\n ";
            double temp = chosenPoint->approxFunctionValue;
            double temp2 = chosenPoint->longestSide;
            cout << "WYSYLAM: " << temp << " " << temp2 << endl;
            MPI_Send(&temp, 1, MPI_DOUBLE, RANK_MASTER, tagApprox, MPI_COMM_WORLD);
            MPI_Send(&temp2, 1, MPI_DOUBLE, RANK_MASTER, tagSide, MPI_COMM_WORLD);
            //showWeightPoint(*chosenPoint);
            /*
            clock_t endTime = clock();
            double timeElapsed = (double) (endTime - startTime) / CLOCKS_PER_SEC;

            cout << "Time elapsed: " << timeElapsed << endl;*/
            /* delete point1;
        delete point2;
        delete point3;
        delete point4;
        delete point5;
        delete point6;
        delete point7;
        delete point8;*/
            /*delete point9;
        delete point10;
        delete point11;
        delete point12;
        delete point13;
        delete point14;
        delete point15;
        delete point16;*/
            //delete block1;
            //delete block2;
            delete chosenPoint;
            clearWeightPointVec(approxFuncArray);
            clearWeightPointVec(devidedBlocks);
            endCondition = chosenPoint->longestSide;
            // wait for new data to continue computations
            MPI_Recv(&tab2, 24, MPI_DOUBLE, RANK_MASTER, tagGoOn, MPI_COMM_WORLD, &status);

        }
    }

    MPI_Finalize();

    delete examplePoint1;
    delete examplePoint2;
    return 0;
}
