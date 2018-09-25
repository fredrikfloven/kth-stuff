#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm> 
#include <math.h>

using namespace std;

int transitionRows; //The amount of states
int transitionColumns; //The transition probabilities for each state
int emissionRows; //The amount of emissions
int emissionColumns; //The probabilities for each emission
int initialStateRows; //Should only be 1
int initialStateColumns; // The initial probability for which state we start in

vector < vector < double > > createMatrix(int rows, int columns) {
    vector < vector <double> > tempMatrix(rows);
    for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                double element;
                cin >> element;
                tempMatrix[i].push_back(element);
            }
        }
    return tempMatrix;
}

vector < double > createVector(int columns) {
    vector < double > tempVector;
    for (int i = 0; i < columns; i++)
        {
            double element;
            cin >> element;
            tempVector.push_back(element);
        }
    return tempVector;
}

vector < double > currentState (vector < double > initialStateMatrix, vector < vector < double > > transitionMatrix)
{
    vector < double > currentStateMatrix;
    //Calculate probabilities of being in different states after 1 step
    for(int i = 0; i < transitionColumns; i++){
        double newValue = 0.0;
        for(int j = 0; j < transitionRows; j++)
            newValue += initialStateMatrix[j]*transitionMatrix[j][i];
        currentStateMatrix.push_back(newValue);
    }
    return currentStateMatrix;
}

vector < double > createEmissionVector(int columns, int rows, vector < double > currentStateMatrix, vector < vector <double> > emissionMatrix){
    vector < double > tempEmissionVector;
    //Calculate probabilities of being in different emissions after 1 step
    for(int i = 0; i < columns; i++){
        double newValue = 0.0;
        for(int j = 0; j < rows; j++)
            newValue += currentStateMatrix[j]*emissionMatrix[j][i];
        tempEmissionVector.push_back(newValue);
    }
    return tempEmissionVector;
    
}


int main()
{
    //TRANSITION MATRIX
    cin >> transitionRows >> transitionColumns;
    vector < vector <double> > transitionMatrix = createMatrix(transitionRows, transitionColumns); //Create A

    //EMISSION MATRIX
    cin >> emissionRows >> emissionColumns;
    vector < vector <double> > emissionMatrix = createMatrix(emissionRows, emissionColumns); //Create B

    //INITIAL STATE MATRIX
    cin >> initialStateRows >> initialStateColumns;
    vector < double > initialStateMatrix = createVector(initialStateColumns); //Create pi


    //EMISSION SEQUENCE
    int emissionLengthOfSequence;
    cin >> emissionLengthOfSequence;
    vector < double > emissionVector = createVector(emissionLengthOfSequence); //Create the emission sequence vector


    //CALCULATIONS
    vector < double > currentStateMatrix = initialStateMatrix;
    vector < double > currentEmissionMatrix;
    currentEmissionMatrix = createEmissionVector(emissionColumns, emissionRows, currentStateMatrix, emissionMatrix);

    
    //VITERBI
    vector < double > delta;
    vector < vector <double> > deltaMatrix(emissionLengthOfSequence);
    vector < vector <double> > backtrackMatrix(emissionLengthOfSequence - 1);

    //INITIALIZE DELTA
    for (int i = 0; i < transitionRows ; i++){
        delta.push_back(log(initialStateMatrix[i]*emissionMatrix[i][emissionVector[0]]));
    }
    if(emissionLengthOfSequence >= 1)
        deltaMatrix[0] = delta;
 
    //SUBSEQUENT DELTA
    for(int steps = 1; steps < emissionLengthOfSequence; steps++){
        vector < double > deltaNext(transitionRows,0);
        vector < double > deltaState(transitionRows,0);
        for(int i = 0; i < transitionRows; i++){

            double maxValue = -INFINITY;
            int index = 0;
            for(int j = 0; j < transitionRows; j++){
                double tempValue = delta[j] + log(transitionMatrix[j][i]) + log(emissionMatrix[i][emissionVector[steps]]);
                if(tempValue > maxValue){
                    maxValue = tempValue;
                    index = j;
                }
            }
            deltaNext[i] = maxValue;
            deltaState[i] = index;
        }
        deltaMatrix[steps] = deltaNext;
        backtrackMatrix[steps-1] = deltaState;
        delta = deltaNext;
    }


    vector < double > mostPossibleState;
    // CREATE DUMMY VECTOR WITH CORRECT LAST STATE
    for (int row = 0; row < emissionLengthOfSequence; row++){
        double maxValue = -INFINITY;
        int index = 0;
        for (int column = 0; column < transitionRows; column++ ){
            if (deltaMatrix[row][column] > maxValue){
                maxValue = deltaMatrix[row][column];
                index = column;
            }
        }
        mostPossibleState.push_back(index);
    }    
    //BACKTRACK FOR CORRECT STATES
    for(int t = 1; t < emissionLengthOfSequence ; t++){
        int row = mostPossibleState[emissionLengthOfSequence-t];
        mostPossibleState[emissionLengthOfSequence-t-1] = backtrackMatrix[emissionLengthOfSequence-t-1][row];
    }
    //PRINT MOST PROBABLE STATE SEQUENCE
    for(int t = 0; t < emissionLengthOfSequence; t++){
        cout << mostPossibleState[t] << " ";
    }
    // // PRINT DELTA MATRIX 
    // for (int row = 0; row < emissionLengthOfSequence; row++){
    //     for (int column = 0; column < transitionRows; column++ ){
    //         cout << deltaMatrix[row][column] << " ";
    //     }
    //     cout << "\n";
    // }

    // cout << "\n" << "\n";

    // // PRINT BACKTRACK MATRIX I.E. DELTAIDX MATRIX 
    // for (int row = 0; row < emissionLengthOfSequence-1; row++){
    //     for (int column = 0; column < transitionRows; column++ ){
    //         cout << backtrackMatrix[row][column] << " ";
    //     }
    //     cout << "\n";
    //}

}