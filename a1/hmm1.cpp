#include <cstdio>
#include <vector>
#include <iostream>

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

    
    //ALPHA-PASS
    vector < double > alphaCurrent;
    vector < double > alphaNext (transitionRows,0);

    for (int i = 0; i < transitionRows ; i++){
        alphaCurrent.push_back(initialStateMatrix[i]*emissionMatrix[i][emissionVector[0]]);
    }
 

    for(int k = 1; k < emissionLengthOfSequence; k++){
        for (int i = 0; i < transitionRows ; i++){
            double sum = 0;
            for(int j = 0; j < transitionColumns ; j++){
                sum += (transitionMatrix[j][i] * alphaCurrent[j]);
            }
            sum *= emissionMatrix[i][emissionVector[k]];
            alphaNext[i] = sum;
        }
        alphaCurrent = alphaNext;
    }
    double sum;
    for(int i = 0; i < alphaNext.size(); i++)
        sum += alphaNext[i];
    cout << sum;
    
    // for (int row = 0; row < 4; row++){
    //     for (int column = 0; column < 4; column++ )
    //         cout << transitionMatrix[row][column] << " ";
    //     cout << "\n";
    // }

}