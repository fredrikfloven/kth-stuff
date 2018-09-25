#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

using namespace std;

int N;         //The amount of states
int M;         //The amount of emissions
int T;         //Amount of emissions observed

// Helper function for creating a matrix with elements from input
vector<vector<double> > createMatrix(int rows, int columns)
{
    vector<vector<double> > tempMatrix(rows);
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

// Helper function for creating a vector with elements from input
vector<double> createVector(int columns)
{
    vector<double> tempVector;
    for (int i = 0; i < columns; i++)
    {
        double element;
        cin >> element;
        tempVector.push_back(element);
    }
    return tempVector;
}

vector<double> currentState(vector<double> PI, vector<vector<double> > A)
{
    vector<double> currentStateMatrix;
    //Calculate probabilities of being in different states after 1 step
    for (int i = 0; i < N; i++)
    {
        double newValue = 0.0;
        for (int j = 0; j < N; j++)
            newValue += PI[j] * A[j][i];
        currentStateMatrix.push_back(newValue);
    }
    return currentStateMatrix;
}

vector<double> createEmissionVector(int columns, int rows, vector<double> currentStateMatrix, vector<vector<double> > B)
{
    vector<double> tempEmissionVector;
    //Calculate probabilities of being in different emissions after 1 step
    for (int i = 0; i < columns; i++)
    {
        double newValue = 0.0;
        for (int j = 0; j < rows; j++)
            newValue += currentStateMatrix[j] * B[j][i];
        tempEmissionVector.push_back(newValue);
    }
    return tempEmissionVector;
}

vector<vector<double> > createAlphaMatrix(vector<vector<double> > A, vector<vector<double> > B, vector<double> O, vector<double> PI)
{

    //ALPHA-PASS, FORWARD ALGORITHM
    vector<double> alphaCurrent;
    vector<double> alphaNext(N, 0);
    vector<vector<double> > alphaMatrix(T);

    //Initialize alpha0(i)
    for (int i = 0; i < N; i++)
    {
        alphaCurrent.push_back(PI[i] * B[i][O[0]]);
    }
    alphaMatrix[0] = alphaCurrent;

    //For observations alphat(i) t=1 to T
    for (int t = 1; t < T; t++)
    {
        //For alphat(i) i=0 to N
        for (int i = 0; i < N; i++)
        {
            double sum = 0;
            //For alphat(i) j=0 to N
            for (int j = 0; j < N; j++)
            {
                sum += (A[j][i] * alphaCurrent[j]);
            }
            sum *= B[i][O[t]];
            alphaNext[i] = sum;
        }
        alphaCurrent = alphaNext;
        alphaMatrix[t] = alphaCurrent;
    }

    return alphaMatrix;
}

vector<vector<double> > createBetaMatrix(vector<vector<double> > A, vector<vector<double> > B, vector<double> O)
{
    //BETA-PASS
    vector<vector<double> > betaMatrix(T);
    //beta initialization
    for (int i = 0; i < N; i++)
        betaMatrix[T - 1].push_back(1);

    //betat(i) t<T
    for (int t = T - 2; t >= 0; t--)
    {
        for (int i = 0; i < N; i++)
        {
            double beta = 0;
            for (int j = 0; j < N; j++)
            {
                beta += (A[i][j] * B[j][O[t + 1]] * betaMatrix[t + 1][j]);
            }
            betaMatrix[t].push_back(beta);
        }
    }

    return betaMatrix;
}

double createProbabilityOfSequence(int N, int T, vector < vector <double> > alphaMatrix){
    //Return probability for given emission sequence to happen
    double probabilityOfSequence;
    for (int i = 0; i < N; i++)
        probabilityOfSequence += alphaMatrix[T-1][i];
    return probabilityOfSequence;
}

vector < vector <double> > createNewA(vector <vector < vector <double> > > diGammaMatrix, vector < vector <double> > gammaMatrix){
    vector < vector <double> > newA(N, vector <double>(N, 0));
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                double nominator = 0;
                double denominator = 0;
                for(int t = 0; t < T-2; t++){
                    nominator += diGammaMatrix[t][i][j];
                    denominator += gammaMatrix[t][i];
                }
                newA[i][j] = nominator/denominator;
            }
        }
    return newA;
}

vector < vector <double> > createNewB(vector <double> O, vector < vector <double> > gammaMatrix){
    vector < vector <double> > newB(N, vector <double>(N, 0));
    for(int j = 0; j < N; j++){
        for(int k = 0; k < M; k++){
            double nominator = 0;
            double denominator = 0;
            for(int t = 0; t < T-1; t++){
                if(O[t] == k)
                    nominator += gammaMatrix[t][j];
                denominator += gammaMatrix[t][j];
            }
            newB[j][k] = nominator/denominator;
        }
    }
    return newB;
}

vector< vector < vector <double> > > createDiGammaMatrix (vector < vector <double> > alphaMatrix, vector < vector <double> > A, vector < vector <double> > B, vector < vector <double> > betaMatrix, vector <double> O, double probabilityOfSequence){
    vector< vector < vector <double> > > diGammaMatrix (T-1, vector < vector <double> >(N, vector <double> (N, 0)));
        for(int t = 0; t < T-1; t++){
            for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                    diGammaMatrix[t][i][j] = alphaMatrix[t][i] * A[i][j] * B[j][O[t+1]] * betaMatrix[t+1][j] / probabilityOfSequence;
                }
            }
        }
    return diGammaMatrix;
}

vector < vector <double> >  createGammaMatrix (vector < vector <double> > alphaMatrix, vector < vector <double> > betaMatrix, double probabilityOfSequence){
    vector < vector <double> >  gammaMatrix (T, vector <double>(N,0));
    //GAMMA MATRIX
    for(int t = 0; t < T; t++){
        for (int i = 0; i < N; i++){
            gammaMatrix[t][i] = alphaMatrix[t][i]*betaMatrix[t][i]/probabilityOfSequence;
        }
    }
    return gammaMatrix;
}

void printMatrix(vector<vector<double> > matrix)
{
    for (int i = 0; i < matrix.size(); i++)
    {
        double rowStochastic = 0.0;
        for (int j = 0; j < matrix[0].size(); j++)
        {
            rowStochastic += matrix[i][j];
            cout << matrix[i][j] << " ";
        }
        cout << "\n" << "Row sum: " << rowStochastic;
        cout << "\n";
    }
    cout << "\n";
}

void printVector(vector < double > v){
    for(int i = 0; i < v.size(); i++){
        cout << v[i] << " ";
    }
    cout << "\n";
}

int main()
{

    // SETUP START SETUP START SETUP START SETUP START SETUP START SETUP START 

    //TRANSITION MATRIX
    int N2; //The transition probabilities for each state
    cin >> N >> N2;
    vector<vector<double> > A;
    if (N == N2)
        A = createMatrix(N, N); //Create A

    //EMISSION MATRIX
    int M2; //The probabilities for each emission
    cin >> M >> M2;
    vector<vector<double> > B;
    if (M2 == N)
        B = createMatrix(N, M); //Create B

    //INITIAL STATE MATRIX
    int piRows; //Should only be 1
    int piColumns; // Should be N
    cin >> piRows >> piColumns;
    vector<double> PI;
    if (piRows == 1 && piColumns == N)
        PI = createVector(N); //Create pi

    //EMISSION SEQUENCE
    cin >> T;
    vector<double> O = createVector(T); //Create the emission sequence vector

    // SETUP END SETUP END SETUP END SETUP END SETUP END SETUP END SETUP END 


    // EMISSION PROBABILITY DISTRIBUTION
    vector<double> currentEmissionMatrix = createEmissionVector(N, M, PI, B);




    //ALPHA-PASS MATRIX TxN
    vector<vector<double> > alphaMatrix = createAlphaMatrix(A, B, O, PI);
    //printMatrix(alphaMatrix);

    //Return probability for given emission sequence to happen
    double probabilityOfSequence = createProbabilityOfSequence(N, T, alphaMatrix);
    cout << "Probability of given sequence: " << probabilityOfSequence << "\n";

    //BETA-PASS MATRIX TxN
    vector<vector<double> > betaMatrix = createBetaMatrix(A, B, O);
    //printMatrix(betaMatrix);

    //DI-GAMMA 3DMATRIX
    vector< vector < vector <double> > > diGammaMatrix = createDiGammaMatrix(alphaMatrix, A, B, betaMatrix, O, probabilityOfSequence);

    // for(int i = 0; i < T-1; i++){
    //     printMatrix(diGammaMatrix[i]);
    //     cout << "\n";
    // }

    //GAMMA MATRIX
    vector < vector <double> >  gammaMatrix = createGammaMatrix(alphaMatrix, betaMatrix, probabilityOfSequence);
    cout << "Gamma matrix" << "\n";
    //printMatrix(gammaMatrix);


    vector < vector <double> > newA(N, vector <double>(N, 0));
    vector < vector <double> > newB(N, vector <double>(M, 0));
    int iter = 0;
    int maxIter = 20;
    double previousPOS = probabilityOfSequence;
    while(iter < 20){

        // NEW INITIAL STATE MATRIX
        vector < double > newPI (T, 0);
        for(int i = 0; i < T-1; i++){
            newPI[i] = gammaMatrix[0][i];
        }
        // NEW TRANSITION MATRIX A
        newA = createNewA(diGammaMatrix, gammaMatrix);
        cout << "New A" << "\n";
        printMatrix(newA);
        // NEW OBSERVATION MATRIX B
        newB = createNewB(O, gammaMatrix);
        cout << "New B" << "\n";
        printMatrix(newB);

        alphaMatrix = createAlphaMatrix(newA, newB, O, newPI);
        probabilityOfSequence = createProbabilityOfSequence(N, T, alphaMatrix);
        betaMatrix = createBetaMatrix(newA, newB, O);
        diGammaMatrix = createDiGammaMatrix(alphaMatrix, newA, newB, betaMatrix, O, probabilityOfSequence);
        gammaMatrix = createGammaMatrix(alphaMatrix, betaMatrix, probabilityOfSequence);
        if(probabilityOfSequence <= previousPOS)
            break;
        iter++;
        previousPOS = probabilityOfSequence;
    }
    cout << "Iterations: " << iter << "\n";
    cout << "New A final" << "\n";
    printMatrix(newA);
    cout << "New B final" << "\n";
    printMatrix(newB);

    





    //VITERBI MATRIX
    vector<double> delta;
    vector<vector<double> > deltaMatrix(T);
    vector<vector<double> > backtrackMatrix(T - 1);

    //INITIALIZE DELTA
    for (int i = 0; i < N; i++)
    {
        delta.push_back(log(PI[i] * B[i][O[0]]));
    }
    if (T >= 1)
        deltaMatrix[0] = delta;

    //SUBSEQUENT DELTA
    for (int steps = 1; steps < T; steps++)
    {
        vector<double> deltaNext(N, 0);
        vector<double> deltaState(N, 0);
        for (int i = 0; i < N; i++)
        {

            double maxValue = -INFINITY;
            int index = 0;
            for (int j = 0; j < N; j++)
            {
                double tempValue = delta[j] + log(A[j][i]) + log(B[i][O[steps]]);
                if (tempValue > maxValue)
                {
                    maxValue = tempValue;
                    index = j;
                }
            }
            deltaNext[i] = maxValue;
            deltaState[i] = index;
        }
        delta = deltaNext;
        deltaMatrix[steps] = delta;
        backtrackMatrix[steps - 1] = deltaState;
    }

    vector<double> mostPossibleState;
    // CREATE DUMMY VECTOR WITH CORRECT LAST STATE
    for (int row = 0; row < T; row++)
    {
        double maxValue = -INFINITY;
        int index = 0;
        for (int column = 0; column < N; column++)
        {
            if (deltaMatrix[row][column] > maxValue)
            {
                maxValue = deltaMatrix[row][column];
                index = column;
            }
        }
        mostPossibleState.push_back(index);
    }
    //BACKTRACK FOR CORRECT STATES
    for (int t = 1; t < T; t++)
    {
        int row = mostPossibleState[T - t];
        mostPossibleState[T - t - 1] = backtrackMatrix[T - t - 1][row];
    }
    //PRINT MOST PROBABLE STATE SEQUENCE
    cout << "Most probable state sequence: " << "\n";
    printVector(mostPossibleState);
    // // PRINT DELTA MATRIX
    // for (int row = 0; row < T; row++){
    //     for (int column = 0; column < N; column++ ){
    //         cout << deltaMatrix[row][column] << " ";
    //     }
    //     cout << "\n";
    // }

    // cout << "\n" << "\n";

    // // PRINT BACKTRACK MATRIX I.E. DELTAIDX MATRIX
    // for (int row = 0; row < T-1; row++){
    //     for (int column = 0; column < N; column++ ){
    //         cout << backtrackMatrix[row][column] << " ";
    //     }
    //     cout << "\n";
    //}
}