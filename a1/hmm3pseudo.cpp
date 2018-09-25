#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

using namespace std;

int N = 0;         //The amount of states
int M = 0;         //The amount of emissions
int T = 0;         //Amount of emissions observed

// Helper function for creating a matrix with elements from input
vector<vector<double> > createMatrix(int rows, int columns)
{
    vector<vector<double> > tempMatrix(rows, vector <double>(columns, 0));
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            double element;
            cin >> element;
            tempMatrix[i][j] = element;
        }
    }
    return tempMatrix;
}

// Helper function for creating a vector with elements from input
vector<double> createVector(int columns)
{
    vector<double> tempVector(columns, 0);
    for (int i = 0; i < columns; i++)
    {
        double element;
        cin >> element;
        tempVector[i] = element;
    }
    return tempVector;
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

void printMatrixAsRows(vector<vector<double> > matrix)
{
    for (int i = 0; i < matrix.size(); i++)
    {
        double rowStochastic = 0.0;
        for (int j = 0; j < matrix[0].size(); j++)
        {
            rowStochastic += matrix[i][j];
            cout << matrix[i][j] << " ";
        }
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
    int N2 = 0; //The transition probabilities for each state
    cin >> N >> N2;
    vector<vector<double> > A;
    if (N == N2)
        A = createMatrix(N, N); //Create A

    //EMISSION MATRIX
    int M2 = 0; //The probabilities for each emission
    cin >> M2 >> M;
    vector<vector<double> > B;
    if (M2 == N)
        B = createMatrix(N, M); //Create B

    //INITIAL STATE MATRIX
    int piRows = 0; //Should only be 1
    int piColumns = 0; // Should be N
    cin >> piRows >> piColumns;
    vector<double> PI;
    if (piRows == 1 && piColumns == N)
        PI = createVector(N); //Create pi

    //EMISSION SEQUENCE
    cin >> T;
    vector<double> O = createVector(T); //Create the emission sequence vector

    // SETUP END SETUP END SETUP END SETUP END SETUP END SETUP END SETUP END 



    // step 1, initialize
    int maxIters = 50;
    int iters = 0;
    double oldLogProb = -INFINITY;
    while(iters<maxIters){
        // Step 2, alpha-pass

        // compute alpha0(i)
        vector <double> c(T, 0);
        vector<vector<double> > alphaMatrix(T, vector<double>(N, 0));
        for(int i = 0; i < N; i++){
            alphaMatrix[0][i] = PI[i]*B[i][O[0]];
            c[0] += alphaMatrix[0][i];
        }
        // scale the a0(i)
        c[0] = 1/c[0];
        for(int i = 0; i < N; i++)
            alphaMatrix[0][i] *= c[0];

        // compute alphat(i)
        for(int t = 1; t < T; t++){
            for(int i = 0; i < N; i++){
                for(int j = 0; j < N; j++)
                    alphaMatrix[t][i] += alphaMatrix[t-1][j] * A[j][i];
                alphaMatrix[t][i] *= B[i][O[t]];
                c[t] += alphaMatrix[t][i];
            }
            // scale alphat(i)
            c[t] = 1/c[t];
            for(int i = 0; i < N; i++)
                alphaMatrix[t][i] *= c[t];
        }


        // Step 3, The beta-pass

        vector<vector<double> > betaMatrix(T, vector<double>(N, 0));
        // Let betaT-1(i) = 1, scaled by cT-1
        for(int i = 0; i < N; i++)
            betaMatrix[T-1][i] = c[T-1];

        // beta-pass
        for(int t = T-2; t>= 0; t--){
            for(int i = 0; i < N; i++){
                for(int j = 0; j < N; j++){
                    betaMatrix[t][i] += A[i][j]*B[j][O[t+1]]*betaMatrix[t+1][j];
                }
                // scale betat(i) with same factor as alphat(i)
                betaMatrix[t][i] *= c[t];
            }
        }


        // Step 4, Compute digammat(i,j) and gammat(i)

        vector<vector<double> > gammaMatrix(T, vector<double>(N, 0));
        vector <vector<vector<double> > > diGammaMatrix(T, vector < vector <double> >(N, vector<double>(N, 0)));

        for(int t = 0; t < T-1; t++){
            double denom = 0;
            for(int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                    denom += alphaMatrix[t][i] * A[i][j] * B[j][O[t+1]] * betaMatrix[t+1][j];
                }
            }
            for(int i = 0; i < N; i++){
                for(int j = 0; j < N; j++){
                    diGammaMatrix[t][i][j] = (alphaMatrix[t][i]*A[i][j]*B[j][O[t+1]]*betaMatrix[t+1][j])/denom;
                    gammaMatrix[t][i] += diGammaMatrix[t][i][j];
                }
            }
        }

        // special case for gammaT-1(i)
        double denom = 0;
        for (int i = 0; i < N; i++){
            denom += alphaMatrix[T-1][i];
        }
        for (int i = 0; i < N; i++){
            gammaMatrix[T-1][i] = alphaMatrix[T-1][i]/denom;
        }

        //Step 5, re-estimate A, B, pi

        // re-estimate pi
        for(int i = 0; i < N; i++)
            PI[i] = gammaMatrix[0][i];

        // re-estimate A
        for(int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                double numer = 0;
                double denom = 0;
                for(int t = 0; t < T-1; t++){
                    numer += diGammaMatrix[t][i][j];
                    denom += gammaMatrix[t][i];
                }
                A[i][j] = numer/denom;
            }
        }

        // re-estimate B
        for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++){
                double numer = 0;
                double denom = 0;
                for(int t = 0; t < T; t++){
                    if(O[t] == j)
                        numer += gammaMatrix[t][i];
                    denom += gammaMatrix[t][i];
                }
                B[i][j] = numer/denom;
            }
        }

        // Step 6 Compute log(P(O|lambda))
        double logProb = 0;
        for(int i = 0; i < T; i++){
            logProb += log(c[i]);
        }
        logProb = -logProb;

        //step 7
        iters++;
        if(iters < maxIters && logProb > oldLogProb){
            oldLogProb = logProb;

        }else{
            break;
        }

    }
    cout << N << " " << N << " ";
    printMatrixAsRows(A);
    cout << N << " " << M << " ";
    printMatrixAsRows(B);
}