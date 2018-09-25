#include <cstdio>
#include <vector>
#include <iostream>

using namespace std;
int main()
{
    //TRANSITION MATRIX
    int transitionRows;    //The states
    int transitionColumns; //The transition probabilities for each state
    cin >> transitionRows >> transitionColumns;
    vector<vector<double>> transitionMatrix(transitionRows);
    for (int i = 0; i < transitionRows; i++)
    {
        for (int j = 0; j < transitionColumns; j++)
        {
            double element;
            cin >> element;
            transitionMatrix[i].push_back(element);
        }
    }

    //EMISSION MATRIX
    int emissionRows;
    int emissionColumns;
    cin >> emissionRows >> emissionColumns;
    vector<vector<double>> emissionMatrix(emissionRows);
    for (int i = 0; i < emissionRows; i++)
    {
        for (int j = 0; j < emissionColumns; j++)
        {
            double element;
            cin >> element;
            emissionMatrix[i].push_back(element);
        }
    }

    //INITIAL STATE MATRIX
    int initialStateRows;
    int initialStateColumns;
    cin >> initialStateRows >> initialStateColumns;
    vector<double> initialStateMatrix;
    for (int i = 0; i < initialStateColumns; i++)
    {
        double element;
        cin >> element;
        initialStateMatrix.push_back(element);
    }

    vector<double> currentStateMatrix;

    if (initialStateColumns != transitionRows || transitionColumns != emissionRows)
    {
        cout << "Incorrect sizes \n";
    }
    else
    {
        vector<double> currentStateMatrix;
        vector<double> currentEmissionMatrix;
        //cout << "Correct sizes \n";
        //Calculate pi1
        for (int i = 0; i < transitionColumns; i++)
        {
            double newValue = 0.0;
            for (int j = 0; j < transitionRows; j++)
                newValue += initialStateMatrix[j] * transitionMatrix[j][i];
            currentStateMatrix.push_back(newValue);
        }
        //for(int i = 0; i < transitionColumns; i++)
        //cout << currentStateMatrix[i] << '\n';
        //Calculate B1
        for (int i = 0; i < emissionColumns; i++)
        {
            double newValue = 0.0;
            for (int j = 0; j < emissionRows; j++)
                newValue += currentStateMatrix[j] * emissionMatrix[j][i];
            currentEmissionMatrix.push_back(newValue);
        }
        cout << 1 << ' ' << emissionColumns << ' ';
        for (int i = 0; i < emissionColumns; i++)
            cout << currentEmissionMatrix[i] << ' ';
    }
}