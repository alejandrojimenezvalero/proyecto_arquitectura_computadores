#include "progargs.hpp"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int ProgArgs (const vector<string> &arguments){
    size_t n_args = arguments.size();
    ifstream inputFile(arguments[1]);
    ofstream outputFile(arguments[2]);
    if (n_args != 3){
        cerr << "Error: Invalid number of arguments: " << n_args << "." << endl;
        return -1;
    }try {
        int numSteps = stoi(arguments[0]);
        if (numSteps < 0) {
            cerr << "Error: Invalid number of time steps." << endl;
            return -2;
        }
    }catch (const invalid_argument &) {
        cerr << "Error: time steps must be numeric." << endl;
        return -1;
    }if (!inputFile) {
        cerr << "Cannot open "<< arguments[1] <<" for reading" << endl;
        return -3;
    }else if (!outputFile){
        cerr << "Cannot open "<< arguments[2] <<" for writing" << endl;
        return -4;
    }return 0;
}
//si el archivo de salida no tiene permisos de lectura y si de escritura, se mete en el error -4