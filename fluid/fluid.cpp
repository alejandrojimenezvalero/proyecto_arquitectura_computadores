//
#include <cstdlib>
#include "../sim/progargs.cpp"
#include <span>
#include <vector>

using namespace std;

int main(int argc, char *argv[]){
    span const args_view{argv, static_cast<std::size_t>(argc)};
    vector<string> const arguments{args_view.begin() + 1, args_view.end()};
    validateParameters(arguments);
    initiateSimulation(arguments[1]);
    return 0;
}