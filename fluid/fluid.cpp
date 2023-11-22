//
#include "../sim/initSimulation.hpp"
#include "../sim/progargs.hpp"

#include <span>
#include <vector>

//#include <chrono>

using namespace std;

int main(int argc, char *argv[]){
    span const args_view{argv, static_cast<std::size_t>(argc)};
    vector<string> const arguments{args_view.begin() + 1, args_view.end()};
    validateParameters(arguments);
    initiateSimulation({arguments[0], arguments[1], arguments[2]});
}