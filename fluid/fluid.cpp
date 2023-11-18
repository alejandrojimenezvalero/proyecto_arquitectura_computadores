//
#include "../sim/initSimulation.hpp"
#include "../sim/progargs.hpp"

#include <span>
#include <vector>

#include <chrono>

using namespace std;

int main(int argc, char *argv[]){
    auto start = std::chrono::high_resolution_clock::now();

    span const args_view{argv, static_cast<std::size_t>(argc)};
    vector<string> const arguments{args_view.begin() + 1, args_view.end()};
    validateParameters(arguments);
    initiateSimulation(arguments[0], arguments[1]);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << duration.count()/1000000 << '\n';
}