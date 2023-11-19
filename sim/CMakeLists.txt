# Add to this list all files related to sim library
add_library(sim
        progargs.hpp
        progargs.cpp
        grid.cpp
        grid.hpp
        block.cpp
        block.hpp
        exceptionHandler.cpp
        exceptionHandler.hpp
        constants.hpp
        fileManager.cpp
        fileManager.hpp
        initSimulation.cpp
        initSimulation.hpp
        particle.cpp
        particle.hpp
        processSimulation.cpp
        processSimulation.hpp
        simulationResults.cpp
        simulationResults.hpp
)
# Use this line only if you have dependencies from stim to GSL
target_link_libraries(sim PRIVATE Microsoft.GSL::GSL)