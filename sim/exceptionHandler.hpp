//

#ifndef FLUID_EXCEPTIONHANDLER_HPP
#define FLUID_EXCEPTIONHANDLER_HPP
#include <iostream>
#include <string>

int exceptionHandler(bool condition, const std::string& errorMessage);
void throwException(const std::string& errorMessage);
#endif
