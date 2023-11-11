//
#include <iostream>
#include <string>

void throwException(const std::string& errorMessage){
  throw std::runtime_error(errorMessage);
}
int exceptionHandler(bool condition, const std::string& errorMessage) {
  if (condition) {throwException(errorMessage);}
  return 0;
}

