//
#include <iostream>
#include <string>

void throwException(const std::string& errorMessage, int errorCode){
  std::cerr <<  errorMessage << std::endl;
  exit(errorCode);

};
int exceptionHandler(bool condition, const std::string& errorMessage, int errorCode) {
  if (condition) {throwException(errorMessage, errorCode);}
  return 0;
}

