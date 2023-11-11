//

#include "gtest/gtest.h"
#include "./sim/progargs.hpp"  // Reemplaza con el nombre de tu archivo que contiene validateParameters

#include <stdexcept>
#include <filesystem>

//Valid Test
TEST(ValidateParametersTest, ValidArguments) {
  std::vector<std::string> args = {"2000", "../../small.fld", "../../final.fld"};
  EXPECT_NO_THROW(validateParameters(args));
}

//Invalid Test: number of arguments must be 3
TEST(ValidateParametersTest, InvalidNumberOfArgumentsLessThanThree) {
  std::vector<std::string> args = {};
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::string expected_error_message = "Error: Invalid number of arguments: " + std::to_string(args.size()) + ".";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());
  }
}

//Invalid Test: number of arguments must be 4
TEST(ValidateParametersTest, InvalidNumberOfArgumentsMoreThanThree) {
  std::vector<std::string> args = {"2000", "../../small.fld", "../../final.fld", "45"};
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::string expected_error_message = "Error: Invalid number of arguments: " + std::to_string(args.size()) + ".";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());
  }
}

//Invalid Test: invalid number of time steps
TEST(ValidateParametersTest, NegativeTimeSteps) {
  std::vector<std::string> args = {"-5",  "../../small.fld", "../../final.fld"};
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::string expected_error_message = "Error: Invalid number of time steps.";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());
  }
}

//Invalid Test: Time steps must be numeric
TEST(ValidateParametersTest, InvalidTimeSteps) {
  std::vector<std::string> args = {"timeSteps",  "../../small.fld", "../../final.fld"};
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::string expected_error_message = "Error: time steps must be numeric.";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());
  }
}

//Invalid Test: invalid input file.
TEST(ValidateParametersTest, InvalidInputFile) {
  std::vector<std::string> args = {"2000",  "input.txt", "../../final.fld"};
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::string expected_error_message = "Cannot open " + args[1]+  " for reading";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());
  }
}

//Invalid Test: invalid input file.
TEST(ValidateParametersTest, InvalidOutputFile) {
  std::vector<std::string> args = {"2000",  "../../small.fld", "output.txt"};
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::string expected_error_message = "Cannot open " + args[1]+  " for reading";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());
  }
}

//Invalid Test: input file valid but without reading permissions
TEST(ValidateParametersTest, InputFileWithoutPermission) {
  std::vector<std::string> args = {"2000", "../../small.fld", "../../final.fld"};
  std::filesystem::permissions(args[1], std::filesystem::perms::owner_read, std::filesystem::perm_options::remove);
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::filesystem::permissions(args[1], std::filesystem::perms::owner_read, std::filesystem::perm_options::add);
    std::string expected_error_message = "Cannot open " + args[1]+  " for reading";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());

  }

}

//Invalid Test: output file valid but without writing permissions
TEST(ValidateParametersTest, OutputFileWithoutPermission) {
  std::vector<std::string> args = {"2000", "../../small.fld", "../../final.fld"};
  std::filesystem::permissions(args[2], std::filesystem::perms::owner_write, std::filesystem::perm_options::remove);
  try {
    validateParameters(args);
  }
  catch (const std::runtime_error& e){
    std::filesystem::permissions(args[2], std::filesystem::perms::owner_write, std::filesystem::perm_options::add);
    std::string expected_error_message = "Cannot open " + args[1]+  " for writing";
    EXPECT_STREQ(expected_error_message.c_str(), e.what());
  }

}