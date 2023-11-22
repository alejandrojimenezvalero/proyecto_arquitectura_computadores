#include "gtest/gtest.h"
#include "./sim/progargs.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <sstream>

using namespace testing;

std::string getFileContent(const std::string &filename) {
    const std::ifstream file(filename);
    if (file) {
        std::ostringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }
    return "";
}
//Valid Test
TEST(ValidateParametersTest, ValidArguments) {
    const std::vector<std::string> args = {"2000", "../../small.fld", "../../output.fld"};
    const std::stringstream captured_stderr;
    std::streambuf *old_stderr = std::cerr.rdbuf(captured_stderr.rdbuf());

    const int exit_code = validateParameters(args);
    EXPECT_EQ(exit_code, 0);

    std::cerr.rdbuf(old_stderr);
    const std::string captured_message = captured_stderr.str();
    EXPECT_TRUE(captured_message.empty());
}

// Invalid Test: number of arguments must be 3
TEST(ValidateParametersTest, InvalidNumberOfArgumentsLessThanThree) {
    const std::vector<std::string> args = {};

    // Create a temporary file to capture stderr
    const std::string temp_filename = "temp_stderr.txt";
    std::ofstream temp_file(temp_filename);

    // Redirect stderr to the temporary file
    std::streambuf *original_stderr = std::cerr.rdbuf(temp_file.rdbuf());

    // Check exit code
    EXPECT_EXIT(validateParameters(args), ::testing::ExitedWithCode(255), ".*");

    // Restore stderr
    std::cerr.rdbuf(original_stderr);

    temp_file.close();

    // Read the captured stderr from the temporary file
    const std::string expected_error_message =
            "Error: Invalid number of arguments: " + std::to_string(args.size()) + ".";
    std::string captured_message = getFileContent(temp_filename);
    captured_message.erase(std::remove(captured_message.begin(), captured_message.end(), '\n'), captured_message.end());

    // Check if the expected error message is equal to the captured message
    EXPECT_EQ(expected_error_message, captured_message);

    // Remove the temporary file
    std::filesystem::remove(temp_filename);
}

//Invalid Test: number of arguments must be 3
TEST(ValidateParametersTest, InvalidNumberOfArgumentsMoreThanThree) {
    const std::vector<std::string> args = {"2000", "../../small.fld", "../../output.fld", "45"};
    // Create a temporary file to capture stderr
    const std::string temp_filename = "temp_stderr.txt";
    std::ofstream temp_file(temp_filename);

    // Redirect stderr to the temporary file
    std::streambuf *original_stderr = std::cerr.rdbuf(temp_file.rdbuf());

    // Check exit code
    EXPECT_EXIT(validateParameters(args), ::testing::ExitedWithCode(255), ".*");

    // Restore stderr
    std::cerr.rdbuf(original_stderr);

    temp_file.close();

    // Read the captured stderr from the temporary file
    const std::string expected_error_message =
            "Error: Invalid number of arguments: " + std::to_string(args.size()) + ".";
    std::string captured_message = getFileContent(temp_filename);
    captured_message.erase(std::remove(captured_message.begin(), captured_message.end(), '\n'),
                           captured_message.end());

    // Check if the expected error message is equal to the captured message
    EXPECT_EQ(expected_error_message, captured_message);

    // Remove the temporary file
    std::filesystem::remove(temp_filename);
}

//Invalid Test: invalid number of time steps
TEST(ValidateParametersTest, NegativeTimeSteps) {
    const std::vector<std::string> args = {"-5", "../../small.fld", "../../final.fld"};
    // Create a temporary file to capture stderr
    const std::string temp_filename = "temp_stderr.txt";
    std::ofstream temp_file(temp_filename);

    // Redirect stderr to the temporary file
    std::streambuf *original_stderr = std::cerr.rdbuf(temp_file.rdbuf());

    // Check exit code
    EXPECT_EXIT(validateParameters(args), ::testing::ExitedWithCode(254), ".*");

    // Restore stderr
    std::cerr.rdbuf(original_stderr);

    temp_file.close();
    const std::string expected_error_message = "Error: Invalid number of time steps.";
    std::string captured_message = getFileContent(temp_filename);
    captured_message.erase(std::remove(captured_message.begin(), captured_message.end(), '\n'),
                           captured_message.end());

    // Check if the expected error message is equal to the captured message
    EXPECT_EQ(expected_error_message, captured_message);

    // Remove the temporary file
    std::filesystem::remove(temp_filename);

}
//Invalid Test: invalid input file because it not exists
TEST(ValidateParametersTest, CannotOpenInputFile) {
    const std::vector<std::string> args = {"2000", "../../input.fld", "../../output.fld"};
    try {
        validateParameters(args);
    }
    catch (const std::ifstream::failure &e) {
        const std::string expected_error_message = "Cannot open " + args[1] + " for reading";
        EXPECT_STREQ(expected_error_message.c_str(), e.what());

    }
}

//Invalid Test: input file valid but without reading permissions
TEST(ValidateParametersTest, InputFileWithoutPermission) {
    const std::vector<std::string> args = {"2000", "../../small.fld", "../../output.fld"};
    std::filesystem::permissions(args[1], std::filesystem::perms::owner_read, std::filesystem::perm_options::remove);
    try {
        validateParameters(args);
    }
    catch (const std::ifstream::failure &e) {
        std::filesystem::permissions(args[1], std::filesystem::perms::owner_read, std::filesystem::perm_options::add);
        const std::string expected_error_message = "Cannot open " + args[1] + " for reading";
        EXPECT_STREQ(expected_error_message.c_str(), e.what());

    }

}
//Invalid Test: output file valid but without writing permissions
TEST(ValidateParametersTest, OutputFileWithoutPermission) {
    const std::vector<std::string> args = {"2000", "../../small.fld", "../../output.fld"};
    std::filesystem::permissions(args[2], std::filesystem::perms::owner_write, std::filesystem::perm_options::remove);
    try {
        validateParameters(args);
    }
    catch (const std::ifstream::failure &e) {
        std::filesystem::permissions(args[2], std::filesystem::perms::owner_write, std::filesystem::perm_options::add);
        const std::string expected_error_message = "Cannot open " + args[1] + " for writing";
        EXPECT_STREQ(expected_error_message.c_str(), e.what());
    }

}