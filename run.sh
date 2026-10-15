#!/bin/bash

EXECUTABLE="./Lab1"
CHECKER="./SolutionChecker"

echo "=================================================="
echo "Test start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================="

echo "Cleaning old executables..."
make clean

echo "Compiling the latest source code..."
make

if [ $? -ne 0 ]; then
    echo "Compilation failed. Please check your source code."
    exit 1
fi

echo "Compilation successful."
echo "--------------------------------------------------"

echo "Running test cases..."

TEST_CASES=("case1.in" "case2.in" "case3.in" "case4.in")

for INPUT_FILE in "${TEST_CASES[@]}"; do
    echo "=================================================="
    
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Input file $INPUT_FILE not found. Skipping this test case."
        continue
    fi

    OUTPUT_FILE="${INPUT_FILE%.in}.out"
    
    echo "Running test case: $INPUT_FILE"
    echo "Output file: $OUTPUT_FILE"
    echo "Current time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "--------------------------------------------------"
    echo "Execution time and results:"

    time $EXECUTABLE $INPUT_FILE $OUTPUT_FILE

    if [ $? -eq 0 ]; then
        echo "Execution completed. Validating results with SolutionChecker:"
        if [ -f "$CHECKER" ]; then
            $CHECKER $INPUT_FILE $OUTPUT_FILE
        else
            echo "Warning: SolutionChecker executable not found in the current directory. Please ensure it exists and has execute permission."
        fi
    else
        echo "An error occurred during execution."
    fi
done

echo "=================================================="
echo "All test cases have been processed."
echo "Test end time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================="