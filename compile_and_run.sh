#!/bin/bash

# Define source files
SOURCE_FILES="main.c AnalogFilters.c DigitalFilters.c BillinearTransform.c Utils.c Filtering.c Tests.c"

# Define output folder and executable name
OUTPUT_DIR="output"
OUTPUT_EXECUTABLE="$OUTPUT_DIR/main"

# Create the output directory (if it doesn't exist already)
mkdir -p $OUTPUT_DIR

# Compile the source files and place them in the output directory
gcc -fdiagnostics-color=always -g $SOURCE_FILES -o $OUTPUT_EXECUTABLE

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Executable: $OUTPUT_EXECUTABLE"

    # Run the executable from the output directory
    ./$OUTPUT_EXECUTABLE
else
    echo "Compilation failed."
fi
