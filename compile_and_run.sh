
#!/bin/bash

# Define source files
SOURCE_FILES="main.c AnalogFilters.c DigitalFilters.c BillinearTransform.c"

# Define output executable name
OUTPUT_EXECUTABLE="main"

# Compile the source files
gcc -fdiagnostics-color=always -g $SOURCE_FILES -o $OUTPUT_EXECUTABLE

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Executable: $OUTPUT_EXECUTABLE"
    
    # Run the executable
    ./$OUTPUT_EXECUTABLE

else
    echo "Compilation failed."
fi
