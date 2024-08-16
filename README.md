## Tasks:

## Implement 1 filter analog design 
## Implement transformation to different filters
## Implement digital filter design using bilinear transform
## Implement actual topologies of filter 
## Implement fixed point implementation of that filter 

## Code of conduct.

### 1. **Use Descriptive Names:**
   - Choose names that clearly describe the purpose or role of the variable, function, or type.
   - Example: Use `buffer_size` instead of `bs`, `calculate_area()` instead of `ca()`.

### 2. **Variable Naming:**
   - Use `snake_case` for variable names.
   - Keep the names concise but descriptive.
   - Example: `int buffer_length;`, `char file_name[256];`

### 3. **Function Naming:**
   - Use `snake_case` for function names.
   - Start function names with a verb to indicate action.
   - Example: `read_file()`, `calculate_sum()`

### 4. **Constants and Macros:**
   - Use `UPPER_CASE` with underscores for constants and macros.
   - Example: `#define MAX_BUFFER_SIZE 1024`, `const int DEFAULT_TIMEOUT = 30;`

### 5. **Type Definitions:**
   - Use `PascalCase` (also called `CamelCase` with an initial uppercase letter) for `typedef` names, structs, and enums.
   - Example: `typedef struct FileHeader FileHeader;`, `enum Color { Red, Green, Blue };`

### 6. **Prefixing:**
   - Use consistent prefixes to group related functions and variables, especially in large codebases or libraries.
   - Example: `file_read()`, `file_write()`, `file_close()` (all related to file operations).

### 7. **Avoid Hungarian Notation:**
   - Avoid prefixing variable names with type information (`iCount`, `strName`). This is less common in modern C practice.
   - Instead, rely on meaningful names and comments to clarify usage.

### 8. **Abbreviations and Acronyms:**
   - Avoid unnecessary abbreviations. If you must use them, ensure they are widely understood.
   - Example: `calculate_average()` is better than `calc_avg()` unless the context is very clear.

### 9. **Global Variables:**
   - Use a prefix (often related to the module or context) for global variables to avoid name collisions.
   - Example: `g_config_file_path`, `app_error_code`

### 10. **Error Handling:**
   - Use consistent naming for error codes and handling functions.
   - Example: `int open_file();`, `int close_file();`, `return FILE_ERROR_OPEN_FAILED;`
