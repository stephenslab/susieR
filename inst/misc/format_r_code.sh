#!/bin/bash

# Function to display error messages and exit the script
display_error() {
    echo "Error: $1"
    exit 1
}

# Check if the R file is provided as an argument
if [ $# -eq 0 ]; then
    display_error "Please provide the path to the R file as an argument."
fi

# Get the R file path from the argument
r_file="$1"

# Check if the R file exists
if [ ! -f "$r_file" ]; then
    display_error "The specified R file '$r_file' does not exist."
fi

echo "Formatting R code in file: $r_file"

# Format the R code using formatR and capture the output
# output=$(echo "library(formatR); tryCatch(tidy_source(\"$r_file\", file = \"$r_file\", indent = 2, args.newline = TRUE, arrow = TRUE), error = function(e) {message(\"Error formatting R code:\"); print(e)})" | R --slave --no-save 2>&1)
# Format the R code using styler and capture the output
output=$(echo "tryCatch(styler::style_file(\"$r_file\"), error = function(e) {message(\"Error formatting R code:\"); print(e)})" | R --slave --no-save 2>&1)

# Check if the formatting was successful
if echo "$output" | grep -q "Error formatting R code:"; then
    echo "Formatting failed. Please check the R code for syntax errors."
    echo -e "\033[0;31m$(echo "$output" | head -n -3)\033[0m"
    echo -e "\033[1;31m$(echo "$output" | tail -n 3)\033[0m"
    exit 1
else
    echo "R code formatting completed successfully."
fi