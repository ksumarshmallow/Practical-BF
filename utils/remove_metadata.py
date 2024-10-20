import os
import sys

def remove_schedule_string(file_path):
    # Read file content
    with open(file_path, 'r') as file:
        content = file.read()

    # Define the start and end strings
    start_string = "');}/* [&lt;] Schedule"
    end_string = "minmax(0, 1fr));  }}"

    # Find the start and end indexes
    start_index = content.find(start_string)
    end_index = content.find(end_string) + len(end_string)
        
    # If both start and end strings are found, remove the substring
    if start_index != -1 and end_index != -1:
        new_content = content[:start_index] + content[end_index:]

        # Write the cleaned content back to the file
        with open(file_path, 'w') as file:
            file.write(new_content)

# Check if the directory is passed as a command line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <directory>")
    sys.exit(1)

# Get the directory path from command line arguments
filename = sys.argv[1]

# Loop through all markdown files in the directory and apply the changes
if filename.endswith(".html"):
    remove_schedule_string(filename)

print("Processing complete.")