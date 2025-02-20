import os
import re
import shutil
import subprocess

def replace_in_file(file_path, old_str, new_str):
    """
    Function to replace a string in the specified file.
    """
    try:
        # Open the file in read mode and read its content
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Perform the replacement
        if old_str in content:
            content = content.replace(old_str, new_str)
            # Open the file in write mode and save the updated content
            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(content)
            print(f"Processed: {file_path}")
        else:
            print(f"No match found in: {file_path}")
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

def process_directory(directory, file_extension, old_str, new_str):
    """
    Function to recursively search for files in the specified directory and perform replacements.
    """
    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        return

    # Recursively search for all files in the directory
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_extension):
                file_path = os.path.join(root, file)
                replace_in_file(file_path, old_str, new_str)

def copy_directory(src_dir, dst_dir):
    """
    Recursively copy a directory and its contents to a new location.
    """
    # Check if the destination directory already exists
    if os.path.exists(dst_dir):
        print(f"Destination directory '{dst_dir}' already exists. Skipping copy.")
        return

    if not os.path.isdir(src_dir):
        print(f"Error: Source directory '{src_dir}' does not exist.")
        return
    try:
        # Copy the directory and its contents recursively
        shutil.copytree(src_dir, dst_dir)
        print(f"Successfully copied '{src_dir}' to '{dst_dir}'.")
    except FileExistsError:
        print(f"Error: Destination directory '{dst_dir}' already exists.")
    except Exception as e:
        print(f"An error occurred: {e}")

def replace_make(file_path):
    """
    'ma_XXXX.mod' -> 'mbar_XXXX.mod'
    """
    pattern = r"ma_([^\s]+)\.mod"
    replacement = r"mbar_\1.mod"

    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    new_content = re.sub(pattern, replacement, content)

    if content != new_content:
        with open(file_path, 'w', encoding='utf-8') as file:
            file.write(new_content)
        print(f"Updated content in file: {file_path}")
    else:
        print(f"No changes needed in file: {file_path}")

def main():
    src_directory = "../../free_energy/mbar_analysis/"
    target_directory = "../mbar_analysis"
    # Copy the mbar_analysis source code
    copy_directory(src_directory, target_directory)

    # Settings for the target directory and replacement strings
    file_extension = ".fpp"
    old_string = "ma_"
    new_string = "mbar_"

    # Execute the processing
    process_directory(target_directory, file_extension, old_string, new_string)

    # replace Makefile.depends
    replace_make(target_directory + "/Makefile.depends")

    print("All .fpp files have been processed.")

if __name__ == "__main__":
    main()
