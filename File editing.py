import os


def truncate_files(input_folder, output_folder, line_limit=1600):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Process each file in the input folder
    for filename in os.listdir(input_folder):
        input_file_path = os.path.join(input_folder, filename)

        # Check if it's a file
        if os.path.isfile(input_file_path):
            # Read the file contents
            with open(input_file_path, 'r') as file:
                lines = file.readlines()

            # Truncate lines after the line limit
            truncated_lines = lines[:line_limit]

            # Save the truncated content to the output folder
            output_file_path = os.path.join(output_folder, filename)
            with open(output_file_path, 'w') as file:
                file.writelines(truncated_lines)

            print(f"Processed {filename} and saved to {output_folder}")


# Folder paths
input_folder = r"C:\Users\ASUS\Desktop\Research\Results\Raman\13.12.2024 - GSH Series pH control\trial 3\raw"
output_folder = r"C:\Users\ASUS\Desktop\Research\Results\Raman\13.12.2024 - GSH Series pH control\trial 3"
# Run the function
truncate_files(input_folder, output_folder)
