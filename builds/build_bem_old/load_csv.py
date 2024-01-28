import csv
import math

# Define the file path, adjust this as needed
file_path = 'SPHERIC_TestCase12/Flap.dat'


def parse_file(file_path):
    # Detect the file format
    if file_path.endswith('.csv'):
        delimiter = ','
    elif file_path.endswith('.dat'):
        delimiter = '\t'  # Changed to tab delimiter

    # Open the CSV file
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimiter)

        # Define the output string
        output = "const std::vector<Tdd> sample = {\n"

        i = 0
        # Loop through each row in the CSV
        for row in csv_reader:
            i += 1
            if i % 20 is not 0:
                continue

            # If this row is a comment or empty, skip it
            if not row or row[0].startswith('#'):
                continue

            try:
                # Append this row to the output string
                output += '          {{{}, {}}},\n'.format(
                    float(row[0]), float(row[1]))
            except (ValueError, IndexError) as e:
                print(f"Error processing line: {row}. Error: {e}")
                continue

        # Close the brackets for the vector declaration and add semicolon
        output += '};'

        # Print the output string
        print(output)


parse_file(file_path)
