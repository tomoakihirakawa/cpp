import csv

# Define the file path, adjust this as needed
# file_path = 'bspline_sample_data.dat'
file_path = 'lag_data.dat'

# Open the CSV file
with open(file_path, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')

    # Define the output string
    output = "const std::vector<Tdd> sample = {\n"

    # Loop through each row in the CSV
    for row in csv_reader:
        # If this row is a comment, skip it
        if row[0].startswith('#'):
            continue

        # Append this row to the output string
        output += '          {{{}, {}}},\n'.format(row[0], row[1])

    # Close the brackets for the vector declaration and add semicolon
    output += '};'

    # Print the output string
    print(output)
