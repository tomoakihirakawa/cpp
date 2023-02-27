import json
import sys
red = '\033[91m'
green = '\033[92m'
magenta = '\033[95m'
coloroff = '\033[0m'


input_directory = str(sys.argv[1])
print('input_directory:', input_directory)

input_files = []


def readJson(setting, color=red):
    global input_files
    f = open(setting)
    data = json.load(f)
    print('------------------------------------')
    for key, value in data.items():
        print(f'{key: <{20}}', '\t', color, value, coloroff)
        if "input_files" == key:
            input_files = value
    print('------------------------------------')
    f.close()


readJson(input_directory + '/setting.json')

# #@ -------------------------------------------------------- #
# #@           その他，water.json,tank.json などを出力           #
# #@ -------------------------------------------------------- #
for INPUTS in input_files:
    readJson(input_directory + "/" + INPUTS, green)

# print("The directory for input files :", magenta, input_directory, coloroff)
