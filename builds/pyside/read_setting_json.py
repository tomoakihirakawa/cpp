import json

types = {}

# import json file

dir = './input_files/Tanizwa1996_H0d05_L1d8_piston/'

with open(dir+'setting.json') as f:
    data = json.load(f)
    # print(data)
    for key, value in data.items():
        # print(key, value)
        if key == "input_files":
            for file in value:
                with open(dir+file) as ff:
                    data = json.load(ff)
                    # print(data)
                    for key, value in data.items():
                        # print(key, value)
                        if key == "type":
                            # add if value is not in types
                            if value not in types:
                                types[value] = []
                            types[value].append(data)
                            # print(types)
                
# display types
for key, value in types.items():
    print(key, len(value))
    for v in value:
        print(v)
    print()