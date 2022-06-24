import sys, json;
import time;

# Reading input.json
inputFile = open(sys.argv[1] + '/input.json')
data = json.load(inputFile)

intIn = data['num']
print(intIn)
faces = [':)', ':(', ':P', ':/', 'xD', ':O']

idx = int(intIn % len(faces))
print(idx)
result = {
    "your_face": faces[len(faces) - 1 if idx == 0 else idx - 1]
}
json_object = json.dumps(result, indent = 2)

# Writing output.json
with open(sys.argv[1] + '/output.json', "w") as outfile:
    outfile.write(json_object)
