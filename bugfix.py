import json

# Read data from the first file
with open('list1.txt', 'r') as file1:
    list1 = json.loads(file1.read())

# Read data from the second file
with open('list2.txt', 'r') as file2:
    list2 = json.loads(file2.read())

# Extract the second elements from the first list
list1_second_elements = [item[1] for item in list1]

# Extract the second elements from the 'coordinates' in the second list
list2_second_elements = [item['coordinates'][1] for item in list2]

# Compare the elements
comparison_results = [a == b for a, b in zip(list1_second_elements, list2_second_elements)]

# Print the comparison results
for i, (a, b, result) in enumerate(zip(list1_second_elements, list2_second_elements, comparison_results)):
    print(f"Index {i}: List1 second element = {a}, List2 'coordinates' second element = {b}, Match = {result}")
