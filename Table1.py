import csv
import os

gamma_list = ['0', '0.25', '0.5', '0.75', '1.0']
type_list = ['nonconvex', 'convex']

print("Table 1")
for the_type in type_list:
    print(the_type)
    message = ""
    for gamma in gamma_list:
        filename = "Table1_{}_{}.out".format(the_type, gamma)
        if os.path.exists(filename):
            data = []
            with open(filename, 'r') as the_file:
                csvreader = csv.reader(the_file, delimiter=',', quotechar='"')
                for row in csvreader:
                    if row[4].strip() == "True":
                        data.append(float(row[3]))
            data = sorted(data)
            # Break into intervals
            ranges = []
            state = 0 # Looking to start a new range
            low = 0
            high = 0
            last = 0
            for item in data:
                if state == 0:
                    low = item
                    last = item
                    high = item
                    state = 1
                elif state == 1:
                    if (item-last) > 0.0015:
                        # Ending an interval
                        ranges.append([low, item])
                        low = item
                        last = item
                        high = item
                    else:
                        last = item
                        high = item
            if high-low > 0.0015:
                ranges.append([low,high])

            message += "{}: ".format(gamma)
            for sub_range in ranges:
                message += "{} ".format(sub_range)
    print(message)
