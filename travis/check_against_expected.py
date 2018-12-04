import csv
import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("--test-errors", type=str, required=True)
parser.add_argument("--test-cross-section", type=str, required=True)

args = parser.parse_args()

# Set test tolerance.
tolerance = 1e-10

# Read in expected values
expected_values = []
with open("travis/smoke_tests_expected.out", "r") as infile:
    csvreader = csv.reader(infile, delimiter=",", quotechar='"')
    for row in csvreader:
        if len(row) > 0:
            expected_values.append(float(row[-1]))

# Read in computed values
computed_values = []
with open("smoke_tests.out", "r") as infile:
    csvreader = csv.reader(infile, delimiter=",", quotechar='"')
    for row in csvreader:
        if len(row) > 0:
            computed_values.append(float(row[-1]))

# Check that the right number of tests are there.
if len(expected_values) != len(computed_values):
    print("Wrong number of experiments were calculated!")
    sys.exit(1)

# Compare expected to computed values
failed = False
for i in range(len(expected_values)):
    if (abs(computed_values[i]-expected_values[i])/expected_values[i] > tolerance):
        print("Test {} Failed! Expected {} got {}!".format(i, expected_values[i], computed_values[i]))
        failed = True

if failed:
    print("Test base values fail!")
    sys.exit(1)

print("Test base values pass!")

# Read in expected values
expected_values = []
with open("travis/smoke_tests_expected_1.out", "r") as infile:
    csvreader = csv.reader(infile, delimiter=" ", quotechar='"')
    for row in csvreader:
        if len(row) > 0:
            expected_values.append(float(row[-1]))

# Read in computed values
computed_values = []
with open("smoke_tests_1.out", "r") as infile:
    csvreader = csv.reader(infile, delimiter=" ", quotechar='"')
    for row in csvreader:
        if len(row) > 0:
            computed_values.append(float(row[-1]))

# Check that the right number of tests are there.
if len(expected_values) != len(computed_values):
    print("Wrong number of cross section values were calculated!")
    sys.exit(1)

# Compare expected to computed values
failed = False
for i in range(len(expected_values)):
    if (abs(computed_values[i]-expected_values[i])/expected_values[i] > tolerance):
        print("Value {} Failed! Expected {} got {}!".format(i, expected_values[i], computed_values[i]))
        failed = True

if failed:
    print("Test cross section values fail!")
    sys.exit(1)

print("Test cross section values pass!")
