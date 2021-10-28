#!/usr/bin/env python3
"""Split mummer delta file by chromosome

Usage:
    <program> input_delta output_name_stub
"""

# Modules
import gzip
import sys

# Classes
class Delta(object):
    """Delta object with chromosome name and associated lines
    """

    def __init__(self, name, lines):
        self.name = name
        self.lines = lines

    def write_to_file(self, handle):
        handle.write("\n".join(self.lines) + "\n")

    def __repr__(self):
        return self.name + "\n" + "\n".join(self.lines)

# Defining functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def delta_iterator(input_file):
    """Takes a delta file input_file and returns a delta iterator
    """
    with myopen(input_file) as f:
        lines = []
        name = ""
        begun = False

        for line in f:
            line = line.strip()

            if line.startswith("/"):
                continue

            elif line.startswith(">"):
                if begun:
                    yield Delta(name, lines)

                name = line[1:].split(" ")[0]
                lines = [line]
                begun = True

            else:
                lines.append(line)

        if name != "":
            yield Delta(name, lines)

# Parsing user input
try:
    input_delta = sys.argv[1]
    output_name_stub = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Iterate over input_delta
deltas = delta_iterator(input_delta)
output_files = dict()

for d in deltas:
    if d.name not in output_files:
        output_files[d.name] = myopen(output_name_stub + "_" + d.name + ".deltas.gz", "wt")

    d.write_to_file(output_files[d.name])
