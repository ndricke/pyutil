#!/usr/bin/env python

import sys

f = open(sys.argv[1])

while True:
    line = f.readline()
    if not line:
        break
    if "Standard Nuclear Orientation" in line:
        f.readline()
        f.readline()
        lines = []
        while True:
            line = f.readline()
            if not line or "---------------------------------------" in line:
                break
            lines.append(" ".join(line.split()[1:]))
        print len(lines)
        print
        for l in lines:
            print l

