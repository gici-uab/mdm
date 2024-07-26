#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Tools to read file with a list of circles (x,y,radius), and remove
as many circles as necessary so that none is overlapping.

@author Miguel Hernández-Cabronero <mhernandez@deic.uab.cat> (http://deic.uab.cat/~mhernandez)
@date 15/05/2012
"""
#--------------------------------------------------------------------------------
#Copyright
#--------------------------------------------------------------------------------

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import math

############################ Begin configurable part

# Be verbose?
be_verbose = True
be_superverbose = False

############################ End configurable part

def remove_overlapping_circles_from_file(input_csv_path, output_csv_path):
    """Process the input csv file containing lines describing circles with format

    x y radius

    then remove the overlapping circles and then output the corresponding list
    to the output path with the same format.
    """
    fi = open(input_csv_path, "r")
    circle_list = []

    # Read and process lines
    while True:
        line = fi.readline().strip()
        if len(line) == 0:
            break
        if line.startswith("#"):
            # Ignore lines starting with '#'; they are comments
            continue

        elements = line.split(" ")
        tuple = (float(elements[0]), float(elements[1]), float(elements[2]))
        circle_list.append(tuple)

    if be_verbose == True:
        print "Obtained", len(circle_list), "raw spots"

    # Check that no spots are overlapping. Remove them if configured to do so
    circle_list = remove_overlapping_circles(circle_list)

    # Output the non-overlapping circles to the output file
    fo = open(output_csv_path, "w")
    fo.write("# Removed overlapping spots\n")
    for circle in circle_list:
        fo.write("%f %f %f\n" % (circle[0], circle[1], circle[2]))

def remove_overlapping_circles(spot_positions):
    """Read a list of (x,y,radius) circles and return a new list with those that are
    not overlapping. Two circles are considered overlapping when the distance
    between their centers (x1,y1) (x2,y2) is not larger than radius_1+radius2.
    The original list is not modified with this method.

    WARNING: this method can be very time-consuming!
    """
    remove_list = []

    segmentation = list(spot_positions)
    for position in segmentation:
        x1,y1,r1 = int(round(position[0])), int(round(position[1])), int(round(position[2]))
        for other_position in segmentation:
            if position is other_position:
                continue
            x2,y2,r2 = int(round(other_position[0])), int(round(other_position[1])), int(round(other_position[2]))
            distance = math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
            if r1 + r2 >= distance:
                remove_list.append(other_position)

                if be_verbose == True:
                    sys.stdout.write(".")
                    sys.stdout.flush()

    for removable in remove_list:
        try:
            segmentation.remove(removable)
            if be_superverbose == True:
                print ">> Ok removing", removable
        except ValueError as ex:
            if be_superverbose == True:
                print ">> Error removing", removable, ":", ex
            pass
    if be_verbose == True:
        print "\nErased", len(remove_list), "overlapping spots"

    return segmentation

############################ Begin main executable part

def show_help(message=""):
    message = message.strip()
    if message != "":
        print "-"*len(message)
        print message
        print "-"*len(message)
    print "Usage:", os.path.basename(sys.argv[0]), "<input_csv>", "<output_csv>"
    print
    for line in __doc__.split("\n"):
        line = line.strip()
        if len(line) == 0 or line[0] == "@":
            continue
        print line


if __name__ == "__main__":
    if len(sys.argv) != 3:
        show_help("Incorrect argument count")
	exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    if os.path.exists(input_path) == False:
        show_help("input_csv must exists")
        exit(1)

    remove_overlapping_circles_from_file(input_path, output_path)
