#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Implementation of the Microarray Distortion Metric (MDM)
described in

Miguel Hernández-Cabronero, Victor Sanchez, Michael W. Marcellin, Joan Serra-Sagristà,
"A distortion metric for the lossy compression of DNA microarray images",
In proceedings of the IEEE International Data Compression Conference, DCC 2013.

This implementation requires that the image has been segmented into background
and spots, and assumes a segmentation template consisting of circles (x,y,radius).

--------------------------------------------------------------------------------
Usage
--------------------------------------------------------------------------------

This script can be invoked as a standalone program with the following syntax:

./distortion_metric.py <image1.raw> <image2.raw> <width> <height> <segmentation.csv>

It can also be included in other python scripts and used as follows:

distortion_metric = DistortionMetric(
    original_image_path,
    distorted_version_path,
    images_width,
    images_height,
    segmentation_template_path)
print distortion_metric.get_distortion()

--------------------------------------------------------------------------------
Copyright
--------------------------------------------------------------------------------

@author Miguel Hernández-Cabronero <mhernandez@deic.uab.cat> (http://deic.uab.cat/~mhernandez)
@date 15/05/2012

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import sys
import os
import time
import math

from binary_stream import *

############################ Begin configurable part

# Default alpha value to apply in the MDM definition
default_alpha = 3

# Be verbose?
be_verbose = True
be_superverbose = False

############################ Advanced configuration

# Given a circle center and radius produced by imfindcircles after rounding,
# these ofsets determine which pixels are considered foreground, edge and background.
# Example values of 0,1,2 mean that the circle with the original radius
# is considered foreground, a 1-pixel wide crown is considered edge and another
# 2-pixel wide crown is considered local background. To obtain the outermost
# crown, a circle the largest radius (radius+default_bg_offset) is created,
# and then a circle with the middle radius is subtracted. The same is done
# to produce the middle crown.
default_fg_offset = -1
default_edge_offset = 1
default_bg_offset = 5

# Default values for the foreground and background pixels, when handling
# artificial circles inside squares (used to build spot templates)
default_fg_value = 1
default_bg_value = 0

# Minimum and maximum radius values for the preloaded circle template offsets
min_template_radius = 2
max_template_radius = 25

# Check for overlapping spots in the segmentation templates?
# It can be very time consuming, so use only when the templates
# have not ben previously processed
remove_overlapping_spots = False
if remove_overlapping_spots == True:
    from remove_overlapping_circles import *

############################ End configurable part


class DistortionMetric:
    """Store the different distortion metrics measured between two DNA
    microarray images, and provide the get_distortion method to blend all
    data together into a single float value.
    """
    def __init__(self, image1_raw_path, image2_raw_path, width, height, segmentation_csv_path):
        """Obtain a DistortionMetric instance that represents the distortion present
        between two DNA microarray images. Images are processed inside the constructor,
        so that the get_distortion method can be called.
        """
        self.image1_raw_path = image1_raw_path
        self.image2_raw_path = image2_raw_path
        self.width = width
        self.height = height
        self.segmentation_csv_path = segmentation_csv_path
        # Build the circle template offsets and store them indexed by circle radius.
        # Each entry has format (spot_offsets, edge_offsets, background_offsets)
        self.template_offsets = build_template_offsets(min_template_radius, max_template_radius)
        # Read the segmentation from the CSV file created with the MATLAB script
        # One (x, y, r) tuple is generated per
        self.circle_positions = read_segmentation_csv(segmentation_csv_path, force_remove_overlapping=remove_overlapping_spots)
        # These are to be filled in the __calculate_distortion_metrics method.
        ## Indexed by the tuples in circle_positions, this dictionary contains
        ## tuples with the format
        ## (spot_mse, spot_ratio, edge_mse, edge_diff_var, local_bg_mse, local_bg_ratio)
        ## 0*spot_mse - MSE of the pixels tagged as spot
        ## 1*spot_ratio - ratio of the average intensity inside the spot for the original and distorted image
        ## 2*edge_mse - MSE of the pixels tagged as edge
        ## 3*edge_diff_var - variance of the differences in the pixels tagged as edge
        ## 4*local_bg_mse - MSE of the pixels tagged as local background
        ## 5*local_bg_ratio - ratio of the average intensity of the piuxels tagged as local background
        self.circle_distortions = dict()
        # MSE of all pixels not tagged as spot nor edge (but yes local background)
        self.background_mse = None
        # Ratio of the average intensities of all pixels not tagged as spot nor edge (but yes local background)
        self.background_ratio = None
        # Number of pixels not tagged as spot
        self.background_count = None
        # Global MSE, PSNR(16 bit) and SNR (20*log(Energy/MSE)) of the image
        self.global_mse = None
        self.global_psnr = None
        self.global_snr = None
        # Energy of each of the images
        self.img1_energy = None
        self.img2_energy = None
        # Intensity sum of the images
        self.img1_intensity_sum = None
        self.img2_intensity_sum = None
        # Maximum intensity in the two images
        self.max_value = None
        # Average spot and local background ratio
        self.avg_spot_ratio = None
        self.avg_local_bg_ratio = None
        self.min_spot_ratio = None
        self.max_spot_ratio = None
        self.min_local_bg_ratio = None
        self.max_local_bg_ratio = None
        self.max_spot_mse = None
        self.max_spot_mse_power_rspot_rlocalbg = None

        # Obtain the different distortion metrics
        self.__calculate_distortion_metrics(
            self.image1_raw_path,
            self.image2_raw_path,
            self.width,
            self.height,
            self.circle_positions)

    def __str__(self):
        return "[DistortionMetric(%s,%s):DM=%.3f,PSNR=%.3f]" % (
            os.path.basename(self.image1_raw_path.replace(".tif","").replace(".pgm","").replace(".raw","")),
            os.path.basename(self.image2_raw_path.replace(".tif","").replace(".pgm","").replace(".raw","")),
            self.get_distortion(),
            self.global_psnr)


    def get_distortion(self, alpha = default_alpha):
        """Blend the different distortion metrics into a single float value
        that represent the distortion among the two microarray images.

        The used formula is described in

        Miguel Hernández-Cabronero, Victor Sanchez, Michael W. Marcellin, Joan Serra-Sagristà,
        "A distortion metric for the lossy compression of DNA microarray images",
        In proceedings of the IEEE International Data Compression Conference, DCC 2013.

        MDM = 10 log10 (max_val)^2 / ME
        ME = (max_val)^p - max_val + min(max_val, MSE_global)
        p = 2 / (1 + exp(-alpha * (r_spot + r_localbg + r_global - 3)))
        """
        r_spot = max(self.max_spot_ratio, 1/float(self.min_spot_ratio))
        r_lbg = max(self.max_local_bg_ratio, 1/float(self.min_local_bg_ratio))
        r_global = max(self.img1_intensity_sum / float(self.img2_intensity_sum), self.img2_intensity_sum / float(self.img1_intensity_sum))

        t = r_spot + r_lbg + r_global
        p = 2 / (1 + math.exp(-alpha*(t-3)))
        ME = (self.max_value ** p) - self.max_value + min(self.max_value, self.global_mse)

        if ME == 0:
            MDM = float('inf')
        else:
            MDM = 10 * math.log10( ((2**16 - 1)**2) / float(ME))

        if be_verbose == True:
            print "@@@ MDM=%.2f for %s:\n\tspot(%.2f) localbg(%.2f) global(%.2f)" % (
                    MDM,
                    os.path.basename(self.image2_raw_path),
                    r_spot,
                    r_lbg,
                    r_global)

        return MDM

    def __calculate_distortion_metrics(self, image1_raw_path, image2_raw_path, width, height, circle_list):
        """Compare two dna microarray images and fill the variables of this instance
        """
        bin1 = BinaryStream.open(image1_raw_path, "r", precision_in_bytes=2, signed=False)
        bin2 = BinaryStream.open(image2_raw_path, "r", precision_in_bytes=2, signed=False)
        matrix1 = bin1.readMatrix(width, height)
        matrix2 = bin2.readMatrix(width, height)
        # 1 in pixels inside spots, 0 for pixels inside background
        spot_mask = [ [0 for y in range(height)] for x in range(width) ]

        mse = 0
        img1_intensity_sum = 0
        img2_intensity_sum = 0
        img1_energy = 0
        img2_energy = 0
        max_value = 0
        for x in range(self.width):
            for y in range(self.height):
                mse += (matrix1[x][y] - matrix2[x][y])*(matrix1[x][y] - matrix2[x][y])
                img1_energy += matrix1[x][y]*matrix1[x][y]
                img2_energy += matrix2[x][y]*matrix2[x][y]
                img1_intensity_sum += matrix1[x][y]
                img2_intensity_sum += matrix2[x][y]
                if matrix1[x][y] > max_value:
                    max_value = matrix1[x][y]
                if matrix2[x][y] > max_value:
                    max_value = matrix2[x][y]
        mse /= float(self.width*self.height)
        self.img1_energy = img1_energy
        self.img2_energy = img2_energy
        self.img1_intensity_sum = img1_intensity_sum
        self.img2_intensity_sum = img2_intensity_sum
        self.global_mse = mse
        self.max_value = max_value
        if self.global_mse != 0:
            self.global_psnr = 10*math.log10(((2**16-1)**2)/float(self.global_mse))
            self.global_snr = 20*math.log10(min(self.img1_energy, self.img2_energy)/float(self.global_mse))
        else:
            self.global_psnr = float("inf")
            self.global_snr = float("inf")

        # Process each of the spots to calculate the needed metrics
        # (spot_mse, spot_ratio, edge_mse, edge_diff_var, local_bg_mse, local_bg_ratio)
        for circle_tuple in circle_list:
            x,y,radius = int(round(circle_tuple[0])), int(round(circle_tuple[1])), int(round(circle_tuple[2]))
            (spot_offsets, edge_offsets, bg_offsets) = self.template_offsets[radius]

            # Calculate spot MSE and spot ratio
            spot_mse = 0
            avg_img1 = 0
            avg_img2 = 0
            skipped_offsets = 0
            for offset in spot_offsets:
                try:
                    value1 = matrix1[x+offset[0]][y+offset[1]]
                    value2 = matrix2[x+offset[0]][y+offset[1]]
                    spot_mse += (value1-value2)*(value1-value2)
                    avg_img1 += value1
                    avg_img2 += value2
                    spot_mask[x+offset[0]][y+offset[1]] = 1
                except IndexError:
                    skipped_offsets += 1
            spot_mse /= float(len(spot_offsets) - skipped_offsets)
            avg_img1 /= float(len(spot_offsets) - skipped_offsets)
            avg_img2 /= float(len(spot_offsets) - skipped_offsets)
            if avg_img1 == 0 or avg_img2 == 0:
                # Very rare case: divide by zero. The ratio is undefined so it is ignored
                spot_ratio = 1
            elif avg_img1 >= avg_img2:
                spot_ratio = avg_img1 / avg_img2
            else:
                spot_ratio = avg_img2 / avg_img1

            # Calculate edge MSE, edge ratio and edge variance
            edge_mse = 0
            avg_img1 = 0
            avg_img2 = 0
            skipped_offsets = 0
            for offset in edge_offsets:
                try:
                    value1 = matrix1[x+offset[0]][y+offset[1]]
                    value2 = matrix2[x+offset[0]][y+offset[1]]
                    edge_mse += (value1-value2)*(value1-value2)
                    avg_img1 += value1
                    avg_img2 += value2
                    if spot_mask[x+offset[0]][y+offset[1]] == 1:
                        if be_verbose == True and be_superverbose == True:
                            print ">>>", x+offset[0], y+offset[1], "is spot and edge !!!"
                except IndexError:
                    skipped_offsets += 1
            edge_mse /= float(len(edge_offsets) - skipped_offsets)
            avg_img1 /= float(len(edge_offsets) - skipped_offsets)
            avg_img2 /= float(len(edge_offsets) - skipped_offsets)
            if avg_img1 >= avg_img2:
                edge_ratio = avg_img1 / avg_img2
            else:
                edge_ratio = avg_img2 / avg_img1
            # X := differences ==> Var(X) = E[X^2] - E[X]^2 = MSE - E[X]^2
            edge_diff_var = edge_mse - (avg_img1-avg_img2)*(avg_img1-avg_img2)


            # Calculate local background MSE and local background ratio
            local_bg_mse = 0
            avg_img1 = 0
            avg_img2 = 0
            skipped_offsets = 0
            for offset in bg_offsets:
                try:
                    value1 = matrix1[x+offset[0]][y+offset[1]]
                    value2 = matrix2[x+offset[0]][y+offset[1]]
                    local_bg_mse += (value1-value2)*(value1-value2)
                    avg_img1 += value1
                    avg_img2 += value2
                except IndexError:
                    skipped_offsets += 1
            local_bg_mse /= float(len(bg_offsets) - skipped_offsets)
            avg_img1 /= float(len(bg_offsets) - skipped_offsets)
            avg_img2 /= float(len(bg_offsets) - skipped_offsets)
            if avg_img1 >= avg_img2:
                local_bg_ratio = avg_img1 / avg_img2
            else:
                local_bg_ratio = avg_img2 / avg_img1

            self.circle_distortions[circle_tuple] = \
                (spot_mse, spot_ratio, edge_mse, edge_diff_var, local_bg_mse, local_bg_ratio)


        # Calculate background MSE, background ratio and number of background pixels
        bg_mse = 0
        bg_count = 0
        avg_img1 = 0
        avg_img2 = 0
        for x in range(width):
            for y in range(height):
                if spot_mask[x][y] == 0:
                    bg_count += 1
                    value1 = matrix1[x][y]
                    value2 = matrix2[x][y]
                    bg_mse += (value1-value2)*(value1-value2)
                    avg_img1 += value1
                    avg_img2 += value2
        bg_mse /= float(bg_count)
        avg_img1 /= float(bg_count)
        avg_img2 /= float(bg_count)

        self.background_mse = bg_mse
        if avg_img1 >= avg_img2:
            self.background_ratio = avg_img1 / avg_img2
        else:
            self.background_ratio = avg_img2 / avg_img1
        self.background_count = bg_count

        # Calculate average values
        distortion = 0
        distortion_spots = 0
        distortion_edges = 0
        distortion_local_bg = 0
        ratio_spots = 0
        ratio_local_bg = 0
        self.min_spot_ratio = None
        self.max_spot_ratio = None
        self.min_local_bg_ratio = None
        self.max_local_bg_ratio = None
        self.max_spot_mse = None
        self.max_spot_mse_power_rspot_rlocalbg = None
        for circle_tuple in self.circle_positions:
            circle_distortions = self.circle_distortions[circle_tuple]
            # Avg values
            distortion_spots += circle_distortions[0]
            distortion_edges += circle_distortions[2]
            distortion_local_bg += circle_distortions[4]
            ratio_spots += circle_distortions[1]
            ratio_local_bg += circle_distortions[5]
            # Min/man ratios
            if self.min_spot_ratio == None:
                self.min_spot_ratio = circle_distortions[1]
                self.max_spot_ratio = circle_distortions[1]
                self.min_local_bg_ratio = circle_distortions[5]
                self.max_local_bg_ratio = circle_distortions[5]
                self.max_spot_mse = circle_distortions[0]
                r_spot = max(circle_distortions[1], 1/float(circle_distortions[1]))
                r_lbg = max(circle_distortions[5], 1/float(circle_distortions[5]))
                self.max_spot_mse_power_rspot_rlocalbg = circle_distortions[0]**((r_spot+r_lbg)/2.0)
            else:
                self.min_spot_ratio = min(self.min_spot_ratio, circle_distortions[1])
                self.max_spot_ratio = max(self.max_spot_ratio, circle_distortions[1])
                self.min_local_bg_ratio = min(self.min_local_bg_ratio, circle_distortions[5])
                self.max_local_bg_ratio = max(self.max_local_bg_ratio, circle_distortions[5])
                self.max_spot_mse = max(self.max_spot_mse, circle_distortions[0])
                self.max_spot_mse_power_rspot_rlocalbg = max(self.max_spot_mse_power_rspot_rlocalbg, circle_distortions[0]**((r_spot+r_lbg)/2.0))
        distortion = distortion_spots + distortion_edges + distortion_local_bg
        ratio_spots /= float(len(self.circle_positions))
        ratio_local_bg /= float(len(self.circle_positions))
        self.avg_spot_ratio = ratio_spots
        self.avg_local_bg_ratio = ratio_local_bg

processed_segmentations = dict()
def read_segmentation_csv(csv_path, force_remove_overlapping=False):
    """Read a CSV file containing lines with format

    x y radius

    (spaces are used as separators) and return a list of tuples (x,y,radius).
    
    If the remove_overlapping_spots configuration flag is True, overlapping
    circles will be removed from the returned list.
    """
    global processed_segmentations
    if csv_path in processed_segmentations:
        if be_verbose == True:
            print ">>> Skipping reading of", os.path.basename(csv_path), "[", len(processed_segmentations[csv_path]),"]", "(already processed)"
        return processed_segmentations[csv_path]

    f = open(csv_path, "r")
    segmentation = []

    # Read and process lines
    while True:
        line = f.readline().strip()
        if len(line) == 0:
            break
        if line.startswith("#"):
            # Ignore lines starting with '#'; they are comments
            continue

        elements = line.split(" ")
        tuple = (float(elements[0]), float(elements[1]), float(elements[2]))
        segmentation.append(tuple)

    if be_verbose == True:
        print "Obtained", len(segmentation), "raw spots"

    # Check that no spots are overlapping. Remove them if configured to do so
    remove_list = []
    if remove_overlapping_spots == True or force_remove_overlapping == True:
        segmentation = remove_overlapping_circles(segmentation)

    processed_segmentations[csv_path] = segmentation
    return segmentation

def build_template_offsets(min_radius, max_radius):
    """Build the template offsets for different radius values
    """
    if (min_radius > max_radius):
        raise Exception("Min radius > max_radius")
    offsets = [None]*(max_radius+1)

    for radius in range(min_radius, max_radius+1):
        # Get the circle template matrices
        spot_matrix = get_circle_in_square(radius + default_fg_offset, 2*(radius+default_bg_offset) + 1)
        edge_matrix = get_circle_in_square(radius + default_edge_offset, 2*(radius+default_bg_offset) + 1)
        bg_matrix = get_circle_in_square(radius + default_bg_offset, 2*(radius+default_bg_offset) + 1)
        # Get the offset lists
        spot_offsets = get_circle_offsets(spot_matrix)
        edge_offsets = get_circle_offsets(edge_matrix)
        bg_offsets = get_circle_offsets(bg_matrix)
        # Prune the offset lists (e.g., remove the spot_offsets from the edge_offsets)
        for offset in spot_offsets:
            edge_offsets.remove(offset)
            bg_offsets.remove(offset)
        for offset in edge_offsets:
            bg_offsets.remove(offset)

        offsets[radius] = (spot_offsets, edge_offsets, bg_offsets)

    return offsets

def get_circle_in_square(radius, square_side=None, fg_value=default_fg_value, bg_value=default_bg_value):
    """Return a square matrix of radius square_side with a circle centered in the middle.
    Pixels inside the circle will have value fg_value and pixels outside
    will have bg_value. If square_side is None, 2*radius will be chosen
    """
    if square_side == None:
        square_side = 2*radius
    matrix = [ [bg_value for y in range(square_side)] for y in range(square_side) ]

    for x in range(square_side):
        for y in range(square_side):
            corrected_x = x+1-(square_side/2)
            corrected_y = y+1-(square_side/2)
            if corrected_x*corrected_x + corrected_y*corrected_y < radius*radius:
                matrix[x][y] = fg_value

    return matrix

def get_circle_offsets(matrix, fg_value=default_fg_value, bg_value=default_bg_value):
    """Return a list of tuples (x_offset, y_offset) with the offsets calculated
    from the (square) matrix center."""
    offset_list = []
    square_side = len(matrix)

    for x in range(square_side):
        for y in range(square_side):
            corrected_x = x-(square_side/2)
            corrected_y = y-(square_side/2)
            if matrix[x][y] == fg_value:
                offset_list.append((corrected_x, corrected_y))

    return offset_list


############################ Begin main executable part

def show_help(message=""):
    message = message.strip()
    if message != "":
        print "-"*len(message)
        print message
        print "-"*len(message)
    print "Usage:", os.path.basename(sys.argv[0]), "<image1.raw> <image2.raw> <width> <height> <segmentation.csv>"
    print "\tCalculate distortion between image1.raw and image2.raw using the given segmentation file"


if __name__ == "__main__":
    if len(sys.argv) != 6:
        show_help("Error! Incorrect argument count")
        exit(1)
    image1 = sys.argv[1]
    image2 = sys.argv[2]
    try:
        width = int(sys.argv[3])
        height = int(sys.argv[4])
        if width <= 0 or height <= 0:
            raise ValueError("Non-positive width or height")
    except ValueError:
        show_help("Error! width and height values must be positive integers")
        exit(1)
    csv_path = sys.argv[5]
    if os.path.exists(image1) == False or \
            os.path.exists(image2) == False or \
            os.path.exists(csv_path) == False:
        show_help("Error! All input files must exist")
        exit(1)

    distortion = DistortionMetric(image1, image2, width, height, csv_path)
    print "Distortion:", distortion.get_distortion()
