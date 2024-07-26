#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Binary stream implementation in python: read chars/bytes.

Example usage:

# Open unsigned, 16bit raw image and get a matrix (list of lists) with the pixel values
bin = BinaryStream.open(filePath, "r", precision_in_bytes=2, signed=False)
matrix = bin.readMatrix(width, height)

# Read a value the fastest way (useful when creating using BinaryStream.open)
bin.readValue()

# Open a signed, 16bit raw image for writing and put the first value
bout = BinaryStream.open(filePath, "w", precision_in_bytes=2, signed=True)
bin.writeValue(-4)

# Read a PPM/PGM file
def readPBM(file_path):  --> (matrix, mode, precision)

# Write a PPM/PGM or RAW file:
def writePBM(output_file, matrix, precision, signed, write_header=True):
def writeRAW(output_file, matrix, precision, signed)

@author http://stackoverflow.com/questions/442188/readint-readbyte-readstring-etc-in-python
@author tuned and modified by Miguel Hernández <mhernandez@deic.uab.es>
@date 21/09/2011
"""

import struct
from struct import *

class BinaryStream:
    @staticmethod
    def open(file_path, mode, precision_in_bytes = 0, signed = False):
        """Return a binary stream from the given path and mode ("r" or "w").
        After that, the readValue() and writeValue(v) functions can be used
        for fastest access to the data."""
        stream = open(file_path, mode)
        return BinaryStream(base_stream = stream, precision_in_bytes = precision_in_bytes, signed = signed)
        
    @staticmethod
    def readPBM(file_path):
        """Read a whole PBM (either greymap PGM(P5) or color PPM(P6)) image file and 
        return:
            -a matrix with the values 
            -the selected mode "pgm"|"ppm"
            -the precision, in bytes.
        """
        
        ## Read header
        bin = BinaryStream.open(file_path, "r", precision_in_bytes=1)
        lines_read = 0
        # mode
        mode = None
        type_string = bin.readLine(include_end_char=False).strip()
        lines_read += 1
        if type_string == "P6":
            mode = "ppm"
        elif type_string == "P5":
            mode = "pgm"
        else:
            raise Exception("Only P6 (PPM) and P5 (PGM) are accepted in the file header")
        # Size
        while True:
            size_string = bin.readLine()
            lines_read += 1
            if size_string.strip()[0] == "#":
                continue
            else:
                break
        (width, height) = [int(s) for s in size_string.strip().split(" ")]
        # color count
        while True:
            count_string = bin.readLine()
            lines_read += 1
            if size_string.strip()[0] == "#":
                continue
            else:
                break
        color_count = int(count_string.strip())
        if color_count == 2**8 - 1:
            precision = 1
        elif color_count == 2**16 - 1:
            precision = 2
        else:
            raise Exception("Color count " + str(color_count) + " not valid for PBM file")
        bin.close()
        
        ## Now skip header and read data
        bin = BinaryStream.open(file_path, "r", precision_in_bytes=precision)
        for _ in range(lines_read):
            bin.readLine()
        matrix = [ [0 for y in range(height)] for x in range(width)]
        for y in range(height):
            for x in range(width):
                if mode == "ppm":
                    r = bin.readValue()
                    g = bin.readValue()
                    b = bin.readValue()
                    matrix[x][y] = (r,g,b)
                elif mode == "pgm":
                    matrix[x][y] = bin.readValue()
        
        return matrix, mode, precision
        
    @staticmethod
    def writeRAW(output_file, matrix, precision, signed):
        """Write a RAW image from the pixel matrix"""
        BinaryStream.writePBM(output_file, matrix, precision, signed, write_header=False)
        
    @staticmethod
    def writePBM(output_file, matrix, precision, signed, write_header=True):
        """Write a PBM file with the pixels in the matrix passed as arguments.
        If write_header is false, then a RAW image is created.
        - If the matrix contains integers, then a PGM file or a simple RAW file
        is created. 
        - If the matrix contains tuples, then a PPM file is created, or if RAW
        mode is selected, then three components are written (first R, then G, then B)
        """
        # Check validity
        if precision not in [1,2]:
            raise Exception("Precision " + str(precision) + " not allowed!")
        
        # Determine image mode
        if isinstance(matrix[0][0], int):
            mode = "pgm"
        elif isinstance(matrix[0][0], tuple):
            mode = "ppm"
        else:
            raise Exception("Type " + str(type(matrix[0][0])) + " not valid!")
        
        bout = BinaryStream.open(output_file, "w", precision, signed=signed)
        width = len(matrix)
        height = len(matrix[0])

        if write_header == True:
            # Write a PGM or PPM file
            if mode == "pgm":
                bout.writeString("P5\n", write_length=False)
            elif mode == "ppm":
                bout.writeString("P6\n", write_length=False)
            bout.writeString(str(width) + " " + str(height) + "\n", write_length=False)
            bout.writeString(str(2**(8*precision)-1) + "\n", write_length=False)
            if mode == "pgm":
                for y in range(height):
                    for x in  range(width):
                        bout.writeValue(matrix[x][y])
            elif mode == "ppm":
                for y in range(height):
                    for x in  range(width):
                        for band in range(3):
                            bout.writeValue(matrix[x][y][band])
        else:
            # Write a RAW file
            if mode == "pgm":
                for y in range(height):
                    for x in  range(width):
                        bout.writeValue(matrix[x][y])
            elif mode == "ppm":
                for band in range(3):
                    for y in range(height):
                        for x in  range(width):
                            bout.writeValue(matrix[x][y][band])
                            
        bout.close()
            
    def readMatrix(self, width, height):
        """Read values from the binary stream using readValue to fill
        a matrix with the specified dimensions"""
        matrix = [ [0 for y in range(height)] for x in range(width) ]
        try:
            for y in range(height):
                for x in range(width):
                    matrix[x][y] = self.readValue()
                    
            return matrix
        except struct.error:
            raise Exception("Error reading " + str(width) + "x" + str(height) + " matrix. Correct dimensions?")
        
    
    def __init__(self, base_stream, precision_in_bytes = 2, signed = False):
        """Create a new binary stream from a file stream, 
        possibly specifying precision and signedness"""
        self.base_stream = base_stream
        self.signed = signed
        self.precision_in_bytes = precision_in_bytes
        
        # Function aliases, tuned for speed
        self.readBytes = self.base_stream.read
        self.writeBytes = self.base_stream.write
        self.readValue = self.quick_unpack
        self.writeValue = self.quick_pack
        self.quick_length = precision_in_bytes
        if precision_in_bytes == 1:
            if signed == False:
                self.quick_format = "B"
            else:
                self.quick_format = "b"
        elif precision_in_bytes == 2:
            if signed == False:
                self.quick_format = ">H"
            else:
                self.quick_format = ">h"
        elif precision_in_bytes == 4:
            if signed == False:
                self.quick_format = ">I"
            else:
                self.quick_format = ">i"
        else:
            raise("Precision '"+precision_in_bytes+"' not supported")
            
    def close(self):
        self.base_stream.close()
        
        
    def readValue(self):
        """Read a value from the stream using the proper function considering
        the precision. If the precision is 0, this method will throw an exception
        
        NOTE: This function is assigned with the constructor"""
        raise Exception("NOTE: This function should be assigned with the constructor")
        
    def writeValue(self, value):
        """Write a value into the stream using the proper function considering
        the precision. If the precision is 0, this method will throw an exception
        
        NOTE: This function is assigned with the constructor"""
        raise Exception("NOTE: This function should be assigned with the constructor")
    
    def readLine(self, end_char="\n", include_end_char=False):
        """Read characters from the stream until end_char is found.
        Return them as a string, including the end_char or not
        depending on the include_end_char parameter"""
        string = []
        while True:
            c = self.readUChar()
            if c != ord(end_char[0]):
                string.append(chr(c))
            else:
                if include_end_char == True:
                    string.append(chr(c))
                return "".join(string)
        

    def readByte(self):
        return self.base_stream.read(1)

    def readBytes(self, length):
        """
        NOTE: This function is assigned with the constructor"""
        pass 

    def readChar(self):
        return self.unpack('b')

    def readUChar(self):
        return self.unpack('B')

    def readBool(self):
        return self.unpack('?')

    def readInt16(self):
        return self.unpack('>h', 2)

    def readUInt16(self):
        return self.unpack('>H', 2)

    def readInt32(self):
        return self.unpack('>i', 4)

    def readUInt32(self):
        return self.unpack('>I', 4)

    def readInt64(self):
        return self.unpack('>q', 8)

    def readUInt64(self):
        return self.unpack('>Q', 8)

    def readFloat(self):
        return self.unpack('f', 4)

    def readDouble(self):
        return self.unpack('d', 8)

    def readString(self):
        length = self.readUInt16()
        return self.unpack(str(length) + 's', length)

    def writeBytes(self, value):
        """
        NOTE: This function is assigned with the constructor"""
        pass

    def writeChar(self, value):
        self.pack('b', value)

    def writeUChar(self, value):
        self.pack('B', value)

    def writeBool(self, value):
        self.pack('?', value)

    def writeInt16(self, value):
        self.pack('>h', value)

    def writeUInt16(self, value):
        self.pack('>H', value)

    def writeInt32(self, value):
        self.pack('>i', value)

    def writeUInt32(self, value):
        self.pack('>I', value)

    def writeInt64(self, value):
        self.pack('>q', value)

    def writeUInt64(self, value):
        self.pack('>Q', value)

    def writeFloat(self, value):
        self.pack('f', value)

    def writeDouble(self, value):
        self.pack('d', value)

    def writeString(self, value, write_length=True):
        length = len(value)
        if write_length == True:
            self.writeUInt16(length)
        self.pack(str(length) + 's', value)

    def pack(self, fmt, data):
        return self.writeBytes(struct.pack(fmt, data))

    def unpack(self, fmt, length = 1):
        return struct.unpack(fmt, self.readBytes(length))[0]
        
    # Quick, aliased versions
    def quick_pack(self, data):
        return self.writeBytes(struct.pack(self.quick_format, data))
    def quick_unpack(self):
        return struct.unpack(self.quick_format, self.readBytes(self.quick_length))[0]
