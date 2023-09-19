#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:25:00 2023

@authors:    
Ryan Boukrouche (RB)
"""
import numpy as np

# Specify the file name
file_name = "/proj/bolinc/users/x_ryabo/spectral_files/Reach/Reach"

x = None
y = None

# Read the file and process lines
with open(file_name, 'r') as file:
    lines = file.readlines()
    
    # Initialize flag to indicate if we are inside BLOCK 1
    inside_BLOCK_1 = False

    for i, line in enumerate(lines):
        # Check if the line contains the section header
        if "Limits of spectral intervals" in line:
            inside_BLOCK_1 = True
        elif inside_BLOCK_1:
            # Check if the line contains "Lower limit"
            if "Lower limit" in line:
                x = i + 1  # Set x to the index of the line following "Lower limit"
            elif "*END" in line:
                if x is not None:
                    y = i  # Set y to the index of the next occurrence of "*END" after "Lower limit"
                    break  # Exit the loop

# Check if x and y were found and extract the desired lines
if x is not None and y is not None:
    relevant_lines = lines[x:y]

    # Extract the lower limits from the relevant lines
    lower_limits = []
    for line in relevant_lines:
        line_parts = line.strip().split()
        lower_limit = float(line_parts[1])
        # Add the upper limit of the first line
        if line_parts[0] == '1':
            upper_limit = float(line_parts[2])

        lower_limits.append(lower_limit)

    bands = np.array(lower_limits)
    # Insert the upper limit 
    bands = np.append(upper_limit, bands)
    # Print the resulting array
    print(bands)
else:
    print("Could not find 'Lower limit' and '*END' within the specified section.")