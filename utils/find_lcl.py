import numpy as np

def find_intersection(curve1, curve2, tolerance=1e-1):
    # Starting from the TOA, find the last index where the absolute difference is within the tolerance
    intersection_indices = np.where(np.abs(curve1 - curve2) < tolerance)[-1]
    if len(intersection_indices) > 0:
        # Find the first occurrence of intersection
        first_intersection_index = intersection_indices[-1]
        first_intersection_value = curve1[first_intersection_index]
        
        #print("LCL index:", first_intersection_index)
        #print("LCL value:", first_intersection_value)
        return first_intersection_index, first_intersection_value
    else:
        print("No LCL found")
        return None