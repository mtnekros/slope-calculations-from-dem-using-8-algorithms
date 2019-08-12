import math
import numpy as np 


# file suffix according to algrithms
def get_slope_MaximumMax(z, cell_size):
    # calculating slope in N,E,W,S direction from center
    indices = [1, 3, 5, 7]
    slopes = [(abs(z[4] - z[i])) / float(cell_size) for i in indices]
    # calculating the slope sin NE, NW, SE, SW directions
    indices = [0, 2, 6, 8]
    slopes += [abs(z[4] - z[i]) / float(cell_size * math.sqrt(2)) for i in indices]
    max_slopes = max( slopes )
    return math.atan(max(slopes))


def get_flattened_moving_window(dem, i, j):
    # flattening and reversing the moving window for convenience
    return dem[i-1:i+2, j-1:j+2].flatten()[::-1]


def set_slope_dataset_with( nRows, nCols, cell_size, dem_dataset, slope_dataset, getSlope ):
    # nRows and nCols are of dem_dataset not of slope_dataset
    for i in range(nRows-2):
        for j in range(nCols-2):
            # extracting the moving window values
            window_z = get_flattened_moving_window(dem_dataset, i+1, j+1)
            slope_dataset[i, j] = math.degrees( getSlope(window_z, cell_size) ) # changed degrees alter if needed later


dem = np.array(
    [(1,1,1,1), 
    (1,2,1,1), 
    (34,1,1,1),
    (13,14,15,16)
    ])

slope = np.zeros(( 2, 2 ))


set_slope_dataset_with( 4, 4, 1, dem, slope, get_slope_MaximumMax )

print( slope )