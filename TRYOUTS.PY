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
    return math.atan(max(slopes))

def get_slope_ConstrainedQuadSurface( z,cellsize ):
    # define Z just transpose of the normalZ
    # caculate X = ( ( inverse_of( A_transpose*A ) ) * A_transpose ) * Z
    # slopeX = X at [0,3] 
    # slopeY = X at [0,4] 
    g = cellsize # for simplicity in repitition
    matA = np.array([
        ( g**2,  g**2, -g**2,   -g,    g,    1 ), #1
        (    0,  g**2,     0,    0,    g,    1 ), #2
        ( g**2,  g**2,  g**2,    g,    g,    1 ), #3 
        ( g**2,     0,     0,   -g,    0,    1 ), #4
        (    0,     0,     0,    0,    0,    1 ), #5
        ( g**2,     0,     0,    g,    0,    1 ), #6
        ( g**2,  g**2,  g**2,   -g,   -g,    1 ), #7
        (    0,  g**2,     0,    0,   -g,    1 ), #8
        ( g**2,  g**2,  g**2,    g,   -g,    1 ), #9
    ])
    matZ = np.array([z]).T
    matA_transpose = matA.T
    # matX is the matrice of coefficients
    # caculate X = ( ( inverse_of( A_transpose*A ) ) * A_transpose ) * Z in two stepsd
    temp = np.linalg.inv( np.dot( matA_transpose,matA ) )
    matX = np.dot( np.dot( temp,matA_transpose ),matZ )
    slope_x = matX[3,0]
    slope_y = matX[4,0]
    return math.atan( math.sqrt(slope_x**2 + slope_y**2) )


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
    [(1,2,1,1), 
    (1,2,2,1), 
    (1,2,1,0),
    (1,1,0,0)
    ])

slope_MaximumMax = np.zeros(( 2, 2 ))
set_slope_dataset_with( 4, 4, 1, dem, slope_MaximumMax, get_slope_MaximumMax )

slope_quadSurf = np.zeros(( 2, 2 ))
set_slope_dataset_with( 4, 4, 1, dem, slope_quadSurf, get_slope_ConstrainedQuadSurface )

print( slope_MaximumMax )
print( slope_quadSurf )
