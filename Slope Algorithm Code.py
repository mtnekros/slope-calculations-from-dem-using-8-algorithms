import numpy as np
import math
import rasterio
import matplotlib.pyplot as plt


def get_slope_2FD(z, cell_size):
    slope_x = (z[5]-z[3])/float(2*cell_size)
    slope_y = (z[7]-z[1])/float(2*cell_size)
    return math.atan(math.sqrt(slope_x**2 + slope_y**2))


def get_slope_MaximumMax(z, cell_size):
    # calculating slope in N,E,W,S direction from center
    indices = [1, 3, 5, 7]
    slopes = [(abs(z[4] - z[i])) / float(cell_size) for i in indices]
    # calculating the slope sin NE, NW, SE, SW directions
    indices = [0, 2, 6, 8]
    slopes += [abs(z[4] - z[i]) / float(cell_size * math.sqrt(2)) for i in indices]
    return math.atan(max(slopes))


def get_slope_SimpleD(z, cell_size):
    slope_x = (z[4]-z[3])/float(cell_size)
    slope_y = (z[4]-z[1])/float(cell_size)
    return math.atan(math.sqrt(slope_x**2 + slope_y**2))


def get_slope_AvgNeighbourhood(z, cell_size):
    slope_x = (z[2] - z[0] + 2 * (z[5] - z[3]) + z[8] - z[6]) / float(8*cell_size)
    slope_y = (z[6] - z[0] + 2 * (z[7] - z[1]) + z[8] - z[2]) / float(8*cell_size)
    return math.atan(math.sqrt(slope_x**2 + slope_y**2))


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


def saveDatasetToFile( dataset,filename ):
    np.savetxt(f"{filename}.csv", dataset, delimiter=",")

def main():
    # opening the dem file with rasterio
    # input for path + filename in one variable # r'C:/Users/Diwas/Desktop/Studies/8th Semester/Final Year Project/Data/Clipped Dem/steeptry1'
    input_file = input("Enter the filepath + name: ")
    with rasterio.open( input_file ) as dem:
        data_set = dem.read(1)
    cell_size = dem.meta["transform"][0] # the first index of transform is cell_size x and 4th index is cell_size y
    print("File opened!")

    # initializing the slope matrix with np array of zeroes
    nRows, nCols = dem.shape
    slope_data_set = [ np.zeros((nRows-2, nCols-2)) for i in range(4) ]

    #  names and the methods of algorithms ordered in the same way for setting slope_data_set and saving files later
    algorithmNames = [ "SimpleD", "AvgNeighbourhood", "MaximumMax", "2FD" ]
    get_slope_methods = [ get_slope_SimpleD, get_slope_AvgNeighbourhood, get_slope_MaximumMax, get_slope_2FD ]
    
    # calculating slope in the predefined order
    for i, get_slope in enumerate( get_slope_methods ):
        set_slope_dataset_with( nRows, nCols, cell_size, data_set, slope_data_set[i], get_slope )
        print( "slope calculation done for", algorithmNames[i] )

    # extracting metadata as kwargs and updating it
    kwargs = dem.meta.copy()
    # getting the bounds
    left, bottom, right, top = dem.bounds
    transform = rasterio.transform.from_bounds(
        left + cell_size,
        bottom + cell_size,
        right - cell_size,
        top - cell_size,
        nCols - 2,
        nRows - 2)

    kwargs.update(
        transform=transform,
        height=nRows-2,
        width=nCols-2,
        driver='GTiff'
    )
    
    # writing slope datasets to csv files
    for i, name in enumerate( algorithmNames ):
        saveDatasetToFile( slope_data_set[i],name ) # we can use index becoz slope data set and algorithm names index matches

    # writing to all files with corresponding algorithm name in order
    output_filename_suffix = input("Enter the output filename suffix: ")
    for i, name in enumerate( algorithmNames ):
        with rasterio.open(f'DTM Data/steep/try1/{output_filename_suffix}_{name}.geotiff', 'w', **kwargs) as out_file:
            out_file.write(slope_data_set[i].astype(rasterio.float32), 1)
        print(name, "File saved!")


if __name__ == '__main__':
    main()





















# import rasterio
# import math
# import matplotlib.pyplot as plt


# def get_slope(dem, cellsize, j, i): # average neighbour
#     k = []
#     for row in dem[j-1:j+2]:
#         for el in row[i-1:i+2]:
#             k.append(el)
#     assert(len(k) == 9)
#     slope_x = k[0] + 2*k[3] + k[6] - (k[2] + 2*k[5] + k[8])/(8*cellsize)
#     slope_y = k[0] + 2*k[1] + k[2] - (k[6] + 2*k[7] + k[8])/(8*cellsize)
#     return math.atan(math.sqrt(slope_x**2 + slope_y**2))
#
#
# with rasterio.open(r'data/dem') as dem:
#     data_set = dem.read(1)
# cell_s = abs(dem.bounds[0]-dem.bounds[2])/dem.width
#
# slope_data_set = []
#
# slope = get_slope
# for j, row in enumerate(data_set[1:-1]):
#     temp_row = []
#     for i, el in enumerate(row[1:-1]):
#         temp_row.append(slope(data_set, cell_s, j+1, i+1))
#     slope_data_set.append(temp_row)
#
# # setting pyplot axis scale
# plt.axis('equal')
#
# # showing the plot
# plt.imshow(slope_data_set)
# plt.show()
#
# # extracting metadata as kwargs and updating it
# kwargs = dem.meta
#
# # writing to file
# with rasterio.open('data/slope_nearestNeighbourhood','w', **kwargs) as out_file:
#     out_file.write(slope_data_set, 1)
