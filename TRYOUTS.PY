import rasterio

with rasterio.open(r'DTM Data/plain/plain_terrrain_project_dtm.tif') as dem:
    data_set = dem.read(1)
cell_size = abs(dem.bounds[0]-dem.bounds[2])/dem.width

print(cell_size)
print(dem.meta["transform"][0])