import glob
import os

import geopandas as gpd
import pandas as pd

#######
# Vietnam Boundaries
#######

sf = gpd.read_file("./vnm_adm_gov_20201027/vnm_admbnda_adm1_gov_20201027.shp")
print(sf.crs)  # Check EPSG / CRS -- epsg:4326 = WGS 84


sf = sf.to_crs(epsg=3405)  # Convert to WGS 84


sf.columns = [x.lower() for x in sf.columns]  # lower case col names
print(sf.dtypes)

# Check geometries
print(len(sf))
sf.geometry.type.value_counts()

sf.head()

viet_provinces = sf[["adm1_en", "adm1_vi", "shape_leng", "shape_area", "geometry"]]


#######
# Pollution data
# Link for data (just tweak dates and vars):
# Data Source: https://giovanni.gsfc.nasa.gov/giovanni/#service=TmAvMp&starttime=&endtime=&bbox=102,8,110.25,24.25&variableFacets=dataFieldMeasurement%3ANO2%2CSO2%3B
#
# Notes:
#
# the unit for the data I downloaded is Dobson unit.
# for NO2 data, there are two measures available. One is NO2 Tropospheric Column (30% Cloud Screened), the other is NO2 Total Column (30% Cloud Screened). I am not sure which one you wanna go with so I kept both of them in the final data
# for SO2 data, there is only one kind available that is SO2 Column Amount – so there’s no cloud screened
#######

import earthpy as et

# Load and validate
import rioxarray as rxr

YEAR = 2024

dirs = dict()
prefix = f"./{YEAR}/GIOVANNI-g4.timeAvgMap."
suffix = f".{YEAR}0101-{YEAR}1231.102E_8N_110E_24N.tif"


dirs["NO2"] = prefix + "OMNO2d_003_ColumnAmountNO2CloudScreened" + suffix
dirs["NO2_trop"] = prefix + "OMNO2d_003_ColumnAmountNO2TropCloudScreened" + suffix
dirs["SO2"] = prefix + "OMSO2e_003_ColumnAmountSO2" + suffix


# View generate metadata associated with the raster file
NO2 = rxr.open_rasterio(dirs["NO2"], masked=True)
print("The bounds of your data are:", NO2.rio.bounds())
print("The crs of your data is:", NO2.rio.crs)
print("The shape of your data is:", NO2.shape)
print("The spatial resolution for your data is:", NO2.rio.resolution())
print("The metadata for your data is:", NO2.attrs)


# View generate metadata associated with the raster file
NO2_trop = rxr.open_rasterio(dirs["NO2_trop"], masked=True)
print("The bounds of your data are:", NO2_trop.rio.bounds())
print("The crs of your data is:", NO2_trop.rio.crs)
print("The shape of your data is:", NO2_trop.shape)
print("The spatial resolution for your data is:", NO2_trop.rio.resolution())
print("The metadata for your data is:", NO2_trop.attrs)


# View generate metadata associated with the raster file
SO2 = rxr.open_rasterio(dirs["SO2"], masked=True)
print("The bounds of your data are:", SO2.rio.bounds())
print("The crs of your data is:", SO2.rio.crs)
print("The shape of your data is:", SO2.shape)
print("The spatial resolution for your data is:", SO2.rio.resolution())
print("The metadata for your data is:", SO2.attrs)

import rasterio
from rasterio.plot import show

# Load the GeoTIFF file of NO2 observations
# raster_path = 'path_to_your_geotiff/NO2_observations.tif'
NO2 = rasterio.open(dirs["NO2"])

# Optionally, visualize the raster data
show(NO2)


import rasterio
from rasterio.warp import Resampling, calculate_default_transform, reproject


def reprojection(input_raster, output_raster, new_crs):
    with rasterio.open(input_raster) as src:
        transform, width, height = calculate_default_transform(
            src.crs, new_crs, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update(
            {"crs": new_crs, "transform": transform, "width": width, "height": height}
        )

        with rasterio.open(output_raster, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=new_crs,
                    resampling=Resampling.nearest,
                )


for file in glob.glob(f"./{YEAR}/*.tif"):
    reprojection(file, f"./{YEAR}_reprojected/" + file.split("/")[-1], "EPSG:3405")


for file in glob.glob(f"./{YEAR}_reprojected/*.tif"):
    with rasterio.open(file) as src:
        assert src.crs == rasterio.crs.CRS.from_epsg(3405)


print(sf.crs)

import re

for d in dirs.keys():
    dirs[d] = re.sub(f"/{YEAR}/", f"/{YEAR}_reprojected/", dirs[d])
    print(dirs[d])


from rasterstats import zonal_stats

# Use zonal_stats to calculate the average NO2 density for each province
mean_no2 = zonal_stats(viet_provinces, dirs["NO2"], stats="mean")
mean_no2_trop = zonal_stats(viet_provinces, dirs["NO2_trop"], stats="mean")
mean_so2 = zonal_stats(viet_provinces, dirs["SO2"], stats="mean")

# Add the average NO2 density to the provinces GeoDataFrame
viet_provinces["NO2_mean"] = [stat["mean"] for stat in mean_no2]
viet_provinces["NO2_trop_mean"] = [stat["mean"] for stat in mean_no2_trop]
viet_provinces["SO2_mean"] = [stat["mean"] for stat in mean_so2]


# Check the first few rows to verify the NO2 density values
viet_provinces.head()

# All measures are in Dobson unit:
#
# The Dobson unit (DU) is a unit of measurement of the
# amount of a trace gas in a vertical column through the Earth's atmosphere.

viet_provinces_nogeom = viet_provinces.drop(columns=["geometry"])

viet_provinces_nogeom.to_csv(f"vietnam_air_pollution_{YEAR}.csv", index=False)

viet_provinces_reprojected = viet_provinces.to_crs(epsg=4326)


import geopandas as gpd

# plot the NO2 mean values for 2022 and 2023 on a map of Vietnam side by side
import matplotlib.pyplot as plt

for p in ["NO2", "NO2_trop", "SO2"]:
    fig, ax = plt.subplots(figsize=(10, 10))
    viet_provinces_reprojected.plot(
        column=f"{p}_mean",
        cmap="OrRd",
        linewidth=0.8,
        ax=ax,
        edgecolor="0.8",
        legend=True,
    )
    ax.axis("off")
    ax.set_title(f"Mean {p} concentration (DU) in {YEAR}")

    fig.savefig(f"{p}_{YEAR}.pdf", format="pdf")
