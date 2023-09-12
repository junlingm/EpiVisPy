import pandas as pd
import geopandas as gpd
import pyreadr
import matplotlib.pyplot as plt
import numpy as np

from shapely.geometry import Polygon, MultiPolygon
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def generate_animation(population_path, bc_map_path, time_series_path, max_height=130000,frame_rate=200,cylinder_radius=5,camera_position=(20,-90),map_color_func=None,color_by_case=False,camera_func=None, map_fill_color='YlOrRd', cylinder_color='YlOrRd'):

    #Section 1
    #====================================================
    
    #1) Reading and Processing Population Data
    population = pd.read_excel(population_path)
    population.columns = ["English_name", "Population"]
    
    #2) Reading Geospatial Data and Assigning ID
    BCmap = gpd.read_file(bc_map_path)
    BCmap['new_ID_edited'] = range(len(BCmap))

    
    #3) Coordinate System Adjustments
    #Based on the latitude and longitude, it determines the UTM zone and constructs a new EPSG code.
    #converts the BCmap to this UTM coordinate system using the to_crs function if it's not UTM system.
    if not BCmap.crs.to_string().startswith('EPSG:326') and not BCmap.crs.to_string().startswith('EPSG:327'):
        lon, lat = BCmap.geometry.unary_union.centroid.xy
        zone_number, zone_letter = utm.from_latlon(lat[0], lon[0])[:2]
        epsg_code = 32600 + zone_number if lat[0] >= 0 else 32700 + zone_number
        BCmap = BCmap.to_crs(epsg=f"EPSG:{epsg_code}")
        
        
    
    #4) Creating map date frame (merging)
    map_data = pd.DataFrame({
        'id': BCmap['new_ID_edited'],
        'Code': BCmap['ABRVN'],
        'English_name': BCmap['AA_NAME']
    })
    map_data['Code'] = map_data['Code'].str.replace('BC.', '')
    map_data = map_data.merge(population, on='English_name', how='left')
    map_data['Population'] = map_data['Population'].astype(float)

    
    #5) Area Calculation and Population Density Determination,then create another data frame
    BCmap['Area'] = BCmap.area
    region_data = map_data[['id', 'English_name']].drop_duplicates()
    region_data = region_data.merge(population, on='English_name', how='left')
    region_data['Population'] = region_data['Population'].astype(float)

    area_df = pd.DataFrame({'id': BCmap['new_ID_edited'].unique(), 'Area': BCmap['Area']})
    region_data = region_data.merge(area_df, on='id', how='left')
    region_data['PopulationDensity'] = region_data['Population'] / region_data['Area']
    #============================================================
    
    
    
    #Section 2
    #===================================================================================================
    
    # 1) Map color
    # The resultant colors are stored in the colors and cylinder_colors variables, respectively.
    colors = plt.get_cmap(map_fill_color)
    cylinder_colors = plt.get_cmap(cylinder_color)
    
    
    #2) Custom Color Values Calculation
    # if a custom color function (map_color_func) is provided, returned color values generate by the map_color_func  are stored in the color_values variable
    # all values in the color_values list are between 0 and 1.
    if map_color_func is not None:
        color_values = custom_color_func()
        if any([v < 0 or v > 1 for v in color_values]):
            
            raise ValueError("All values returned by color_func should be between 0 and 1.")
        
    #3) Default Color Values Calculation    
    else:
        color_values = region_data['PopulationDensity']
        
    normed_densities= (color_values - color_values.min()) / (color_values.max() - color_values.min())
    #=====================================================================================================  

   

    #Section 3
    #=====================================================================================================
    
    #1)  loads a CSV file
    time_series_data = pd.read_csv(time_series_path)
    
    #2) Identifying Case Columns
    # the dataset has multiple columns indicating the number of cases for different regions.
    case_columns = [col for col in time_series_data.columns if "_Cases_Number" in col]
    
    #3) Extracting Region Names
    # extracts the prefix before the first underscore (_) and region's name associated with the cases.
    regions = [col.split('_')[0] for col in case_columns]
    
    
    #4) Normalizing Case Numbers
    for region in regions:
        col_name = f"{region}_Cases_Number"
        new_col = f"Normalized_{col_name}"
        time_series_data[new_col] = time_series_data[col_name] / time_series_data[col_name].max() * max_height

    
    #5) Plot Initialization
    
    # A figure of size 15x15 is created.
    fig = plt.figure(figsize=(15, 15))
    
    #A 3D subplot is added to the figure.
    ax = fig.add_subplot(111, projection='3d')
    
    # The initial camera/viewing position of the 3D plot
    ax.view_init(elev=camera_position[0], azim=camera_position[1])
    
    # A boundary_linewidth variable is initialized with a value of 0.5
    boundary_linewidth = 0.5
    #=======================================================================================================================
    
    
    
    #Section 4
    #========================================================================================================================
    def update(num, time_series_data, ax):
        
        #1) Clearing the existing plots on the 3D axis before plotting the new data.
        ax.cla()
        
        #2) Helper Function for Plotting Polygon Boundaries
        # Function plots the boundaries of a polygon poly on the 3D axis at a specific height z_value.
        def plot_polygon_boundary(poly, z_value, ax, color='black', linewidth=boundary_linewidth):
            x, y = poly.exterior.xy
            ax.plot(x, y, zs=z_value, zdir='z', color=color, linewidth=linewidth)
            for interior in poly.interiors:
                x, y = interior.xy
                ax.plot(x, y, zs=z_value, zdir='z', color=color, linewidth=linewidth)
        
        
        #3) Plotting Map Geometries
        #This loop iterates over rows in BCmap and plots the geometry of each region
        
        for idx, map_row in BCmap.iterrows():
            color = plt.get_cmap(map_fill_color)(normed_densities[idx])
            
            # If the geometry is a single polygon ('Polygon'), the polygon is added to the 3D plot with its color determined by population density or color function.
            if isinstance(map_row['geometry'], Polygon):
                
                poly = map_row['geometry']
                
                verts = [list(zip(*poly.exterior.xy, np.zeros_like(poly.exterior.xy[0])))]
                poly3d_fill = Poly3DCollection(verts, facecolor=color, edgecolor='black', linewidth=boundary_linewidth)
                ax.add_collection3d(poly3d_fill)
            
            
            # If the geometry consists of multiple polygons (MultiPolygon), each individual polygon is plotted.
            elif isinstance(map_row['geometry'], MultiPolygon):
                for poly in map_row['geometry'].geoms:
                    
                    verts = [list(zip(*poly.exterior.xy, np.zeros_like(poly.exterior.xy[0])))]
                    poly3d_fill = Poly3DCollection(verts, facecolor=color, edgecolor='black', linewidth=boundary_linewidth)
                    ax.add_collection3d(poly3d_fill)
        
        #4) Setting Axis Limits
        ax.set_xlim(BCmap.total_bounds[0], BCmap.total_bounds[2])
        ax.set_ylim(BCmap.total_bounds[1], BCmap.total_bounds[3])
        ax.set_zlim(0, max_height)
        
        row = time_series_data.iloc[num]
        
        
        #5) Plotting Cylinders Representing Cases Number
        # This loop visualizes the number of cases for each region as cylinders
        for region in regions:
            
            # For each region, it retrieves the UTM coordinates and normalized case numbers.
            easting_col = f"{region}_UTM_Easting"
            northing_col = f"{region}_UTM_Northing"
            cases_col = f"Normalized_{region}_Cases_Number"
            
            
            # If color_by_case is set to True, the color of the cylinder is determined by the number of cases; otherwise, it defaults to purple.
            if color_by_case is True:
                norm_case = row[cases_col] / max_height
                color =  cylinder_colors(norm_case)
                
            else:
                color = 'purple'
            
            
            # The height of the cylinder corresponds to the normalized cases number.
            segments = [((row[easting_col], row[northing_col], 0), (row[easting_col], row[northing_col], row[cases_col]))]
            radius = cylinder_radius
            lc = Line3DCollection(segments, linewidths=radius, colors=color, capstyle='round')
            ax.add_collection(lc)
            
            
            ax.set_xlim(BCmap.total_bounds[0], BCmap.total_bounds[2])
            ax.set_ylim(BCmap.total_bounds[1], BCmap.total_bounds[3])
            ax.set_zlim(0, max_height)
        
        
        #6)  Setting the Camera Position
        # This code determines the camera viewpoint of the 3D plot
        
        # If a camera_func is provided, it calculates the camera position based on the current time step 'num' (number of frame).
        if camera_func is not None:
            elev, azim = camera_func(num)
            ax.view_init(elev=elev, azim=azim)
            
        # Otherwise, it defaults to a pre-defined 'camera_position'.
        else:
            ax.view_init(elev=camera_position[0], azim=camera_position[1])
            
    #========================================================================================================================
    
    ani = FuncAnimation(fig, update, frames=len(time_series_data), fargs=[time_series_data, ax],interval=frame_rate)
    return ani




#Exampe 1
# Example of custom camera position based on frame number
#This function calculates and returns the camera's elevation and azimuthal angle based on the provided frame number. 
#This allows for dynamic camera positioning as the frames of a temporal visualization progress.

def custom_camera_position(frame_num):
    
    # For frame numbers between 10 and 20 inclusive, the camera's elevation angle is set to 30 degrees.
    # The azimuthal angle starts at -80° for frame 10 (since -90 + 10 = -80) and increases by 1° for each subsequent frame until it reaches -70° for frame 20.
    if 10 <= frame_num <= 20:
        
        return (30, -90 + frame_num)
    
    
    # For frame numbers between 30 and 40 inclusive, the camera's elevation angle is set to 50 degrees.
    # The azimuthal angle starts at -120° for frame 30 (since -90 - 30 = -120) and decreases by 1° for each subsequent frame, ending at -130° for frame 40.
    elif 30 <= frame_num <= 40:
        
        return (50, -90 - frame_num)
    
    
    #For any frame number outside the ranges of 10-20 and 30-40, the camera's elevation is set to 20 degrees, and the azimuthal angle is set to -90°.
    else:
        
        return (20, -90)

ani = generate_animation(
    'C:/Users/Dell/Desktop/STAT457/Population_data.xlsx',
    'C:/Users/Dell/Desktop/STAT457/bc_map_utm.shp',
    'C:/Users/Dell/Desktop/STAT457/Time_series_data.csv',
    frame_rate=200,
    cylinder_radius=5,
    camera_position=(20, -90),color_by_case=True,
    camera_func=custom_camera_position,
    cylinder_color='viridis'
)
    
HTML(ani.to_html5_video())





#Example 2
#Example of custom color function
#==========================================================================
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt


population = pd.read_excel('C:/Users/Dell/Desktop/STAT457/Population_data.xlsx')


population.columns = ["English_name", "Population"]

BCmap = gpd.read_file("C:/Users/Dell/Desktop/STAT457/bc_map_utm.shp")

map_data = pd.DataFrame({
    'id': BCmap['new_ID1'],
    'Code': BCmap['ABRVN'],
    'English_name': BCmap['AA_NAME']
})

map_data['Code'] = map_data['Code'].str.replace('BC.', '')



map_data = map_data.merge(population, on='English_name', how='left')

map_data['Population'] = map_data['Population'].astype(float)


boundaries = []
for poly in BCmap.geometry:
    if isinstance(poly, Polygon):
        boundaries.append(list(poly.exterior.coords))
    elif isinstance(poly, MultiPolygon):
        multi_coords = []
        for sub_poly in poly.geoms:
            multi_coords.append(list(sub_poly.exterior.coords))
        boundaries.append(multi_coords)

BCmap['Area'] = BCmap.area

region_data = map_data[['id', 'English_name']].drop_duplicates()

region_data = region_data.merge(population, on='English_name', how='left')

region_data['Population'] = region_data['Population'].astype(float)


#The entire expression essentially normalizes the population values.
#Use normalization to change the range of values in a dataset to be between 0 and 1.
#When visualizing data using color gradients, it can ensure the lowest value gets one end of the color spectrum and the highest value gets the other end.
def custom_color_func():
    values = (region_data['Population'] - region_data['Population'].min()) / \
             (region_data['Population'].max() - region_data['Population'].min())
    return values

#===================================================================================================

def custom_camera_position(frame_num):
    if 1 <= frame_num <= 20:
        
        return (30, -90 + frame_num)
    elif 20 <= frame_num <=45:
        
        return (60, -90 - frame_num)
    else:
        
        return (20, -90)

ani = generate_animation(
    'C:/Users/Dell/Desktop/STAT457/Population_data.xlsx',
    'C:/Users/Dell/Desktop/STAT457/bc_map_utm.shp',
    'C:/Users/Dell/Desktop/STAT457/Time_series_data.csv',
    frame_rate=400,
    cylinder_radius=10,
    camera_position=(30, -80),map_color_func=custom_color_func,color_by_case=True,camera_func=custom_camera_position,map_fill_color='viridis',cylinder_color='YlOrRd')

HTML(ani.to_html5_video())





#Introductions of arguments in the 'generate_animation' function

#1. population_path: Excel file that has the population data
#   bc_map_path:  shapefile
#   time_series_path are: CSV file that has the time series data


#2. max_height
#1) maximum number of cases (what it is?)
#2) Provide a numerical value to set the height limit in the 3D visualization. The default value is 130000.(How to use?)


#3. frame_rate
#1) The speed at which the animation progresses. This represents the number of frames shown per second in the animation.
#2) Set a numerical value to determine the animation's speed. A higher value will make the animation run faster. The default is 200.


#4. cylinder_radius
#1) This specifies the radius of the cylinders in the 3D plot.
#2) Provide a numerical value for the desired radius of the cylinders. The default is 5.


#5. camera_position
#1) A tuple that determines the camera's elevation and azimuthal angle for viewing the 3D visualization.
#2) Provide a tuple, where the first value is the elevation and the second is the azimuthal angle. Default is (20, -90).


#6. map_color_func
#1) An optional function that customizes the color of regions based on certain conditions or data values
#2) If you have a custom function to derive color values (like custom_color_func), you can provide it here. If not provided, the default values will be used.


#7. color_by_case
#1) A boolean flag that determines if the color of cylinders in the 3D plot should be based on case numbers.
#2) Set it to True if you want the color of cylinders to vary based on case numbers, otherwise False. The default is False.


#8. camera_func
#1) An optional function to dynamically change the camera's position during the animation.
#2)  If you have a custom function to adjust camera angles over time (like custom_camera_position), you can provide it here.


#9. map_fill_color
#1) The colormap to be used for filling the regions on the map.
#2)  Provide the name of a matplotlib colormap as a string, such as 'YlOrRd'. The default is 'YlOrRd'.


#10. cylinder_color
#1) The colormap to be used for coloring the cylinders in the 3D plot.
#2)  Provide the name of a matplotlib colormap as a string. The default is 'YlOrRd'.

