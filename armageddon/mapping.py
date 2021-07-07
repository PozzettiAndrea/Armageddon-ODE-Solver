import folium
import numpy as np


def plot_circle(lat, lon, radii, map=None, **kwargs):
    """
    Plot a circle on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    radii: list
        radii of blast zones (we only need to plot the first for our example) 
    map: folium.Map 
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    >>> import folium
    >>> armageddon.plot_circle(52.79, -2.95, 1e3, map=None)
    """
    # lat = 51.5073219
    # lon = -0.1276474

    lat_c = 53.2744122
    lon_c = -9.0490632

    if not map:
    #     map = folium.Map(location=[lat, lon], control_scale=True)
    # folium.Circle([lat, lon], radius, fill=True, fillOpacity=0.6, **kwargs).add_to(map)
        map_UK = folium.Map(location=[lat_c, lon_c], zoom_start=5.5, control_scale=True)

    folium.Circle(
        [lat, lon], 
        radius = radii[0], 
        fill=True, 
        fillOpacity=0.1, 
        fill_color='blue',
        color='blue',parse_html=False).add_to(map_UK)

    map_UK.save('plot_data.html')

    return map_UK

def plot_line(lat, lon, angle, distance, map_UK=None, **kwargs):

    lat_c = 53.2744122
    lon_c = -9.0490632

    if not map_UK:
    #     map = folium.Map(location=[lat, lon], control_scale=True)
    # folium.Circle([lat, lon], radius, fill=True, fillOpacity=0.6, **kwargs).add_to(map)
        map_UK = folium.Map(location=[lat_c, lon_c], zoom_start=5.5, control_scale=True)
    
    dx_lat = distance*np.sin(20)/11132
    dx_lon = distance*np.cos(20)/10000

    ls = folium.PolyLine(locations=[[lat,lon],[lat + dx_lat, lon - dx_lon]],
                    color='black')

    ls.add_to(map_UK)
    map_UK.save('plot_data.html')
    return map_UK

def plot_outcome(lat, lon, radii, bearing, distance, map=None, **kwargs):
    """
    Plot a circle on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    radii: list
        radii of blast zones
    bearing: float
        bearing of trajectory
    distance: float
        distance traversed by meteor before bursting
    map: folium.Map 
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    >>> import folium
    >>> armageddon.plot_circle(52.79, -2.95, 1e3, map=None)
    """
    # lat = 51.5073219
    # lon = -0.1276474

    lat_c = 53.2744122
    lon_c = -9.0490632

    if not map:
    #     map = folium.Map(location=[lat, lon], control_scale=True)
    # folium.Circle([lat, lon], radius, fill=True, fillOpacity=0.6, **kwargs).add_to(map)
        map_UK = folium.Map(location=[lat_c, lon_c], zoom_start=5.5, control_scale=True)

    folium.Circle(
        [lat, lon], 
        radius = radii[3], 
        fill=False, 
        fillOpacity=1.5, 
        weight = 2,
        popup=folium.Popup(max_width=450),
        fill_color='red',
        color='red',parse_html=False).add_to(map_UK)

    folium.Circle(
        [lat, lon], 
        radius = radii[2], 
        fill=False, 
        fillOpacity=1, 
        weight = 2,
        popup=folium.Popup(max_width=450),
        fill_color='yellow',
        color='yellow',parse_html=False).add_to(map_UK)

    folium.Circle(
        [lat, lon], 
        radius = radii[1], 
        fill=True, 
        fillOpacity=0.2, 
        fill_color='green',
        color='green',parse_html=False).add_to(map_UK)

    folium.Circle(
        [lat, lon], 
        radius = radii[0], 
        fill=True, 
        fillOpacity=0.1, 
        fill_color='blue',
        color='blue',parse_html=False).add_to(map_UK)
    
    folium.Marker([lat, lon],
              popup='Place',
              clustered_marker = True,
              color = 'yellow',
              icon=folium.Icon(color = 'beige')).add_to(map_UK)
    
    ls = folium.PolyLine(locations=[[lat,lon],[lat+2,lon-2]],
                    color='black')

    ls.add_to(map_UK)

    map_UK.save('plot_data.html')

    return map_UK


def multiple_plots(latlon_array, map=None, **kwargs):
    """
    Plot many circle on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    latlon_array: list
        list of lists, each containing latitude, longitude and radius of our randomised blasts.
    map: folium.Map 
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    >>> import folium
    >>> armageddon.plot_circle([[52.79, -2.95, 10e3],[52.78, -2.94, 12e3],[52.785, -2.943, 11e3]], map=None)
    """
    # lat = 51.5073219
    # lon = -0.1276474

    lat_c = 53.2744122
    lon_c = -9.0490632

    if not map:
    #     map = folium.Map(location=[lat, lon], control_scale=True)
    # folium.Circle([lat, lon], radius, fill=True, fillOpacity=0.6, **kwargs).add_to(map)
        map_UK = folium.Map(location=[lat_c, lon_c], zoom_start=5.5, control_scale=True)
    length = len(latlon_array)
    for lat, lon, radius in latlon_array:

        folium.Circle(
            [lat, lon], 
            radius = radius, 
            fill=True, 
            fillOpacity=0.1, 
            fill_color='blue',
            color='blue',parse_html=False).add_to(map_UK)
        
        # ls = folium.PolyLine(locations=[[lat,lon],[lat+2,lon-2]],
        #                 color='black')

        # ls.add_to(map_UK)

    map_UK.save('plot_multiple.html')

    return map_UK
