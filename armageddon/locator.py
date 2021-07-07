"""Module dealing with postcode information."""

import os
import sys
import time
import numpy as np
import pandas as pd

# Worth trying later
# project_dir = os.path.dirname(__file__)
postcodes_path = os.path.join('./armageddon/resources/full_postcodes.csv')
population_path = os.path.join('./armageddon/resources/population_by_postcode_sector.csv')

def great_circle_distance(latlon1, latlon2):
    """
    Calculate the great circle distance (in metres) between pairs of 
    points specified as latitude and longitude on a spherical Earth
    (with radius 6371 km).

    Parameters
    ----------

    latlon1: arraylike
        latitudes and longitudes of first point (as [n, 2] array for n points)
    latlon2: arraylike
        latitudes and longitudes of second point (as [m, 2] array for m points)

    Returns
    -------

    numpy.ndarray
        Distance in metres between each pair of points (as an n x m array)

    Examples
    --------

    >>> import numpy
    >>> fmt = lambda x: numpy.format_float_scientific(x, precision=3)
    >>> with numpy.printoptions(formatter={'all', fmt}):
        print(great_circle_distance([[54.0, 0.0], [55, 0.0]], [55, 1.0]))
    [1.286e+05 6.378e+04]
    """
    # If input is one list, make it list of lists
    if np.array(latlon1).ndim == 1:
        latlon1 = [latlon1]
    elif np.array(latlon2).ndim == 1:
        latlon2 = [latlon2]

    # Convert to radians
    latlon1 = np.array(latlon1) * np.pi/180.
    latlon2 = np.array(latlon2) * np.pi/180.

    # Earth radius given
    Rp = 6.371e6

    distance = np.empty((len(latlon1), len(latlon2)), float)

    for i in range(len(latlon1)):
        for j in range(len(latlon2)):
            # We set a criteria for considering 2 points on the same latitude or longitude
            # For example, if they are within 1e-5 radians of each other
            if abs(latlon1[i,1] - latlon2[j,1]) < 1e-5:
                distance[i,j] = abs(latlon1[i,0] - latlon2[j,0]) * Rp

            elif abs(latlon1[i,0] - latlon2[j,0]) < 1e-5:
                distance[i,j] = abs(latlon1[i,1] - latlon2[j,1]) * Rp * np.cos((latlon1[i,0]+latlon2[j,0])/2)

            else:
                # Harvesine formula
                #distance[i,j] = 2 * Rp * np.arcsin(np.sqrt((np.sin(abs(latlon1[i,0]-latlon2[j,0])/2)**2) + 
                 #   np.cos(latlon1[i,0]) * np.cos(latlon2[j,0]) * np.sin(abs(latlon1[i,1]-latlon2[j,1])/2)**2)) 
                
                # Spherical Vincenty formula
                numerator = np.sqrt((np.cos(latlon2[j,0])*np.sin(abs(latlon1[i,1]-latlon2[j,1])))**2 + 
                (np.cos(latlon1[i,0])*np.sin(latlon2[j,0]) - np.sin(latlon1[i,0])*np.cos(latlon2[j,0])*np.cos(abs(latlon1[i,1]-latlon2[j,1])))**2)
                
                denominator = np.sin(latlon1[i,0])*np.sin(latlon2[j,0]) + np.cos(latlon1[i,0])*np.cos(latlon2[j,0])*np.cos(abs(latlon1[i,1]-latlon2[j,1]))

                distance[i,j] = Rp * np.arctan2(numerator, denominator)

    return distance


class PostcodeLocator(object):
    """Class to interact with a postcode database file."""

    def __init__(self, postcode_file=postcodes_path,
                 census_file=population_path,
                 norm=great_circle_distance):
        """
        Parameters
        ----------

        postcode_file : str, optional
            Filename of a .csv file containing geographic
            location data for postcodes.

        census_file :  str, optional
            Filename of a .csv file containing census data by postcode sector.

        norm : function
            Python function defining the distance between points in latitude-longitude space.

        """
        #print(os.getcwd())
        #print(project_dir)

        self.norm = norm

        #Dataframe of all postcodes
        self.postcodes = pd.read_csv(postcode_file)

        #Added a number of units in sector dictionary for later use. The process takes about 0.6 seconds but speeds up postcode population lookup
        #from 0.6-1 seconds each to <0.01s.

        self.nUiSdict = self.postcodes["Postcode"].str[:-2].value_counts().to_dict()
        
        #PC = pd.read_csv(postcode_file)
        #self.postcodes = PC.sort_values(by=['Latitude'])

        #We read the csv file
        self.census = pd.read_csv(census_file)[['geography', 'Variable: All usual residents; measures: Value']]

        #we create a dict for sector:population
        self.census = dict(zip(self.census.geography, self.census.iloc[:, 1]))

        #We strip the keys of their whitespaces
        self.census = {k.replace(' ', ''): v for k, v in self.census.items()}

    def get_postcodes_by_radius(self, X, radii, sector=False):
        """
        Return (unit or sector) postcodes within specific distances of 
        input location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair of centre location 
        radii : arraylike
            array of radial distances from X
        sector : bool, optional
            if true return postcode sectors, otherwise postcode units

        Returns
        -------
        list of lists
            Contains the lists of postcodes closer than the elements of radii to the location X.


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.13e3])
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.4e3, 0.2e3], True)                                                                   
        """
        # If input is one list, make it list of lists
        if np.array(X).ndim == 1:
            X = [X]

        D = []
        for i in range(len(X)):
            # Grab Latidude and Longitude from full_postcodes.csv
            PC = self.postcodes[["Latitude", "Longitude"]].sort_values(by=['Latitude'])
            
            # M
            max_r = max(radii)

            # Radius of Earth
            Rp = 6.371e6

            # Calculate approximate latitude difference between impact and 
            # maximum damage radius to narrow postcode search
            lat_diff = max_r*180/(Rp*np.pi)*1.02

            # Postcode search narrowed to below maximum latitude where impact 
            # may be felt
            PC_range = PC.loc[PC.Latitude < X[i][0]+lat_diff]

            # Postcode search narrowed to above minimum latitude where impact 
            # may be felt
            PC_range = PC_range.loc[(PC_range.Latitude > X[i][0]-lat_diff)]

            # Calculate distances from X coordinates to postcodes within 
            # narrowed range of latitudes
            distances = pd.DataFrame(great_circle_distance(X[i], PC_range).T)
            for r in radii:
                
                post_sort = self.postcodes.sort_values(by=['Latitude'])

                # Same action of narrowing range of latitudes for search, 
                # so 'distances' and postcode DataFrame match dimensions
                p_sort_range = post_sort.Postcode.loc[post_sort.Latitude < X[i][0]+lat_diff]
                p_sort_range = p_sort_range.loc[post_sort.Latitude > X[i][0]-lat_diff]

                # Get postcodes at a distance smaller than radius
                d = p_sort_range.reset_index().loc[distances[0] < r].iloc[:,1]

                # When no postcode is found
                if sector == False:
                    D.append(d.tolist())
                elif sector == True:
                    # Convert units of a same sector into single postcode sector
                    d = d.str.slice(stop=-2).drop_duplicates().tolist()
                    D.append(d)

        return D
        



    def get_population_of_postcode(self, postcodes, sector=False):
        """
        Return populations of a list of postcode units or sectors.

        Parameters
        ----------
        postcodes : list of lists
            list of postcode units or postcode sectors 
        sector : bool, optional
            if true return populations for postcode sectors, otherwise postcode units

        Returns
        -------
        list of lists
            Contains the populations of input postcode units or sectors


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_population_of_postcode([['SW7 2AZ','SW7 2BT','SW7 2BU','SW7 2DD']])
        >>> locator.get_population_of_postcode([['SW7  2']], True)
        """
        # If input is one list, make it list of lists
        if (np.array(postcodes).ndim == 1 and len(postcodes)==1):
            postcodes = [postcodes]

        populations = []
        if sector == True:
            for zone in postcodes:
                pop = []
                for s in zone:
                    # Compare postcode sectors in population_by_postcode_sector.csv with
                    # input postcodes by reformatting postcodes in both removing
                    # all spaces in strings
                    try:
                        find_pop = self.census[s.replace(' ', '')]
                        pop.append(find_pop)
                    except:
                        # When postcode does not match to a population, assume 0
                        pop.append(0)
                populations.append(pop)
        else:
            for zone in postcodes:
                pop = []
                for s in zone:
                    # No data for unit population, so we assume uniform population 
                    # density within each sector, and calculate population in units 
                    # by dividing sector population by number of units. Line below
                    # calculates number of units in a postcode
                    sector = s[:-2]
                    nUis = 0
                    try:
                        nUiS = self.nUiSdict[sector]
                    except:
                        # When no postcode matches, divide sector population by 1
                        nUis = 1
                    # Compare postcode sectors in population_by_postcode_sector.csv
                    # with input postcodes by reformatting postcodes in both removing 
                    # all spaces in strings and changing units into sectors (i.e. remove
                    # last 2 digits of input postcodes)
                    try:
                        find_pop = self.census[s.replace(' ', '')[:-2]]
                        pop.append(find_pop/nUiS)
                    # When postcode does not match to a population, assume 0
                    except:
                        pop.append(0)
                populations.append(pop)
        return populations