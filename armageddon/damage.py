import numpy as np
import pandas as pd
import armageddon
from scipy.optimize import fsolve

def f(k, p):
    r""" Airblast empirical function.
    
    k is substituted in the place of the parethesis:
     .. math::

        k := \begin{equation}
            \frac{r^2+z^2}{E^{2/3}}
            \end{equation}
    So our equation is the following:
     .. math::

        p = \begin{equation}
            3.14*10^{11}*k^{-1.3}+1.8*10^7*k^{-0.565}
            \end{equation}
    
    Parameters
    ----------
    k: float
        aggregated parameter which we want to calculate later on
    p: float
        pressure
    
    Returns
    -------
    eq: float

    Returns
    -------
    val: float
    
    Examples
    --------
    >>> import armageddon
    >>> armageddon.f(49652.59, 43e3)
    0.0008996416727313772
    """
    eq = 3.14*np.power(10,11) * np.power(k,-1.3) + 1.8*np.power(10,7)*np.power(k, -0.565) - p
    return eq


def central_diff(f, k, p, dx=1e-6):
    r""" Calculating the central difference for the given function at the given points.
    
    Parameters
    ----------
    f: function
        The function which we want to differentiate.
    k: float
        Aggregate number for the radius, blast Energy and burst altitude.
    p: float
        Pressure
    dx: float
        :math:`\Delta x`

    Returns
    -------
    der: float
        Value of the derivative of the given function at the given points
    
    Examples
    --------
    >>> import armageddon
    >>> function = armageddon.f(49652.59, 43e3)
    >>> armageddon.central_diff(function, 10000, 43e3)
    -8.720089681446552
    """
    fxph = f(k + dx, p)
    fxnh = f(k - dx, p)
    der = (fxph - fxnh) / (2 * dx)
    return der


def newton(f, k, p, atol=1.0e-6):
    """Solver for nonlinear equation using Newton-Raphson method.
    
    Parameters
    ----------
    f: function
        The function which we want to solve.
    k: float
        Initial guess for the aggregate number which contains the radius, blast Energy and burst altitude.
    p: float
        Pressure
    atol: float
        Tolerance

    Returns
    -------
    x: float
        Solution for k in the given equation.
    
    Examples
    --------
    >>> import armageddon
    >>> armageddon.newton(armageddon.f, 10000, 43e3)
    49652.59168566593
    """
    x = [k]
    while True:
        x.append(x[-1] - f(x[-1], p)/central_diff(f, x[-1], p))
        if abs(x[-1]-x[-2]) < atol:
            return x[-1]


def calculate_r(k, E, z):
    """Calculate r from the substitution.

    Parameters
    ----------
    k: float
        solution for the previous equation
    E: float
        burst energy
    z: float
        burst altitude
    
    Returns
    -------
    r: float
        blast radius
        
    Examples
    --------
    >>> import armageddon
    >>> armageddon.calculate_r(49652.59, 6000., 9000.)
    0
    """
    # We'll check for negative values of r as well.
    # Physically a negative value of r means that there are no points on the ground that experience a certain pressure.
    r2 = k * np.power(E,2./3) - np.power(z, 2)
    if r2 < 0:
        r = 0
    else:
        r = np.sqrt(r2)
    return r


def radii(E, z, ps):
    """Determining the radius for the given zones.
    
    Parameters
    ----------
    E: float
        burst energy
    z: float
        burst altitude
    ps: arraylike
        pressures (Pa)

    Returns
    -------
    rs: arraylike
        radius for the damage zones for the given pressure levels
        
    Examples
    --------
    >>> import armageddon
    >>> armageddon.radii(6000., 9000., ps = [1000., 3500., 27e3, 43e3])
    [105627.75747908822, 33887.521090799775, 0, 0]
    """
    rs = []
    for i in ps:
        rs.append(calculate_r(newton(f, 10000, i),E, z))
    return rs


def damage_zones(outcome, lat, lon, bearing, pressures):
    """
    Calculate the latitude and longitude of the surface zero location and the
    list of airblast damage radii (m) for a given impact scenario.
    
    Parameters
    ----------
    outcome: Dict
        the outcome dictionary from an impact scenario
    lat: float
        latitude of the meteoroid entry point (degrees)
    lon: float
        longitude of the meteoroid entry point (degrees)
    bearing: float
        Bearing (azimuth) relative to north of meteoroid trajectory (degrees) 
    pressures: float, arraylike
        List of threshold pressures to define airblast damage levels

    Returns
    -------
    blat: float
        latitude of the surface zero point (degrees)
    blon: float
        longitude of the surface zero point (degrees)
    damrad: arraylike, float
        List of distances specifying the blast radii for the input damage levels

    Examples
    --------
    >>> import armageddon
    >>> outcome = {'burst_altitude': 8e3, 'burst_energy': 7e3,
                'burst_distance': 90e3, 'burst_peak_dedz': 1e3,
                'outcome': 'Airburst'}
    >>> armageddon.damage_zones(outcome, 52.79, -2.95, 135, pressures=[1e3, 3.5e3, 27e3, 43e3])
    (52.21396905216966, -2.015908861677074, [111312.95097589909, 36033.62609582931, 0, 0])
    """
    Rp = 6.371e6
    r = outcome["burst_distance"]
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    bearing = np.deg2rad(bearing)
    blat = np.arcsin(np.sin(lat)*np.cos(r/Rp)+np.cos(lat)*np.sin(r/Rp)*np.cos(bearing))
    blon = lon + np.arctan((np.sin(bearing)*np.sin(r/Rp)*np.cos(lat))/(np.cos(r/Rp)-np.sin(lat)*np.sin(blat)))
    damrad = radii(outcome["burst_energy"],outcome["burst_altitude"], ps = pressures)
    return float(np.degrees(blat)), float(np.degrees(blon)), damrad



fiducial_means = {'radius': 10, 'angle': 20, 'strength': 1e6,
                'density': 3000, 'velocity': 19e3,
                'lat': 51.5, 'lon': 1.5, 'bearing': -45.}

fiducial_stdevs = {'radius': 1, 'angle': 1, 'strength': 5e5,
                    'density': 500, 'velocity': 1e3,
                    'lat': 0.025, 'lon': 0.025, 'bearing': 0.5}


def impact_risk(planet, means=fiducial_means, stdevs=fiducial_stdevs,
                pressure = 1000, nsamples=100, sector=True):
    """
    Perform an uncertainty analysis to calculate the risk for each affected
    UK postcode or postcode sector

    Parameters
    ----------
    planet: armageddon.Planet instance
        The Planet instance from which to solve the atmospheric entry

    means: dict
        A dictionary of mean input values for the uncertainty analysis. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    stdevs: dict
        A dictionary of standard deviations for each input value. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    pressure: float
        The pressure at which to calculate the damage zone for each impact

    nsamples: int
        The number of iterations to perform in the uncertainty analysis

    sector: logical, optional
        If True (default) calculate the risk for postcode sectors, otherwise
        calculate the risk for postcodes

    Returns
    -------
    risk: DataFrame
        A pandas DataFrame with columns for postcode (or postcode sector) and
        the associated risk. These should be called ``postcode`` or ``sector``,
        and ``risk``.
        """

    # A list to hold all the postcodes or sectors for each blast zone
    BlastZones = []
    for i in range(nsamples):

        #Randomising each of our values, mean + randomshift*stdev
        radius = means["radius"] + np.random.randn()*stdevs["radius"]
        velocity = means["velocity"] + np.random.randn()*stdevs["velocity"]
        density = means["density"] + np.random.randn()*stdevs["density"]
        strength = means["strength"] + np.random.randn()*stdevs["strength"]
        angle = means["angle"] + np.random.randn()*stdevs["angle"]
        lat = means["lat"] + np.random.randn()*stdevs["lat"]
        lon = means["lon"] + np.random.randn()*stdevs["lon"]
        bearing = means["bearing"] + np.random.randn()*stdevs["bearing"]

        #Finding the surface zero location and blast radius from our randomised values
        result = planet.solve_atmospheric_entry(radius, velocity, density, strength, angle, init_altitude=100e3, dt=0.05, radians=False)
        outcome = planet.analyse_outcome(planet.calculate_energy(result))
        d_zones = damage_zones(outcome, lat, lon, bearing, pressures = [pressure])
        BlastZones.append([d_zones[0],d_zones[1],d_zones[2]])
    
    Locator = armageddon.locator.PostcodeLocator()
    AffectedAreas = []
    for bz in BlastZones:
        #Finding the postcodes or sectors in our blast zone
        BZ = Locator.get_postcodes_by_radius([[bz[0],bz[1]]],bz[2],sector=sector)

        #Append the postcodes or sectors within our blast zone to the original array, as a Pandas series
        AffectedAreas.append(pd.Series(BZ[0]))
    
    #Joining our series in a single Pandas series, then obtaining a dataframe with the individual value counts for each postcode or sector
    nameofcolumn = ""
    if sector == True:
        nameofcolumn = "Sector"
    else:
        nameofcolumn = "Postcodes"
    
    df = pd.concat(AffectedAreas).value_counts().rename_axis(nameofcolumn).reset_index(name='risk')

    #Dividing the individual counts for postcodes by the number of samples 
    df['risk'] *= (1/nsamples)

    #Finding the population for each postcode
    Pops =  pd.Series(Locator.get_population_of_postcode([df[nameofcolumn].to_numpy()],sector=sector)[0])

    #Multiplying the risk of each postcode being hit by its population
    df['risk'] = df['risk'] * Pops

    return df