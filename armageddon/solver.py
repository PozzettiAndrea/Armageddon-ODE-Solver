import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
# import timeit

from numba import jit


# t_s = timeit.default_timer()

class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(self, atmos_func='exponential',
                 atmos_filename="data/AltitudeDensityTable.csv",
                 Cd=1., Ch=0.1, Q=1e7, Cl=1e-3, alpha=0.3, Rp=6371e3,
                 g=9.81, H=8000., rho0=1.2):
        """
        Set up the initial parameters and constants for the target planet

        Parameters
        ----------
        atmos_func : string, optional
            Function which computes atmospheric density, rho, at altitude, z.
            Default is the exponential function rho = rho0 exp(-z/H).
            Options are 'exponential', 'tabular' and 'constant'

        atmos_filename : string, optional
            Name of the filename to use with the tabular atmos_func option

        Cd : float, optional
            The drag coefficient

        Ch : float, optional
            The heat transfer coefficient

        Q : float, optional
            The heat of ablation (J/kg)

        Cl : float, optional
            Lift coefficient

        alpha : float, optional
            Dispersion coefficient

        Rp : float, optional
            Planet radius (m)

        rho0 : float, optional
            Air density at zero altitude (kg/m^3)

        g : float, optional
            Surface gravity (m/s^2)

        H : float, optional
            Atmospheric scale height (m)

        """

        # Input constants
        self.Cd = Cd
        self.Ch = Ch
        self.Q = Q
        self.Cl = Cl
        self.alpha = alpha
        self.Rp = Rp
        self.g = g
        self.H = H
        self.rho0 = rho0

        # set function to define atmospheric density
        if atmos_func == 'exponential':
            self.rhoa = lambda z: rho0 * np.exp(-z / H)
        elif atmos_func == 'tabular':

            lines_to_ignore = 6
            ATM_list = open(atmos_filename).read().splitlines()[lines_to_ignore:]
            Altitude_array = np.zeros(len(ATM_list))
            Pressure_array = np.zeros(len(ATM_list))

            # Try Except Block for floating
            #####################################################

            try:
                for index in range(len(ATM_list)):
                    try:
                        Altitude_array[index] = float(ATM_list[index].split()[0])
                        Pressure_array[index] = float(ATM_list[index].split()[1])
                    except:
                        print('Some values cannot be converted into float')
                        raise ValueError
            except:
                print('There seems to be some bad data entries')
                raise ValueError

            ####################################################
            # Using CubicSpline in scipy to do CubicSpline for ATM tabular
            # I note that there is a problem because the data table only goes to 8.6e4m
            # And Cubic Spline extrapolation will result in pressure going into negative
            # When the altitude is above 8.6e4
            # Thus I have set the density to be 0 above 8.6e4
            # Since density there is already about 1 millionth of a Pascal
            ####################################################

            interpol_3 = CubicSpline(Altitude_array, Pressure_array)
            self.rhoa = lambda z: 0 if interpol_3(z) < 0 else interpol_3(z)

        elif atmos_func == 'constant':
            self.rhoa = lambda z: rho0
        else:
            raise NotImplementedError(
                "atmos_func must be 'exponential', 'tabular' or 'constant'")

    ##This solves the Atmospheric Entry

    def solve_atmospheric_entry(
            self, radius, velocity, density, strength, angle,
            init_altitude = 1e5, dt = 0.05, radians=False):
        """
        Solve the system of differential equations for a given impact scenario

        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters

        velocity : float
            The entery speed of the asteroid in meters/second

        density : float
            The density of the asteroid in kg/m^3

        strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2

        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians

        init_altitude : float, optional
            Initial altitude in m

        dt : float, optional
            The output timestep, in s

        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False
            Angles returned in the dataframe will have the same units as the
            input

        Returns
        -------
        Result : DataFrame
            A pandas dataframe containing the solution to the system.
            Includes the following columns:
            'velocity', 'mass', 'angle', 'altitude',
            'distance', 'radius', 'time'
                    
        Examples
        -------
        >>> Earth = Planet()
        >>> dataframe = (Earth.solve_atmospheric_entry(radius = 20., velocity = 17e3, density = 8000.,
                                    strength = 1e7, angle = 45, init_altitude=1e6))
        >>> print(dataframe.iloc[-1, :])
        velocity    9.786400e+02
        mass        5.711703e+07
        angle       3.816044e+01
        altitude    5.012545e+03
        distance    1.045633e+06
        radius      1.781137e+02
        time        8.865000e+01
        Name: 1773, dtype: float64
        
        >>> Earth = Planet()
        >>> dataframe = (Earth.solve_atmospheric_entry(radius = 20., velocity = 17e3, density = 8000.,
                                    strength = 1e7, angle = 6, init_altitude=1e6))
        >>> print(dataframe.iloc[-1, :])
        velocity    1.703099e+04
        mass        3.351032e+07
        angle      -4.591476e-03
        altitude    9.462484e+05
        distance    8.883275e+05
        radius      1.000000e+01
        time        6.020000e+01
        Name: 1204, dtype: float64
        
        """

        # Initilization
        # step 1: setting up new variables
        # prevent change to old variables
        v_old = velocity
        m_old = density * 4 / 3 * np.pi * (radius ** 3)
        z_old = init_altitude
        x_old = 0
        r_old = radius
        t_old = 0

        count = 0 # This counts to ensure that that output is in dt
        atmos_max = 70000 # The maximum height of the ATM
        min_mass = 1e-3 # The minimum mass that we are interested in

        #Step1, extra, Converting to radians if the input is degrees
        #Will convert back from radians

        if radians == False:
            theta_old = np.radians(angle)
        else:
            theta_old = angle

        # Step2 , creating the list to store the values
        # arrays could be created too,
        # but list.append seems ok compared to numpy.append

        vlist = [v_old, ]
        mlist = [m_old, ]
        xlist = [x_old, ]
        zlist = [z_old, ]
        rlist = [r_old, ]
        thetalist = [theta_old, ]
        tlist = [t_old, ]

        # Step3, in case bugs happen
        # Give some info
        if False:
            print(vlist)
            print(mlist)
            print(xlist)
            print(zlist)
            print(rlist)
            print(thetalist) 
            print(tlist)
            print(init_altitude)
            print(dt)
            print(radians)

            print(self.Cd)
            print(self.Ch)
            print(self.Q)
            print(self.Cl)
            print(self.alpha)
            print(self.Rp)
            print(self.g)
            print(self.H)
            print(self.rho0)

        # Step4: solving the ODEs

        # The asteroid path is going to be divided into 3 to 4 stages (deoend on airburst)
        # And different set of ODEs to be solved for each stage

        # Stage I, the asteroid is so far away that the Earth curvature does not matter
        # i.e. z >>>> Rp
        # i.e. At this point, the gravitational effect is also minimal 
        # Using very large timesteps here
        # Stage I is not implemented since such distance will usually result in another
        # Planetary or stellar body eating that up
        # Such as Jupyter or the Sun 

        # Stage II, the asteroid moves closer and becomes closer to the planet
        # But the asteroid is still outside the ATM of Earth, so the pressure = 0
        # i.e. z > 80000  
        # Using kinda large timesteps here

        # Stage III, the asteroid moves through the ATM, here we need full computational power
        # so here the timestep is going to be very very very small 
        # This way we can have accurate simulation 
        # Using very very small timesteps here

        # Stage IV, post breakup 
        # If the asteroid breaksup, then things reach terminal velocity
        # Using large timesteps here

        ## Think about Zeno Phenomenon



        # Begin stage II
        # I am assuming that it is for the Atmosphere of Earth
        # Air density of Earth becomes approx 0 above 70000
        # Since no air, then dmdt = 0 and drdt = 0
        # We will only need to solve 4 ODEs for 
        # v, theta, x ,t

        count = 1
        # Count continues through the stages!!!
        ###########################################
        ###########################################
        # Stage II timestep can remain relatively large
        # Since we are expecting little to no change
        # Need to remember to populate rlist and mlist

        n = 1
        timestep = dt * (1/n) 

        while theta_old > 0:
            if z_old < atmos_max:
                break
            else:
                sintheta = np.sin(theta_old)
                costheta = np.cos(theta_old)

                dzdt = (-1) * v_old * sintheta
                dvdt = self.g * sintheta
                dxdt = v_old * costheta / (1 + z_old/self.Rp)
                dthetadt = self.g * costheta/v_old +\
                           (-1) * v_old * costheta/(self.Rp + z_old)


                z_new = z_old + timestep * dzdt
                v_new = v_old + timestep * dvdt
                x_new = x_old + timestep * dxdt
                theta_new = theta_old + timestep * dthetadt
                t_new = t_old + timestep

                # print('This is new z ',z_new)
                # print('This is new v ',v_new)
                # print('This is new x ',x_new)
                # print('This is new theta ',theta_new)
                
                z_old = z_new
                v_old = v_new
                x_old = x_new
                theta_old = theta_new
                t_old = t_new
                
                if mlist[-1] < min_mass:
                    print('Timestep here is way too large')
                    print('Velocity seems unrealistic')
                    break
                else:
                    pass

                if count == n:
                    zlist.append(z_new)
                    vlist.append(v_new)
                    xlist.append(x_new)
                    thetalist.append(theta_new)
                    tlist.append(t_new)
                    count = 1
                else:
                    count += 1

        # Populating mlist and rlist
        mlist = [m_old] * len(vlist)
        rlist = [r_old] * len(vlist)
    
        # Will it airburst? 
        # Airburst is the comparison between velocity and Strength
        # I think about the best scenario for airburst to happen
        # And if airburst does not happen then
        # No airburst can happen
        # I know that faster velocity leads to airburst
        # g I will assume is a constant
        # For faster velocity, I will assume that there is no drag
        # Then dvdt = gsin(\theta), and sin is periodic, so max sin = 1
        # Then max dvdt = g thus if no drag, straight down, then 
        # acceleration = g, and then I know that the object falls from init_altitude
        # Then I have v_max = np.sqrt(velocity^2 + 2*g*init_altitude)
        # Thus plugging this v_max into the rho0 * v^2 compare to strength
        # And using rho0 = 1.2, if it is still not more than strength
        # Then defo no airburst

        v_max = np.sqrt(velocity ** 2 + 2 * self.g * init_altitude)
        if (v_max**2) * (self.rho0) > strength:
            maybe_airburst = True
        else:
            maybe_airburst = False       

        # print('Length of vlist at Stage II ', len(vlist))
        # print('Length of rlist at Stage II ', len(rlist))
        # print('Will airburst happen? ',maybe_airburst)

        #####################################################
        # Stage III timesteps are going to be very small
        # Since we are expecting rapid change in this stage
        # Airburst could or could not happen in Stage III
        # If airburst cannot happen, then drdt remains 0
        #####################################################

        n = 10
        timestep = dt * (1/n) 
        drdt = 0
        airburst_happening = False
        airburst_happened = False

        if maybe_airburst == False:
            ## No airburst situation

            while theta_old > 0:
                if z_old < 0:
                    break
                else:
                    sintheta = np.sin(theta_old)
                    costheta = np.cos(theta_old)
                    area = np.pi * (r_old ** 2)

                    dzdt = (-1) * v_old * sintheta
                    dvdt = self.g * sintheta + \
                           (-1)*self.Cd*self.rhoa(z_old)*area*(v_old**2)/(2*m_old)
                    dmdt = (-1)*self.Ch*self.rhoa(z_old)*area*(v_old**3)/(2*self.Q)
                    dxdt = v_old * costheta / (1 + z_old/self.Rp)
                    dthetadt = self.g * costheta/v_old +\
                            (-1) * v_old * costheta/(self.Rp + z_old)

                    z_new = z_old + timestep * dzdt
                    v_new = v_old + timestep * dvdt
                    x_new = x_old + timestep * dxdt
                    theta_new = theta_old + timestep * dthetadt
                    t_new = t_old + timestep
                    m_new = m_old + timestep * dmdt

                    # print('This is new z ',z_new)
                    # print('This is new v ',v_new)
                    # print('This is new x ',x_new)
                    # print('This is new theta ',theta_new)
                    
                    z_old = z_new
                    v_old = v_new
                    x_old = x_new
                    theta_old = theta_new
                    t_old = t_new
                    m_old = m_new

                    if mlist[-1] < min_mass:
                        print('Timestep here is way too large')
                        print('Velocity seems unrealistic')
                        break
                    else:
                        pass

                    if count == n:
                        zlist.append(z_new)
                        vlist.append(v_new)
                        xlist.append(x_new)
                        thetalist.append(theta_new)
                        tlist.append(t_new)
                        mlist.append(m_old)
                        count = 1
                    else:
                        count += 1

            # Populating mlist and rlist
            rlist = [r_old] * len(vlist)

        # We understand airburst above 5000
        # But we do not understand airburst below 5000
        # So this simulator does not try to do that
        airburst_altitude = 5000
        if maybe_airburst == True:
            ## Maybe Airburst situation

            while theta_old > 0:
                if z_old < 5000:
                    break
                else:
                    sintheta = np.sin(theta_old)
                    costheta = np.cos(theta_old)
                    area = np.pi * (r_old ** 2)

                    dzdt = (-1) * v_old * sintheta
                    dvdt = self.g * sintheta + \
                           (-1)*self.Cd*self.rhoa(z_old)*area*(v_old**2)/(2*m_old)
                    dmdt = (-1)*self.Ch*self.rhoa(z_old)*area*(v_old**3)/(2*self.Q)
                    dxdt = v_old * costheta / (1 + z_old/self.Rp)
                    dthetadt = self.g * costheta/v_old +\
                            (-1) * v_old * costheta/(self.Rp + z_old)
                    
                    # print('Did the airburst happen?',airburst_happened)
                    if airburst_happened == False:

                        if self.rhoa(z_old) * (v_old**2) < strength:
                            # print('No breakup')
                            r_new = r_old
                            if airburst_happening == True:
                                airburst_happened = True
                                

                        if self.rhoa(z_old) * (v_old**2) > strength:
                            # print('Yes breakup')
                            airburst_happening = True
                            drdt = np.sqrt(3.5 * self.alpha * self.rhoa(z_old) / density) * v_old
                            r_new = r_old + timestep * drdt
                            r_old = r_new
                    else:

                        r_new = r_old

                    # print(r_new)

                    z_new = z_old + timestep * dzdt
                    v_new = v_old + timestep * dvdt
                    x_new = x_old + timestep * dxdt
                    theta_new = theta_old + timestep * dthetadt
                    t_new = t_old + timestep
                    m_new = m_old + timestep * dmdt

                    # print('This is new z ',z_new)
                    # print('This is new v ',v_new)
                    # print('This is new x ',x_new)
                    # print('This is new theta ',theta_new)
                    
                    z_old = z_new
                    v_old = v_new
                    x_old = x_new
                    theta_old = theta_new
                    t_old = t_new
                    m_old = m_new

                    if z_new < airburst_altitude or m_new < min_mass:
                        break

                    if count == n:
                        zlist.append(z_new)
                        vlist.append(v_new)
                        xlist.append(x_new)
                        thetalist.append(theta_new)
                        tlist.append(t_new)
                        mlist.append(m_new)
                        rlist.append(r_new)
                        count = 1
                    else:
                        count += 1

            # Populating mlist and rlist
            # rlist = rlist + rlist[-1] * (len(vlist) - len(rlist))


        # Step 6 Preparing the output
        # Ensuring that degrees in, degrees out
        # Ensures that radians in, radians out
        # Returns pandas. dataframe

        thetaarray = np.array(thetalist)

        if radians == False:
            thetaarray = np.degrees(thetaarray)

        data = {'velocity': vlist,
                'mass': mlist,
                'angle': thetaarray,
                'altitude': zlist,
                'distance': xlist,
                'radius': rlist,
                'time': tlist}

        dataframe = pd.DataFrame(data=data, index=range(len(vlist)))

        return dataframe

    def calculate_energy(self, result):
        """
        Function to calculate the kinetic energy lost per unit altitude in
        kilotons TNT per km, for a given solution.
        Parameters
        ----------
        result : DataFrame
            A pandas dataframe with columns for the velocity, mass, angle,
            altitude, horizontal distance and radius as a function of time
        Returns : DataFrame
            Returns the dataframe with additional column ``dedz`` which is the
            kinetic energy lost per unit altitude
        """
        result_np = result.to_numpy()

        # calculate kinetic energy (0.5*m*v**2)
        kinetic_energy = 0.5 * result_np[:, 1] * result_np[:, 0] ** 2
        # find change in kinetic energy between every timestep
        diff_kinetic_energy = kinetic_energy[:-1] - kinetic_energy[1:]
        # calculate hight difference at every timestep
        diff_altitude = result_np[:-1, 3] - result_np[1:,3]

        # find kinetic energy per length
        kinetic_energy_per_km = diff_kinetic_energy / diff_altitude
        # put zero as first point to ensure it is same length as result
        kinetic_energy_per_km = np.insert(kinetic_energy_per_km, 0, 0)

        # conversion factor 1 j to 1 TNT
        conv_j_ktnt = 2.39006e-13
        # conversion to correct units
        kinetic_energy_per_km = conv_j_ktnt * 1000 * kinetic_energy_per_km


        # Replace these lines with your code to add the dedz column to
        # the result DataFrame
        result = result.copy()
        result = pd.concat([result, pd.DataFrame({'dedz': kinetic_energy_per_km})], axis=1)

        return result

    def analyse_outcome(self, result):
        """
        Inspect a pre-found solution to calculate the impact and airburst stats
        Parameters
        ----------
        result : DataFrame
        pandas dataframe with velocity, mass, angle, altitude, horizontal
        distance, radius and dedz as a function of time
        dt : float
             The output timestep, in s
        Returns
        -------
          outcome : Dict
             dictionary with details of the impact event, which should contain
             the key ``outcome`` (which should contain one of the following strings:
             ``Airburst``, ``Cratering`` or ``Airburst and cratering``), as well as
             the following keys:
             ``burst_peak_dedz``, ``burst_altitude``, ``burst_distance``, ``burst_energy``
        """


        # conversion factor 1 j to 1 TNT
        conv_j_ktnt = 2.39006e-13


        # if no timesteps are run (i.e. in the tests) return a dictornary of zeros
        # this needs to be done to avoid index errors that occour when when lists of
        # size zero are passed in the function
        if len(result['dedz']) <= 1:
            outcome = 'no timestepping happened'
            burst_peak_dedz = 0
            burst_altitude = 0
            burst_distance = 0
            burst_energy = 0

            outcome = {'outcome': outcome,
                   'burst_peak_dedz': burst_peak_dedz,
                   'burst_altitude': burst_altitude,
                   'burst_distance': burst_distance,
                   'burst_energy': 1000 * burst_energy * conv_j_ktnt}

            return  outcome

        # find index of maximum kinetic energy loss per km
        # this will give the index for the burst event
        # or if cratering event it will give the last timestep
        index_burst_peak = result['dedz'].idxmax()


        # check altitude of airburst event to determine the outcome type
        if result.loc[index_burst_peak, 'altitude'] > 5000:
            outcome = 'Airburst'
        elif result.loc[index_burst_peak, 'altitude'] < 300:
            outcome = 'Cratering'
        else:
            outcome = 'Airburst and cratering'


        # if airbust event
        if outcome != 'Cratering':
            # find maximum kinectic energy loss per km
            burst_peak_dedz = result['dedz'].iloc[index_burst_peak]
            # use index of maximum engery to find altitude
            burst_altitude = result['altitude'].iloc[index_burst_peak]
        # if cratering event, there is no peak dedz and altitude
        else:
            # no burst altitude since it is cratering event
            burst_altitude = 0.0
            # burst dEdz is equal to dEdz at impact i.e. last timestep
            burst_peak_dedz = result['dedz'].iloc[-1]

        # calculate kinetic energy for every timestep
        kinetic_energy = 0.5 * result.loc[:, 'mass'] * result.loc[:, 'velocity'] ** 2
        # Find index of break up point by finding the first index
        # where the radius is larger than at the beginning
        index_break_up = result[result['radius'] > result['radius'].iloc[0]].index.values
        # index_break_up is of length zero if the radius doesnt chnage so no break up happes
        # to ensure that there is no error, append zero if length is zero
        if len(index_break_up) == 0:
            index_break_up = np.append(index_break_up, 0)
        # find the kinetic energy lost from break up point until the air burst
        kinetic_energy_lost = kinetic_energy[index_break_up[0]] - kinetic_energy[index_burst_peak]
        # kinetic energy of the asteroid at last time step
        kinetic_energy_remain = 0.5 * (result['mass'].iloc[-1]) * (result['velocity'].iloc[-1]) ** 2
        # find maximum of either two
        burst_energy = np.maximum(kinetic_energy_lost, kinetic_energy_remain)

        # find the distance traveled by the asteriod until burst event or crash
        burst_distance = result.loc[index_burst_peak, 'distance']

        outcome = {'outcome': outcome,
                   'burst_peak_dedz': burst_peak_dedz,
                   'burst_altitude': burst_altitude,
                   'burst_distance': burst_distance,
                   'burst_energy': burst_energy * conv_j_ktnt}

        return outcome

    def find_parameters(self, solver_tol = 1):
        """
        Find asteroid size and strength of asteroid from existing Energy curve
        by finding the minimal differnece between the measured and calculated data

        Parameters
        ----------
        solver_to : float
        tolerance of the optimization

        Returns
        -------
          cond : list
          list of optimized conditions for asteroid values

          min_par : list
          list of error values for the peak and Altitude for the
          optimized solution

        """
        # import Chelyabinsk Energy vs Altitude data and
        # find maxiumum energy and corresponding height
        chelyab_curve = np.loadtxt('./data/ChelyabinskEnergyAltitude.csv',
                                   delimiter=',', skiprows=1)
        chelyab_height = 1000 * chelyab_curve[:,0]
        chelyab_energy = chelyab_curve[:,1]
        chelyab_peak_id = chelyab_energy.argmax()
        chelyab_energy_max = chelyab_energy[chelyab_peak_id]
        chelyab_height_max = chelyab_height[chelyab_peak_id]

        # set up class
        Earth = Planet(atmos_func='exponential')

        # initialize error metric to a high value for first itteration
        tol = 12

        # initialize conditions for asteroid
        # inital radius is put to 18 m because it is an airbust event so not large asteroid
        # inital strength is put to 3.5e6 due to the measured for recovered meteorites
        cond = [18, 19.2e3, 3300, 3.5e6, 18.3, 55000, 0.05, False]

        # do loop until tolerance reached
        # first loop is to optimize the dedz value
        while solver_tol < abs(tol):

            # use parameters (cond) to find error metrics
            min_par = Earth.minimise(cond, chelyab_height_max, chelyab_energy_max)

            # if error is negative, meteorites needs to be bigger
            if tol < 0:
                tol = abs(min_par[0])
                cond[0] = cond[0]*Earth.learing_rate(tol) #use of adaptive learing rate
                tol = (min_par[0])
            # if error is positive, meteorites needs to be smaller
            else:
                tol = abs(min_par[0])
                cond[0] = cond[0]/Earth.learing_rate(tol)
                tol = (min_par[0])

        # initialize error metric for next loop
        tol = abs(min_par[1])

        # do loop until tolerance reached
        # this loop is to optimize the height of the maxiumum dedz
        while abs(min_par[1]) > 800:

            # this will also effect the max dedz so if is gets to large break
            if abs(min_par[0]) > 10:
                break

            # use parameters (cond) to find error metrics
            min_par = Earth.minimise(cond, chelyab_height_max, chelyab_energy_max)

            # if error is negative, meteorites needs to be stronger
            if tol < 0:
                tol = abs(min_par[1])
                cond[3] = cond[3]/Earth.learing_rate(tol * 0.01)
                tol = (min_par[1])
            # if error is postive, meteorites needs to be less stronger
            else:
                tol = abs(min_par[1])
                cond[3] = cond[3]*Earth.learing_rate(tol * 0.01)
                tol = (min_par[1])

        # find percentage error
        min_par[0] = (100 - 100 * abs(chelyab_energy_max - min_par[0]) / chelyab_energy_max)
        min_par[1] = (100 - 100 * abs(chelyab_height_max - min_par[1]) / chelyab_height_max)

        return cond, min_par

    def learing_rate(self, error_est):
        """
        adaptive learing rate that depends on error.
        Large error means larger learning rate

        Parameters
        ----------
        error_est : float
        error estimate from minimise function

        Returns
        -------
          learning_rate : float
          learning rate       
        """
        return 0.008*np.sqrt(error_est)+1

    def minimise(self, condition, chelyab_height_max, chelyab_energy_max):
        """
        Find the error between the peak dedz value from solver and data
        and error between the burst height from solver and data. These two
        values need to minimized

        Parameters
        ----------
        condition : list
        list of values that define the asteroid

        chelyab_height_max : float
        maximum dedz value from data

        chelyab_energy_max : float
        the height of the burst point from data

        Returns
        -------
          minimise_parameter : list
          list of errors that need to minimized
        """

        # use condition to run ODE solver and analyse outcome
        Earth = Planet(atmos_func='exponential')
        result = Earth.solve_atmospheric_entry(*condition)
        result = Earth.calculate_energy(result)
        outcome = Earth.analyse_outcome(result)

        # find difference between values
        minimise_parameter = np.array([outcome['burst_peak_dedz'] - chelyab_energy_max,
                                        outcome['burst_altitude'] - chelyab_height_max])


        return minimise_parameter




#Earth = Planet(atmos_func='exponential')
#a, error= Earth.find_parameters(1)
#print(a)
#print('error')
#print(error)

#Earth1 = Planet(atmos_func='exponential')
#result = Earth1.solve_atmospheric_entry(*a)
#result = Earth.calculate_energy(result)
#outcome = Earth.analyse_outcome(result)
#print(outcome)

#chelyab_curve = np.loadtxt('./data/ChelyabinskEnergyAltitude.csv', delimiter=',',  skiprows=1)
#chelyab_height = 1000*chelyab_curve[:,0]
#chelyab_energy = chelyab_curve[:,1]


#import matplotlib.pyplot as plt
#result.plot(y = 'altitude', x ='dedz',label='Fitted parameter solution', yticks = [0, 10000, 20000, 30000, 40000, 50000])
#plt.plot(chelyab_energy, chelyab_height, label='Chelyabinsk data')
#plt.legend()
#plt.xlabel('energy loss per unit height [kt/km]')
#plt.ylabel('altitude [m]')
#plt.title('Energy loss per unit height for optimization solution')
#plt.show()

# Solve the atmospheric entry problem (for something similar to Chelyabinsk). 

#Earth = Planet(atmos_func='exponential')
#result = Earth.solve_atmospheric_entry(radius=10, angle=20,strength=1e6, density=3000,velocity=19e3, init_altitude = 100000, dt = 0.05)

# # Calculate the kinetic energy lost per unit altitude and add it
# # as a column to the result dataframe
#result = Earth.calculate_energy(result)


# # Determine the outcomes of the impact event
#outcome = Earth.analyse_outcome(result)
#print(outcome)
#import matplotlib.pyplot as plt
#result.plot(y = 'altitude', x ='dedz', yticks = [0, 10000, 20000, 30000, 40000, 50000])
#plt.show()


# Initialise the Planet class
#earth = Planet(atmos_func='exponential')

# Solve the atmospheric entry problem (for something similar to Chelyabinsk). 
#result = earth.solve_atmospheric_entry(radius=10, angle=20,
                                    #  strength=1e6, density=3000,
                                    # velocity=19e3)

# Calculate the kinetic energy lost per unit altitude and add it
# as a column to the result dataframe
#result = earth.calculate_energy(result)

# Determine the outcomes of the impact event
#outcome = earth.analyse_outcome(result)
#print(outcome)


#import matplotlib.pyplot as plt
#result.plot(y = 'altitude', x ='dedz', yticks = [0, 10000, 20000, 30000, 40000, 50000, 60000, 70000])
#plt.show()

# Earth = Planet()

# dataframe = (Earth.solve_atmospheric_entry(radius = 3.0, velocity = 19.2e3, density = 9300,
#                                     strength = 1e7, angle = 25, init_altitude=4.5e5, dt = 0.05))

# t_e = timeit.default_timer()
# print('Total time : {:.2f}'.format(t_e-t_s))

# dataframe.to_csv('altitude_bug.csv')

# zlist = dataframe.loc[:,'altitude'].values
# xlist = dataframe.loc[:,'distance'].values
# vlist = dataframe.loc[:,'velocity'].values
# mlist = dataframe.loc[:,'mass'].values
# thetalist = dataframe.loc[:,'angle'].values
# rlist = dataframe.loc[:,'radius'].values

# print(len(zlist))

# if 1:
#     plt.plot(vlist,zlist)
#     plt.title('Velocity change with Height')

#     plt.grid()
#     plt.savefig('img/fig1.png')
#     plt.show()
#     plt.close()

#     plt.plot(mlist,zlist)
#     plt.title('Mass change with Height')

#     plt.grid()
#     plt.savefig('img/fig2.png')
#     plt.show()
#     plt.close()

#     plt.plot(thetalist,zlist)
#     plt.title('Angle change with Height')
#     plt.grid()
#     plt.savefig('img/fig3.png')
#     plt.show()
#     plt.close()


#if __name__ == "__main__":
#    import doctest
 #   doctest.testmod()
