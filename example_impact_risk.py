#this example is exclusively for the impact_risk funtion within damage.py

import armageddon


#Values for something similar to the Brextinction event

fiducial_means = {'radius': 20, 'angle': 45, 'strength': 10e6,
                'density': 8000, 'velocity': 17e3,
                'lat': 53.79, 'lon': -3.55, 'bearing': 122.}

fiducial_stdevs = {'radius': 1, 'angle': 1, 'strength': 5e5,
                    'density': 50, 'velocity': 1e3,
                    'lat': 0.05, 'lon': 0.05, 'bearing': 1}
earth = armageddon.Planet(atmos_func='exponential')

print(armageddon.impact_risk(earth, means=fiducial_means, stdevs=fiducial_stdevs, pressure = 3000, nsamples=5, sector=True))
