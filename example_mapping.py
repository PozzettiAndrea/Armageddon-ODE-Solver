import armageddon

earth = armageddon.Planet(atmos_func='exponential')

# Solve the atmospheric entry problem (for something similar to Chelyabinsk). 
result = earth.solve_atmospheric_entry(radius=10, angle=20,
                                       strength=1e6, density=3000,
                                       velocity=19e3)

# Calculate the kinetic energy lost per unit altitude and add it
# as a column to the result dataframe
result = earth.calculate_energy(result)

# Determine the outcomes of the impact event
outcome = earth.analyse_outcome(result)

blast_lat, blast_lon, damage_rad = armageddon.damage_zones(outcome, 
                                                           lat=51.2, lon=0.7, bearing=-35.,
                                                           pressures=[1e3, 3.5e3, 27e3, 43e3])

#Plot the outcome on the map, with a certain bearing and distance
damage_map2 = armageddon.plot_outcome(blast_lat, blast_lon, damage_rad, 112, outcome["burst_distance"])
damage_map2.save("damage_map2.html")

#Plot a few imaginary (and imaginarily randomised) meteor impacts on the map
damage_map3 = armageddon.multiple_plots([[52.79, -2.95, 10e3],[52.78, -2.94, 12e3],[52.785, -2.943, 11e3]])
damage_map3.save("damage_map3.html")
