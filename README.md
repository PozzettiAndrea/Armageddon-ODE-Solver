# Armageddon - Asteroid atmospheric entry and impact modelling
This is package simulates the behaviour of an asteroid that enters the atmosphere, breaks up and hits the ground. The dynamics of the asteroid is described by a coupled system of ordinary differential equations. These equations are based on various assumptions on the asteroid and its behaviour such as no rotation, no internal density changes, spherical shape and homogenous strength.  For more information on the project specfication, see the python notebooks: `ProjectDescription.ipynb`, `AirburstSolver.ipynb` and `DamageMapper.ipynb`.


## Results 
The Armageddon package is able to predict whether an incoming asteroid will result in an airburst, cratering or Airburst/cratering event. It will also find the energy released by the asteroid and its distance travelled along the earths surface. This information is than used to determine the impact location and estimate the damage as well as finding the affected post codes (only UK) and amount people living in the post codes. 

## Installation
To use the Armageddon package, clone this repository to your local machine and install all packages listed in requirements.txt via: 
`pip install -r requirements.txt`
Additionally, you need to install the Armageddon module via the setup.py file. If the postcode data is not downloaded, run:  python download_data.py 

## Usage 
An example usage of the code can be found in of the `example.py` or run: 
`Python example.py`

First, one needs to initialize a Planet class, this class incudes all functions that simulate the atmospheric entry. The `solve_atmospheric_entry` functions will solve the ODE system and output a Pandas Dataframe which includes the solution. The `calculate_energy` and `analyse_outcome` functions are used to calculate the kinetic energy lost per unit altitude and further statistics about the asteroid entry. 
The information can than be used to determine the blast radius, location and damage radii for the impact using the `damage_zones` functions.  The result can be visualized using the `plot_circle` function which will show the circle of the lowest damage zone on a map.  To ensure that the affected area can be evacuated, the `PostcodeLocator()` function will find all postcodes that are located in the area. An estimate for the overall population that is affected by the impact for the different damage radii is given by `postcodeLocator()`. 

## Testing
Automated tests are run after every merge or push to the main repository on github. Test can be found in the `Armageddon/test folder`. 

## Extension 
There are multiple functions that extend the functionality of the code.
-	Access to tabulated atmospheric density profiles 
-	Ability to determine asteroid parameters from observed energy deposition curves
-	Plotting the impact radii on a map 
-	A uncertainty analysis which will calculate the risk for affected areas. 




