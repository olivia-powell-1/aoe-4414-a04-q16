# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  Will convert from ecef coordinate frame to sez coordinate frame
# Parameters:
#  o_x_km: sez x position of ecef origin
#  o_y_km: sez y position of ecef origin 
#  o_z_km: sez z position of ecef origin 
#  x_km: ecef x position
#  y_km: ecef y position 
#  z_km: ecef z position
# Output:
#  s_km: s position
#  e_km: e position
#  z_km: z position
# Written by Olivia Powell
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
# e.g., import math # math module
import sys # argv
import math # math ops
import numpy as np # matrix function

# "constants"
r_e_km = 6378.137 # redius of earth
e_e = 0.081819221456 # eccentricity of earth

# helper functions

## calculated denominator
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

## convert ecef to llh coordinate system
def calc_ecef_to_llh(r_x_km:float, r_y_km:float, r_z_km:float):
   # calculate longitude
   lon_rad = math.atan(r_y_km,r_x_km)
   lon_deg = lon_rad*180.0/math.pi
   # initialize lat_rad, r_lon_km, r_z_km
   lat_rad = math.asin(r_z_km/math.sqrt(r_x_km**2+r_y_km**2+r_z_km**2))
   r_lon_km = math.sqrt(r_x_km**2+r_y_km**2)
   prev_lat_rad = float('nan')
   # iteratively find latitude
   c_E = float('nan')
   count = 0
   while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
     denom = calc_denom(e_e,lat_rad)
     c_E = r_e_km/denom
     prev_lat_rad = lat_rad
     lat_rad = math.atan((r_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
     count = count+1
  # calculate hae
   hae_km = r_lon_km/math.cos(lat_rad)-c_E
   return [lat_rad, lon_rad, hae_km]

# initialize script arguments
o_x_km = float('nan')
o_y_km = float('nan')
o_z_km = float('nan')
x_km = float('nan')
y_km = float('nan')
z_km = float('nan')

# parse script arguments
if len(sys.argv)==7:
   o_x_km = float(sys.argv[1])
   o_y_km = float(sys.argv[2])
   o_z_km = float(sys.argv[3])
   x_km = float(sys.argv[4])
   y_km = float(sys.argv[5])
   z_km = float(sys.argv[6])
   ...
else:
   print(\
    'Usage: '\
    'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
   )
   exit()

# write script below this line
# first we must convert from ecef to llh
r_ecef = calc_ecef_to_llh(o_x_km, o_y_km, o_z_km)
r_f = np.array([x_km-o_x_km],[y_km-o_y_km],[z_km-o_z_km])
o_lat_rad = r_ecef[0]
o_lon_rad = r_ecef[1]
hae = r_ecef[2]

# define y axis rotation matrix
ry = np.array([[math.sin(o_lat_rad), 0, -math.cos(o_lat_rad)], [0, 1, 0], [math.cos(o_lat_rad), 0, math.sin(o_lat_rad)]])

# define z axis rotation matrix
rz = np.array([[math.cos(o_lon_rad), math.sin(o_lon_rad), 0], [-math.sin(o_lon_rad), math.cos(o_lon_rad), 0], [0, 0, 1]])

# calculate sez vector
y_rot = np.matmul(ry, r_f) # multiply y axis rotation matrix by r_f vector
z_rot = np.matmul(rz, y_rot) # multiply z axis rotation matrix by results of y axis rotation above
s_km = r_f[0] + z_rot[0]
e_km = r_f[1] + z_rot[1]
z_km = r_f[2] + z_rot[2]

# print
print(s_km)
print(e_km)
print(z_km)