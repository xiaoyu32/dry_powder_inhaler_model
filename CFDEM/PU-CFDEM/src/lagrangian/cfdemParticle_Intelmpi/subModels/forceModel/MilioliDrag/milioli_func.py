#!/usr/bin/python

import subprocess

rhos = 1489.86
rho=1.2
ds=78.e-06
nuf=1.5e-05
grav=9.81

# calc tau_t,v_t
taup = rhos/rho*ds**2/(18.*nuf)
v_t = taup * grav

# 
alps=float

# h_env1 
if alps <= 0.1:
	func_h_env1 == 0.5643*(1.0 + alps)*(alps^0.15)/(0.5766*(alps^0.3) + 0.1997)
elif alps > 0.1 & alps <= 0.54:
	func_h_env1 = 0.8428 + 0.6393*alps - 0.6743*alps*alps
elif alps > 0.54 & alps <= 0.65:
	func_h_env1 = 0.4099*(0.65-alps)^0.25/(alps^(-0.25) - 0.9281)
elif alps >  0.65:
	func_h_env1 = 0.0

# h_env2
func_h_env2 = 0.8428 + 0.6393*alps - 0.6743*alps*alps

# h_env
#  func_h_env = min(func_h_env1,func_h_env2)

# h_1
func_h_1 = alps*(1.6*magUr/v_t + 4.0)/(7.9*magUr/v_t + 0.08) + (0.9394 - 0.22/(0.6*magUr/v_t + 0.01))

# func_f_inf
func_f_inf = 0.882*(2.145 - 7.8*pow(magUr/v_t,1.8)/(7.746*(magUr/v_t)^1.8 + 0.5586) )

if func_h_1 > 0:
	func_h_lin = func_f_inf*func_h_1
else:
	func_h_lin = 0

# Solid volume dependency of model func_h
func_h = min(func_h_env,func_h_lin)

# Infinite grid resolution
func_f = -1.

# calc drag coefficient 
drag_correction = Kii*func_h*func_f
  
  
proc = subprocess.Popen(['gnuplot','-p'], 
                        shell=True,
                        stdin=subprocess.PIPE,
                        )
proc.stdin.write('set xrange [0:10]; set yrange [-2:2]\n')
proc.stdin.write('plot sin(x)\n')
proc.stdin.write('quit\n')
