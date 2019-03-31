import math

width = 1.0e-9   # meters
V = 1.6e-18      # volts
m = 9.11e-31     # kilograms
hbar = 1.0546e-34

alpha = 0.0001
delta = 0.0001
h = 1.0e-12
u = 0.0

firstrun = True
evenfunction = True

# Here we start a big loop to search for roots/solutions
while (alpha < 1.0):
   x = 0.0
   oldu = 0.0 
   u = 0.0

# Based on the odd/even nature of solution, set initial psi and dpsi
   if ( evenfunction ):
       psi = 1.0
       dpsi = 0.0
   else:
       psi = 0.0
       dpsi = 1.0

# The next loop propagates the solution of the funtion to u=2 
   while ( u < 2.0 ):
       x = x + h
       u = x / width
       du = u - oldu
       if ( abs(u) > 0.5 ):
           c = 2 * m * width**2 * V * (1.0 - alpha) / hbar**2
       else:
           c = - 2 * m * width**2 * V * alpha / hbar**2
           
       dpsi = dpsi + c * psi * du
       psi = psi + dpsi * du
       oldu = u

# Store the value of psi for the first time through the while loop
# for each new alpha tested

   if ( firstrun ):
       oldpsi = psi
       firstrun = False
   else:
# if the oldpsi and the new psi are different, we found a solution
       if ( (psi * oldpsi) < 0.0 ):
           print "*** FOUND   SOLUTION  ***,  alpha = ", str(alpha)
           firstrun = True
           evenfunction = not evenfunction
           
   alpha = alpha + delta


