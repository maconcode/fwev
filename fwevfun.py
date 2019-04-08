import math

width = 1.0e-9   # meters
V = 1.6e-18      # volts
m = 9.11e-31     # kilograms
hbar = 1.0546e-34

alpha = 0.02
h = 1.0e-14
u = 0.0
root=0


firstrun = True
evenfunction = True

oldalpha = 0.0 
psi = 0.0
dpsi = 0.0

# Here we start a big loop to search for roots/solutions
while (alpha < 0.04):
   x = 0.0
   oldu = 0.0 
   u = 0.0
   delta = 0.001
   olddelta = delta
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
       oldalpha = alpha
   else:
# if the oldpsi and the new psi are different, we found a solution
       if ( (psi * oldpsi) < 0.0 ):
           print ("*** FOUND INITIAL  SOLUTION  ***,  alpha =  %18.16f" % alpha)
           # Now put the wavefunction data in a file
           root = root + 1
           filename = "wavefunction_" + str(root) + ".dat"
           filehandle = open(filename,"w")

          # improve the solution
           iteration = 0
           oldalpha = alpha
           olddelta = delta
           newh = h
           while (iteration < 15):
             alpha = alpha - 10*delta 
             delta = delta / 10             
             #print ("   current alpha is %18.16f" % alpha)
             solutionFound = False
             while ( not solutionFound  ):
                alpha = alpha + delta;
                if ( evenfunction ):
                  psi = 1.0
                  dpsi = 0.0
                else:
                  psi = 0.0
                  dpsi = 1.0

# The next loop propagates the solution of the funtion to u=2 
                u = 0.0
                x = 0.0
                oldu = 0.0 
                while ( u < 2.0):
                   x = x + newh
                   u = x / width
                   du = u - oldu
                   if ( abs(u) > 0.5 ):
                     c = 2 * m * width**2 * V * (1.0 - alpha) / hbar**2
                   else:
                     c = - 2 * m * width**2 * V * alpha / hbar**2
           
                   dpsi = dpsi + c * psi * du
                   psi = psi + dpsi * du
                   oldu = u

                print "                   alpha = " + str(alpha) + "  psi = " + str(psi) + " oldpsi = " + str(oldpsi)
                if ((psi * oldpsi) < 0 ):
                   solutionFound = True
                   iteration = iteration + 1 
                   #oldpsi = psi
                   newh = newh * .75 
                   print( "    -->  improved alpha is %18.13f on iteration %i"  % (alpha, iteration)) 
                   print "h = " + str(newh)       

           if ( evenfunction ):
              psi = 1.0
              dpsi = 0.0
           else:
              psi = 0.0
              dpsi = 1.0

           u=0.0
           x=0.0
           oldu=-h/width

           filehandle.write(str(u) + " " + str(psi) + "\n");

           print "using alpha = " + str(alpha) + " to generate wavefunction plot"
           while ( u < 2.0 ):
              x = x + newh
              u = x / width
              du = u - oldu
              if ( abs(u) > 0.5 ):
                 c = 2 * m * width**2 * V * (1.0 - alpha) / hbar**2
              else:
                 c = - 2 * m * width**2 * V * alpha / hbar**2
           
              dpsi = dpsi + c * psi * du
              psi = psi + dpsi * du
              filehandle.write(str(u) + " " + str(psi) + "\n");
             
              oldu = u

           filehandle.close()

           firstrun = True
           evenfunction = not evenfunction
   
   print "Checking alpha = " + str(alpha)        
   alpha = alpha + delta


