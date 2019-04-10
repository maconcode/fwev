import math
import sys
import time
start_time = time.time()

# The mass of the electron and hbar are global variables used in all parts of the
# code so they are defined here
m = 9.11e-31     # kilograms
hbar = 1.0546e-34

# This is the function to return the value of the wavefunction at 2*u
def psiAtTwoU( width, V, h, alpha, evenfunction ):

   psi = 0.0    # I think these have to be here so things scope properly
   dpsi = 0.0

# Based on the odd/even nature of solution, set initial psi and dpsi
   if ( evenfunction ):
       psi = 1.0
       dpsi = 0.0
   else:
       psi = 0.0
       dpsi = 1.0
   x = 0.0
   oldu = -h/width
   u = 0.0
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
   return psi    # end of function

################################
### Main program starts here ###
################################

### These should be the only two lines you have to really modify

width = 5.0e-10   # potential width in meters
V     = 1.6e-18   # potential depth in volts


alpha = 0.0
h = 1.0e-15
u = 0.0
root=0

firstrun = True
evenfunction = True

oldalpha = 0.0 
psi = 0.0
dpsi = 0.0
startAlpha = alpha

#This is simply a heading for the program when it runs
print " "
print " ==============================================================="
print " ====  F I N I T E   W E L L   E I G E N F U N C T I O N S  ===="
print " ==============================================================="
print " "

# Here we start a big loop to search for roots/solutions
while (alpha < 1.0):
   x = 0.0
   oldu = 0.0 
   u = 0.0
   delta = 0.001
   olddelta = delta

   if ( evenfunction ):
       psi = 1.0
       dpsi = 0.0
   else:
       psi = 0.0
       dpsi = 1.0


 
# Store the value of psi for the first time through the while loop
# for each new alpha tested

   if ( firstrun ):
       oldpsi = psiAtTwoU( width, V, h, alpha, evenfunction )
       firstrun = False
       oldalpha = alpha
       print "*** USING LINEAR SEARCH TO LOCATE NEIGHBORHOOD OF NEXT EIGENVALUE  ***"
   else:
# if the oldpsi and the new psi are different, we found a solution
       psi = psiAtTwoU( width, V, h, alpha, evenfunction )
       if ( (psi * oldpsi) < 0.0 ):
           print ("*** FOUND INITIAL  SOLUTION  ***,  alpha =  %18.16f" % alpha)
           print  "    --> SWITCHING TO BISECTION ALGORITHM TO IMPROVE SOLUTION BEFORE PLOTTING EIGENFUNCTION <--"

           # start the bisection algorithm.  For this to work properly
           # the size of h has to be much smaller 

           A = priorAlpha
           B = alpha
           fatA = priorPsi
           fatB = psi
        
           iteration = 0 
           bisecth = h / 100
           while ( abs(psi) > 0.00001 and iteration < 65 ):
              iteration = iteration + 1
              half = (A+B)/2
              test = psiAtTwoU( width, V, bisecth, half, evenfunction )

              if ( test*fatA < 0.0):
                 B = half 
                 fatB = test
              else:
                 A = half 
                 fatA = test
              
              alpha = (A+B)/2 
              psi = psiAtTwoU( width, V, bisecth, alpha, evenfunction )
              print( "    -->  Improved alpha is %18.16f on iteration %i with psi = %f"  % (alpha, iteration, psi))
             
           if ( iteration == 65 ):
              print " *** COULD NOT CONVERGE TO SOLUTION AFTER 64 BISECTIONS ***"
              print " *** TERMINATING PROGRAM EXECUTION ***"
              sys.exit() 

           print "    -- CONVERGENCE CRITERIA MET IN " + str(iteration) + " ITERATIONS ---"

           # Now put the wavefunction data in a file
          
           root = root + 1
           filename = "wavefunction_" + str(root) + ".dat"
           filehandle = open(filename,"w")


           if ( evenfunction ):
             psi = 1.0
             dpsi = 0.0
           else:
             psi = 0.0
             dpsi = 1.0

           print "       Writing eigenfunction plot data to file named " + filename
           x = 0.0
           u = 0.0
           filehandle.write(str(u) + " " + str(psi) + "\n");
           bisecth = h / 10
           oldu = -bisecth/width
           pointcount = 0
           while ( u < 2.0 ):
              x = x + bisecth 
              u = x / width
              du = u - oldu
              if ( abs(u) > 0.5 ):
                 c = 2 * m * width**2 * V * (1.0 - alpha) / hbar**2
              else:
                 c = - 2 * m * width**2 * V * alpha / hbar**2
           
              dpsi = dpsi + c * psi * du
              psi = psi + dpsi * du
              if ( pointcount % 10 == 0):
                 filehandle.write(str(u) + " " + str(psi) + "\n");
             
              oldu = u
              pointcount = pointcount + 1

           filehandle.close()
           print  "       Plot data successfully written to " + filename

           firstrun = True
           evenfunction = not evenfunction
       else:
           priorAlpha = alpha
           priorPsi = psi 
            
 
   #print "Checking alpha = " + str(alpha)        
   alpha = alpha + delta

print " "
print " *** LINEAR SEARCH DONE ***" 
print " "
print " *** PROGRAM COMPLETE ***"

elapsed_time = time.time() - start_time
print time.strftime(" *** RUN TIME: %H hours, %M minutes, %S seconds", time.gmtime(elapsed_time))
print " "
