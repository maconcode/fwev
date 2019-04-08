import math

m = 9.11e-31     # kilograms
hbar = 1.0546e-34

def psiAtTwoU( width, V, h, alpha, evenfunction ):


# Based on the odd/even nature of solution, set initial psi and dpsi
   if ( evenfunction ):
       psi = 1.0
       dpsi = 0.0
   else:
       psi = 0.0
       dpsi = 1.0

   x = 0.0
   oldu = 0.0
#   oldu = -h/width
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
   return psi

   
width = 1.0e-9   # meters
V = 1.6e-18      # volts

alpha = 0.0
h = 1.0e-14
u = 0.0
root=0

firstrun = True
evenfunction = True

oldalpha = 0.0 
psi = 0.0
dpsi = 0.0
startAlpha = alpha

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
       oldpsi = psi 
       firstrun = False
       oldalpha = alpha
   else:
# if the oldpsi and the new psi are different, we found a solution
       psi = psiAtTwoU( width, V, h, alpha, evenfunction )
       if ( (psi * oldpsi) < 0.0 ):
           print ("*** FOUND INITIAL  SOLUTION  ***,  alpha =  %18.16f" % alpha)

           # improve the solution via 10 bisections

           A = priorAlpha
           B = alpha
           fatA = priorPsi
           fatB = psi
        
           iteration = 0 
           bisecth = h / 10
           while ( abs(psi) > 0.0001 ):
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
              print( "    -->  improved alpha is %18.16f on iteration %i with psi = %f"  % (alpha, iteration, psi))
              

           # Now put the wavefunction data in a file

          
           root = root + 1
           filename = "wavefunction_" + str(root) + ".dat"
           filehandle = open(filename,"w")


           
           filehandle.write(str(u) + " " + str(psi) + "\n");

           print "using alpha = " + str(alpha) + " to generate wavefunction plot"
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
              filehandle.write(str(u) + " " + str(psi) + "\n");
             
              oldu = u

           filehandle.close()

           firstrun = True
           evenfunction = not evenfunction
       else:
           priorAlpha = alpha
           priorPsi = psi 
            
 
   print "Checking alpha = " + str(alpha)        
   alpha = alpha + delta


