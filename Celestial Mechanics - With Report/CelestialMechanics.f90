PROGRAM CelestialMech 
IMPLICIT NONE

INTEGER:: i, k, t
DOUBLE PRECISION:: time, G, dt, L, o, pi, Au, d1, d2, d3, c1, c2, c3, c4
DOUBLE PRECISION, DIMENSION(3):: Vector, dvdt, dvdt2
DOUBLE PRECISION, DIMENSION(10):: Mass, Albedo, temp, OrbitTime, X1, Y1, Z1, X2, Y2, Z2, Vx, Vy, Vz, r, avgr, rmin, rmax, e
CHARACTER(LEN=5), DIMENSION(10)::Body=(/'Sun', 'Earth', 'Moon', 'Mer', 'Venus', 'Mars', 'Jup', 'Sat', 'Ura', 'Nep'/)
CHARACTER(LEN=10) :: the_time
CHARACTER(LEN=1):: choice
G = 6.67408E-11; L=3.828E26; o=5.6705E-8; pi=4.0*atan(1.0); Au = 149598000E3
d1=1.351207192; d3=d1; d2=-1.702414384; c1=0.675603596; c4=c1; c2=-0.175603596; c3=c2; rmax=0.0; rmin=9E9999
!Defining constants for Yoshida integration, to calculate acceleration and temperature

OPEN(12, file='ExtraData.dat'); OPEN(13, file='IntialConditions.dat')
OPEN(14, file='Earth.dat');OPEN(15, file='Moon.dat');OPEN(16, file='Mercury.dat')
OPEN(17, file='Venus.dat');OPEN(18, file='Mars.dat');OPEN(19, file='Jupiter.dat')
OPEN(20, file='Saturn.dat'); OPEN(21, file='Uranus.dat');OPEN(22, file='Neptune.dat')
!Gives every planet their own file, making plotting, finding the data and formatting much easier
READ(13,*) Albedo, Mass, X1, Y1, Z1, Vx, Vy, Vz
!Reads initial conditions
CLOSE(13)
!Telling the user about the program
WRITE(6,*) 'This program is designed to calculate the positions of the planets'
WRITE(6,*) 'solar systems planets using newtons law of gravity'
WRITE(6,*) 'and various numerical integration techniques, times, and time steps'
WRITE(6,*) 'Would you like to choose a time step, integration time,'
WRITE(6,*) 'and numerical technique or use the default settings?'
WRITE(6,*) 'The default settings are 1 year, 40 seconds and the Yoshida Method'
WRITE(6,*) 'Press y to enter your own settings, n for default.'
READ(5,*) Choice


IF (choice =='y' .or. choice =='Y') THEN
  
!Gives user a way to change settings without changing the code directly
!Which can be easier in some cases
  
  WRITE(6,*) 'Please enter the numerical integration technique you would like to use'
  WRITE(6,*) 'Y for Yoshida, L for leapfrog, E for Euler'
  READ(5,*) Choice
  WRITE(6,*) 'Please enter the integration time.'
  READ(5,*) time
  WRITE(6,*) 'Please enter the time step'
  READ(5,*) dt
  
ELSE IF (choice =='n'.or. choice =='N') THEN
!This choice still allows the option to change the code directly 
 dt=40; time=31536000; choice='y'
 
END IF

WRITE(6,*) 'Coordinates are printed to files labeled with the planets name, such'
WRITE(6,*) 'as ''Earth.dat'', the coordinates are printed, X, Y, Z for 3d plotting'
WRITE(6,*) 'For other data such as orbit time and temperature please see the file'
WRITE(6,*) '''ExtraData.dat'''
!Tells user where to find the data generated

CALL DATE_AND_TIME(TIME=the_time)
WRITE(12,*) 'Start time: ', the_time
!Records the time the calculation started for later use

X2=X1; Y2=Y1; Z2=Z1


IF (choice =='y'.or. choice =='Y') THEN 
  !Yoshida Integration, if loop allows for choice and allows quick testing of all 
  !integration methods without having to change the code
    
  DO t = 0, NINT(time/dt)
    !Do loop over all time for the required number of time steps
    WRITE(14, '(e13.5,e13.5,e13.5)') X2(2), Y2(2), Z2(2);  WRITE(15, '(e13.5,e13.5,e13.5)') X2(3), Y2(3), Z2(3)
    WRITE(16, '(e13.5,e13.5,e13.5)') X2(4), Y2(4), Z2(4);  WRITE(17, '(e13.5,e13.5,e13.5)') X2(5), Y2(5), Z2(5)
    WRITE(18, '(e13.5,e13.5,e13.5)') X2(6), Y2(6), Z2(6);  WRITE(19, '(e13.5,e13.5,e13.5)') X2(7), Y2(7), Z2(7)
    WRITE(20, '(e13.5,e13.5,e13.5)') X2(8), Y2(8), Z2(8);  WRITE(21, '(e13.5,e13.5,e13.5)') X2(9), Y2(9), Z2(9)
    WRITE(22, '(e13.5,e13.5,e13.5)') X2(10), Y2(10), Z2(10)   
    !Prints X & Y coordinates to files

    X1=X2; Y1=Y2; Z1=Z2 
    !Defines new positions X1 as the final positions of last time step, X2
  
    DO k=2,10 !k= planet force acting on, i= planet causing force. 
    
      dvdt=0.0 !Resets acceleration for next planet  
      X2(k)=X1(k)+c1*Vx(K)*dt; Y2(k)=Y1(k)+c1*Vy(k)*dt; z2(k)=z1(k)+c1*Vz(k)*dt
      !X2 is the variable point, and X1 is the initial position of the planets at the 
      !start of the time step

      DO i=1,10  !calculating acceleration, dvdt, at new point  
        !Do loop sums over different planets to add their contributions to the acceleration 
        IF (i /= k)Then !Avoids calculating acceleration of a planet on itself
          r(k) = sqrt((X2(k)-X1(i))**2+(Y2(k)-Y1(i))**2+(Z2(k)-Z1(i))**2) 
          !calculating 'r' from planets new position and other planets initial positions
          !This means that the acceleration for the last planet will still be considering
          !The other planets positions as the start of the time step
          Vector=(/(X2(k)-X1(i)), (Y2(k)-Y1(i)), (Z2(k)-Z1(i))/)/r(k)  !unit Vector for direction of force
          dVdt=-(G*Mass(i)/r(k)**2)*vector+dVdt !Newtons Law  to calculate acceleration due to gravity        
        END IF  
      END DO
     
      Vx(k)=Vx(k)+d1*dvdt(1)*dt; Vy(k)=Vy(k)+d1*dvdt(2)*dt; Vz(k)=Vz(k)+d1*dvdt(3)*dt
      X2(k)=X2(k)+c2*Vx(k)*dt;  Y2(k)=Y2(k)+c2*Vy(k)*dt; Z2(k)=Z2(k)+c2*Vz(k)*dt;
      dvdt=0.0 !Acceleration reset in between as we are now considering it at a new point
           
      DO i=1,10   !calculating acceleration, dvdt, at new point
        IF (i /= k)Then 
          r(k) = sqrt((X2(k)-X1(i))**2+(Y2(k)-Y1(i))**2+(Z2(k)-Z1(i))**2) 
          Vector=(/(X2(k)-X1(i)), (Y2(k)-Y1(i)), (Z2(k)-Z1(i))/)/r(k)  
          dVdt=-(G*Mass(i)/r(k)**2)*vector+dVdt            
        END IF    
      END DO
     
      !Yoshida integration step
      Vx(k)=Vx(k)+d2*dvdt(1)*dt; Vy(k)=Vy(k)+d2*dvdt(2)*dt; Vz(k)=Vz(k)+d2*dvdt(3)*dt
      X2(k)=X2(k)+c3*Vx(k)*dt;  y2(k)=y2(k)+c3*Vy(k)*dt;  z2(k)=z2(k)+c3*Vz(k)*dt
      dvdt=0.0  !Acceleration reset in between as we are now considering it at a new point
     
      DO i=1,10   !calculating acceleration, dvdt,  at new point 
        IF (i /= k)Then 
          r(k) = sqrt((X2(k)-X1(i))**2+(Y2(k)-Y1(i))**2+(Z2(k)-Z1(i))**2) 
          Vector=(/(X2(k)-X1(i)), (Y2(k)-Y1(i)), (Z2(k)-Z1(i))/)/r(k)  
          dVdt=-(G*Mass(i)/r(k)**2)*vector+dVdt             
        END IF
      END DO
      
      !Yoshida integration step
      Vx(k)=Vx(k)+d3*dvdt(1)*dt; Vy(k)=Vy(k)+d3*dvdt(2)*dt; Vz(k)=Vz(k)+d3*dvdt(3)*dt
      x2(k)=x2(k)+c4*vx(k)*dt; y2(k)=y2(k)+c4*vy(k)*dt; z2(k)=z2(k)+c4*vz(k)*dt
      dvdt=0.0 !Acceleration reset in between as we are now considering it at a new point
      
      avgr(k)=sqrt((X2(k)-X1(1))**2+(Y2(k)-Y1(1))**2+(Z2(k)-Z1(1))**2)+avgr(k)  
      !Summing up total distance between the planet and the sun to calculate the mean late
      
      r(k) = sqrt((X2(k)-X1(1))**2+(Y2(k)-Y1(1))**2+(Z2(k)-Z1(1))**2)
      !Considers new distance (r) between planet and sun to calculate eccentricity
       
      IF (r(k) < rmin(k)) THEN
        !If r is smaller than the smallest r so far it is stored and used to calculate eccentricity    
        rmin(k)=r(k)
      ELSE IF (r(k) > rmax(k)) THEN
        rmax(k)=r(k)
        !If r is larger than the largest so far it is stored and used to calculate eccentricity
      END IF

    END DO   
  
  END DO
  
ELSE IF (choice == 'l'.or. choice =='L') THEN !Leapfrog method
  
  DO t = 0, NINT(time/dt) 
    
    !Do loop over all time for the required number of time steps
    WRITE(14, '(e13.5,e13.5,e13.5)') X2(2), Y2(2), Z2(2);  WRITE(15, '(e13.5,e13.5,e13.5)') X2(3), Y2(3), Z2(3)
    WRITE(16, '(e13.5,e13.5,e13.5)') X2(4), Y2(4), Z2(4);  WRITE(17, '(e13.5,e13.5,e13.5)') X2(5), Y2(5), Z2(5)
    WRITE(18, '(e13.5,e13.5,e13.5)') X2(6), Y2(6), Z2(6);  WRITE(19, '(e13.5,e13.5,e13.5)') X2(7), Y2(7), Z2(7)
    WRITE(20, '(e13.5,e13.5,e13.5)') X2(8), Y2(8), Z2(8);  WRITE(21, '(e13.5,e13.5,e13.5)') X2(9), Y2(9), Z2(9)
    WRITE(22, '(e13.5,e13.5,e13.5)') X2(10), Y2(10), Z2(10)   
    !Prints X & Y coordinates to files
    
    X1=X2; Y1=Y2; Z1=Z2 
    !Defines new positions X1 as the final positions of last time step, X2
    
    DO k=2,10 !k= planet force acting on, i= planet causing force. 
      dvdt=0.0; dvdt2=0.0 !Resets acceleration for when a new planet is considered
    
      DO i=1,10  !calculating acceleration of intial point
          
        IF (i /= k) THEN !if command avoids calculating acceleration of a planet on itself
          r(k) = sqrt((X2(k)-X1(i))**2+(Y2(k)-Y1(i))**2+(Z2(k)-Z1(i))**2)
          !calculating 'r' from intial positions of time step
          Vector=(/(X2(k)-X1(i)), (Y2(k)-Y1(i)), (Z2(k)-Z1(i))/)/r(k)  !unit Vector for direction of force
          dVdt=-(G*Mass(i)/r(k)**2)*vector+dVdt !Sums up all contributions to accelerations       
        END IF
         
      END DO
         
      X2(k)=X2(K)+dt*Vx(k)+0.5*(dt**2)*dvdt(1); Y2(k)=Y2(K)+dt*Vy(k)+0.5*(dt**2)*dvdt(2)
      Z2(k)=Z2(K)+dt*Vz(k)+0.5*(dt**2)*dvdt(3)
      !Calculating  new positions
          
      DO i=1,10  !calculating acceleration from new positions, dvdt2
          
        IF (i /= k) THEN
          r(k) = sqrt((X2(k)-X1(i))**2+(Y2(k)-Y1(i))**2+(Z2(k)-Z1(i))**2)
           !calculating 'r' from intial positions of time step
          Vector=(/(X2(k)-X1(i)), (Y2(k)-Y1(i)), (Z2(k)-Z1(i))/)/r(k)  !unit Vector for direction of force
          dVdt2=-(G*Mass(i)/r(k)**2)*vector+dVdt2 
        END IF
       
      END DO
           
      Vx(k)=Vx(k)+0.5*dt*(dvdt(1)+dvdt2(1)); Vy(k)=Vy(k)+0.5*dt*(dvdt(2)+dvdt2(2))
      Vz(k)=Vz(k)+0.5*dt*(dvdt(3)+dvdt2(3))
      !New velocity calculated using the average acceleration between intial and final positions
      
      avgr(k)=sqrt((X2(k)-X1(1))**2+(Y2(k)-Y1(1))**2+(Z2(k)-Z1(1))**2)+avgr(k)  
      !Summing up total distance between the planet and the sun to calculate the mean late
      
      r(k) = sqrt((X2(k)-X1(1))**2+(Y2(k)-Y1(1))**2+(Z2(k)-Z1(1))**2)
      !Calculates new distance (r) between planet and sun to calculate eccentricity
       
      IF (r(k) < rmin(k)) THEN
        !If r is smaller than the smallest r so far it is stored and used to calculate eccentricity    
        rmin(k)=r(k)
      ELSE IF (r(k) > rmax(k)) THEN
        rmax(k)=r(k)
        !If r is larger than the largest so far it is stored and used to calculate eccentricity
      END IF
    
    END DO
    
  END DO
       
ELSE IF (choice == 'e'.or. choice =='E') THEN !Euler Method
  
  DO t = 0, NINT(time/dt)  
    
    !Do loop over all time for the required number of time steps
    WRITE(14, '(e13.5,e13.5,e13.5)') X2(2), Y2(2), Z2(2);  WRITE(15, '(e13.5,e13.5,e13.5)') X2(3), Y2(3), Z2(3)
    WRITE(16, '(e13.5,e13.5,e13.5)') X2(4), Y2(4), Z2(4);  WRITE(17, '(e13.5,e13.5,e13.5)') X2(5), Y2(5), Z2(5)
    WRITE(18, '(e13.5,e13.5,e13.5)') X2(6), Y2(6), Z2(6);  WRITE(19, '(e13.5,e13.5,e13.5)') X2(7), Y2(7), Z2(7)
    WRITE(20, '(e13.5,e13.5,e13.5)') X2(8), Y2(8), Z2(8);  WRITE(21, '(e13.5,e13.5,e13.5)') X2(9), Y2(9), Z2(9)
    WRITE(22, '(e13.5,e13.5,e13.5)') X2(10), Y2(10), Z2(10)   
    !Prints X & Y coordinates to files
    
    X1=X2; Y1=Y2; Z1=Z2
    
    DO k=2,10 !k= planet force acting on, i= planet causing force. 
      dvdt=0
      
      DO i=1,10  !calculating acceleration
        
        IF (i /= k)Then !Avoids calculating acceleration of a planet on itself
          r(k) = sqrt((X2(k)-X1(i))**2+(Y2(k)-Y1(i))**2+(Z2(k)-Z1(i))**2) 
          !calculating 'r' from intial positions of time step
          Vector=(/(X2(k)-X1(i)), (Y2(k)-Y1(i)), (Z2(k)-Z1(i))/)/r(k)  !unit Vector for direction of force
          dVdt=-(G*Mass(i)/r(k)**2)*vector+dVdt !Newtons law of gravitation    
        END IF   
        
      END DO
      
      X2(k) = Vx(k)*dt + X2(k); Y2(k) = Vy(k)*dt + Y2(k); Z2(k) = Vz(k)*dt + Z2(k) 
      !Assume velocity is constant over time step and use this to calculate new positions
      Vx(K) = dVdt(1)*dt + Vx(k); Vy(k) = dVdt(2)*dt + Vy(k); Vz(k) = dVdt(3)*dt + Vz(k) 
      !Assume acceleration is constant over time step and use this to calculate new velocities
      
      avgr(k)=sqrt((X2(k)-X1(1))**2+(Y2(k)-Y1(1))**2+(Z2(k)-Z1(1))**2)+avgr(k)  
      !Summing up total distance between the planet and the sun to calculate the mean late
      
      r(k) = sqrt((X2(k)-X1(1))**2+(Y2(k)-Y1(1))**2+(Z2(k)-Z1(1))**2)
      !Considers new distance (r) between planet and sun to calculate eccentricity
       
      IF (r(k) < rmin(k)) THEN
        !If r is smaller than the smallest r so far it is stored and used to calculate eccentricity    
        rmin(k)=r(k)
      ELSE IF (r(k) > rmax(k)) THEN
        rmax(k)=r(k)
        !If r is larger than the largest so far it is stored and used to calculate eccentricity
      END IF
   
    END DO
    
  END DO
  
END IF
     
avgr = avgr/NINT(time/dt)

Temp = ((L*(1-albedo))/(4.0*pi*o*(avgr)**2))**0.25 -273.15
!Calculating surface temp. from flux, albedo, and average distance of planet to sun
OrbitTime = (2*pi*((avgr**3)/(G*Mass(1)))**0.5)/31536000
!Calculating orbit time from average distance of planet to sun and the mass of the sun, in years
e=(rmax-rmin)/(rmax+rmin) 
!Calculating eccentricity of orbit

WRITE(12,*) 'Body, Avg. dist. to sun (AU Orbital Period, Temp, eccentricity'
!Printing titles to make data easier to interpret.

DO k = 2, 10
  
  WRITE(12,'(a5, f6.3, f8.3, f9.3, f8.4)') Body(k), avgr(k)/au, orbitTime(k), TEMP(K), e(k)
  !Printing data to file

END DO

WRITE(12,'(a10, f6.1)') 'Step time:', dt
WRITE(12,'(a17, f13.2)') 'Integration time:', time
WRITE(12,*) 'Number of steps:', NINT(time/dt)
WRITE(12,'(a11, a2)') 'Method used', Choice
CALL DATE_AND_TIME(TIME=the_time)
WRITE(12,*) 'Final time: ', the_time
!Recording all the final data in a file

CLOSE(12); CLOSE(14); CLOSE(15); CLOSE(16); CLOSE(17); CLOSE(18); CLOSE(19); CLOSE(20)
CLOSE(21); CLOSE(22)

END PROGRAM 

