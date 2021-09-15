#Pennycuick calculations from http://folk.uib.no/nzlkj/bio241/
#source('http://folk.uib.no/nzlkj/bio241/power.curve.r')

#m = mass in kg
#wingspan = span in m
pennypower <- function(m, wingspan, altitude, airspeed)
{
  
  #Frontal area (area of the cross section of the bird's body (m^2)):
  Ab <- 0.00813*m^0.666 #This formula is described on page 50 of Modelling the Flying Bird by C.J. Pennycuick (2008), ISBN: 978-0-12-374299-8.
  
  #Acceleration of gravity
  g <- 9.81
  
  #Drag Coefficient of the body:
  Cb <- 0.366
  
  #Area of the wing disk (m^2) based on given wingspan (m):
  Ad <- pi*(wingspan/2)^2
  
  #Air pressure depending on given altitude and temp:
  pressure <- 101325*(1-((0.0065*altitude)/288.15))^((9.80665-0.0289644)/(8.31447-0.0065))
  temperature <- 25+273.15
  #The corresponding air density:
  p <- (pressure*0.0289644)/(8.31447*temperature)
  
  #Velocity range based on given values (m/s):
  V <- airspeed
  
  #Induced Power:
  Pind <- 1.2*(m*g)^2/(2*Ad*V*p)
  
  #Parasite Power:
  Ppar <- 0.5*p*Ab*Cb*V^3
  
  #Profile Power:
  Sum.Pind.and.Ppar <- Pind+Ppar
  Ppro <- min(Sum.Pind.and.Ppar)*1.2
  
  #Mechanical Power:
  Pmech <- Pind+Ppar+Ppro
  return(Pmech)
  
  #Making a functional expression of the mechanical power curve:
  # Cind <- 1.2*(m*g)^2/(2*Ad*p)
  # Cpar <- 0.5*p*Ab*Cb
  # MP <- function(x)
  # {
  #   Cind/x + Cpar*x^3 + Ppro
  # }
  
  #Finding Vmp:
  #x <- seq(minvel,maxvel,0.0001)
  #y <- MP(x)
  #Vmp <- x[y==min(y)]
  
  #Finding Vmr:
  # yy <- MP(x)/x
  # Vmr <- x[yy==min(yy)]
}
