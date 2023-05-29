
#define wgs84_a 6378137.0  // wgs84.SemimajorAxis
#define wgs84_b 6356752.314245 // wgs84.SemiminorAxis
#define wgs84_f 0.003352810664747 // wgs84.Flattening : f = (a-b)/b
#define wgs84_e -0.006739496742276 // first numerical eccentricity of the Earth : e = f*(2-f)

#define deg2rad 0.017453292519943
#define rad2deg 57.295779513082320


// Declinação magnetica da Terra em Sidney na data: 12/02/2016 [rad] (https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml)
#define mag_dec 12.34*deg2rad

// Sampling time
#define dt 0.05

// GPS sampling time
#define dt_gps 0.05

// Local gravity (m/s^2)
#define g -9.81

// Data vector size
#define N 70524