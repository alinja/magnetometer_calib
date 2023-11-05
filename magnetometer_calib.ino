#include <Wire.h>
#include <Adafruit_LIS3MDL.h>
#include <Adafruit_Sensor.h>

int led = 13; //led pin

//
//Sensor: Adafruit_LIS3MDL
//
Adafruit_LIS3MDL lis3mdl;
void sensor_setup() {
  if (! lis3mdl.begin_I2C()) {
    Serial.println("Failed to find LIS3MDL chip");
    while (1) { delay(10); }
  }
  lis3mdl.setPerformanceMode(LIS3MDL_ULTRAHIGHMODE);
  lis3mdl.setOperationMode(LIS3MDL_CONTINUOUSMODE);
  lis3mdl.setDataRate(LIS3MDL_DATARATE_10_HZ);
  lis3mdl.setRange(LIS3MDL_RANGE_4_GAUSS);
}

//Get samples in uTeslas
void sensor_get(float &x, float &y, float &z) {
  lis3mdl.read();
  sensors_event_t event; 
  lis3mdl.getEvent(&event);
  x = event.magnetic.x;
  y = event.magnetic.y;
  z = event.magnetic.z;
}

//
// Floating-point implementation of sphere fitting based on https://arxiv.org/ftp/arxiv/papers/1506/1506.02776.pdf
//
class magnetometer_calib {
public:
  magnetometer_calib();
  void init(void);
  void accumulate(float x, float y, float z);
  void finalize(int n, float &x, float &y, float &z, float &r);

private:
  float Sx;   float Sy; 	float Sz;
  float Sxx;  float Syy;
  float Szz;  float Sxy;
  float Sxz;  float Syz;
  float Sxxx; float Syyy;
  float Szzz; float Sxyy;
  float Sxzz; float Sxxy;
  float Sxxz; float Syyz;
  float Syzz;
};

magnetometer_calib::magnetometer_calib()
{
}
void magnetometer_calib::init(void)
{
  Sx = 0;   Sy = 0;  Sz = 0;
  Sxx = 0;  Syy = 0;
  Szz = 0;  Sxy = 0;
  Sxz = 0;  Syz = 0;
  Sxxx = 0; Syyy = 0;
  Szzz = 0; Sxyy = 0;
  Sxzz = 0; Sxxy = 0;
  Sxxz = 0; Syyz = 0;
  Syzz = 0;
}

void magnetometer_calib::accumulate(float x, float y, float z)
{
	Sx   = Sx+x;         Sy   = Sy+y; Sz = Sz+z;
	Sxx  = Sxx+x*x;     Syy  = Syy+y*y;
	Szz  = Szz+z*z;     Sxy  = Sxy+x*y;
	Sxz  = Sxz+x*z;     Syz  = Syz+y*z;
	Sxxx = Sxxx+x*x*x; Syyy = Syyy+y*y*y;
	Szzz = Szzz+z*z*z; Sxyy = Sxyy+x*y*y;
	Sxzz = Sxzz+x*z*z; Sxxy = Sxxy+x*x*y;
	Sxxz = Sxxz+x*x*z; Syyz = Syyz+y*y*z;
	Syzz = Syzz+y*z*z;
}

void magnetometer_calib::finalize(int N, float &x, float &y, float &z, float &r)
{
	float A1 = Sxx +Syy +Szz;

	float a = 2*Sx*Sx-2*N*Sxx;
	float b = 2*Sx*Sy-2*N*Sxy;
	float c = 2*Sx*Sz-2*N*Sxz;
	float d = -N*(Sxxx +Sxyy +Sxzz)+A1*Sx;
	float e = 2*Sx*Sy-2*N*Sxy;
	float f = 2*Sy*Sy-2*N*Syy;
	float g = 2*Sy*Sz-2*N*Syz;
	float h = -N*(Sxxy +Syyy +Syzz)+A1*Sy;
	float j = 2*Sx*Sz-2*N*Sxz;
	float k = 2*Sy*Sz-2*N*Syz;
	float l = 2*Sz*Sz-2*N*Szz;
	float m = -N*(Sxxz +Syyz + Szzz)+A1*Sz;

	float delta = a*(f*l - g*k)-e*(b*l-c*k) + j*(b*g-c*f);

	x = (d*(f*l-g*k) -h*(b*l-c*k) +m*(b*g-c*f))/delta;
	y = (a*(h*l-m*g) -e*(d*l-m*c) +j*(d*g-h*c))/delta;
	z = (a*(f*m-h*k) -e*(b*m-d*k) +j*(b*h-d*f))/delta;
	r = sqrt(x*x+y*y+z*z+(A1-2*(x*Sx+y*Sy+z*Sz))/N);
}

//
// Samples collection
//

// Total samples to collect, minimum is about 8 samples
#define SAMPLES_TOTAL  (16)

//Experimentally found reasonable spacing to require between samples, depends on strength of magnetic field!
#define SAMPLES_MIN_DISTANCE_UT  (45.0*6/SAMPLES_TOTAL)

void print_xyz(float x, float y, float z) {
  Serial.print(""); Serial.print(x);
  Serial.print(","); Serial.print(y); 
  Serial.print(","); Serial.print(z);
  Serial.println("");
}

int check_distance_ok(float x1, float y1, float z1, float x2, float y2, float z2) {
  float dist = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
  if (dist > SAMPLES_MIN_DISTANCE_UT*SAMPLES_MIN_DISTANCE_UT) return 1;
  else return 0;
}

int check_min_distance(float x, float y, float z, float *vec_x, float *vec_y, float *vec_z) {
  for(int i=0; i<SAMPLES_TOTAL; i++) {
    if( !(vec_x[i] == 0 && vec_y[i] == 0 && vec_z[i] == 0) ) {
      if(!check_distance_ok(x, y, z, vec_x[i], vec_y[i], vec_z[i])) return -1;
    }
    if(vec_x[i] == 0 && vec_y[i] == 0 && vec_z[i] == 0 ) {
      vec_x[i] = x;
      vec_y[i] = y;
      vec_z[i] = z;
      return i;
    }
  }
  return -2;
}

void calibrate() {
  magnetometer_calib cal;
  float x_vec[SAMPLES_TOTAL];
  float y_vec[SAMPLES_TOTAL];
  float z_vec[SAMPLES_TOTAL];
  for(int i=0;i<SAMPLES_TOTAL;i++) {
    //Init to 0,0,0 meaning empty slot. If we really get a measurement of 0,0,0 it will just get discarded.
    x_vec[i]=0; y_vec[i]=0; z_vec[i]=0;
  }

  cal.init();
  int round=0;
  int sample_c=0;
  int cal_ready=0;
  while(1) {

    delay(100);

    float x, y, z;
    sensor_get(x, y, z);

    if(check_min_distance(x, y, z, x_vec, y_vec, z_vec) >= 0) {
      print_xyz(x, y, z);
      cal.accumulate(x,y,z);
      sample_c++;
    } else if(cal_ready) {
      //We could keep on accumulating points forever, as long as the sensor is rotated to avoid weighting one direction
      //print_xyz(x, y, z);
      cal.accumulate(x,y,z);
      sample_c++;
    }

    int zeros=0;
    for(int i=0;i<SAMPLES_TOTAL;i++)
      if(x_vec[i] == 0 && y_vec[i] == 0 && z_vec[i] == 0 ) zeros++;
    //Serial.print("Zeros left: "); Serial.println(zeros);

    if(zeros == 0) {
      float xc, yc, zc, r;
      cal.finalize(sample_c, xc, yc, zc, r);

      Serial.println("CLIBRATION READY:");
      Serial.print(""); Serial.print(xc);
      Serial.print(","); Serial.print(yc); 
      Serial.print(","); Serial.print(zc);
      Serial.print(","); Serial.print(r);
      Serial.println("");
      Serial.println("");
      cal_ready = 1;
      break;
    }
  }
}



void setup() {
  pinMode(led, OUTPUT);
  Serial.begin(115200);
  sensor_setup();
  delay(3000);
  Serial.println("magnetometer_calib ready to start!");
}

void loop() {
  digitalWrite(led, HIGH); 
  calibrate();
  digitalWrite(led, LOW);
  
  delay(5000);
  delay(5000);
  delay(5000);
}

