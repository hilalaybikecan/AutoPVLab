#include <CNCShield.h>

#define NO_OF_STEPS               200
#define SLEEP_BETWEEN_STEPS_MS    15

CNCShield cnc_shield;
StepperMotor *motor = cnc_shield.get_motor(1);

void setup()
{
  Serial.begin(9600);
  cnc_shield.begin();
  cnc_shield.enable();

 /* motor->set_dir(COUNTER);
  for (int i = 0; i < NO_OF_STEPS; i++) {
    motor->step();
    delay(SLEEP_BETWEEN_STEPS_MS);
  }*/

  //cnc_shield.disable();
  Serial.println("init");
  motor->set_speed(100);
}

int i = 0;
int i2 = 0;
int i3 = 0;
int i4 = 0;
int result = 0;
int pos = 0;
int steps = 0;
int sign = 1;


void loop() {


 
  sign = 1;
 while (Serial.available() > 0) 
 {
    i = Serial.read();
          Serial.println(i);

    delay(50);
    while(i > 57)
    {i = i - 48;}
      {
      if ((i >= 48) && (i <= 57))
        result = result * 10 + i - 48;

      i2 = 1;
      }
      if (i == 45)
          sign = -1;
      if (i == 18)
          sign = -1;
     
      if (i == 95)
      {
        sign = -1;
      }
    
 }
    
    steps = result*sign; //calculate step from previus position
    result = 0;
    if(steps == 0)
    {
      Serial.println("won't move");
      delay(500);
      // See if motor is still turning
      // Serial.println("done");

      
    }
    if (steps > 0)
    {

      Serial.println("back");
      Serial.println(steps);
      motor->set_dir(COUNTER);
      motor->step(steps);
      steps = 0;
    }
    if (steps < 0)
    {
      steps = steps*(-1);    

      
      Serial.println("fwd");
      Serial.println(steps);
      motor->set_dir(CLOCKWISE);
      motor->step(steps);
      steps = 0;
    }
 }
