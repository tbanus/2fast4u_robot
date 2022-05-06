#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <avr/io.h>
#include <util/delay.h> 
#include "USART.h"
#include "USART.c"
#include "pid.c"
#include "kalman.c"
// #include get_sensors




int main()
{
    pid_type heading_controller, velocity_controller;
    measurement_type sensor;
    ekf_t kalman; kalman_init(&kalman);
    DDRD=0b00000000;
    DDRB=0b00000000;
    
    while (1)
    {
        sensor=get_sensor_meas();
        kalman.z[0][0]=sensor.velocity_encoder;
        kalman.z[1][0]=sensor.heading_imu;
        kalman.z[2][0]=sensor.acceleration;
        kalman_update(&kalman);
        //velocity=kalman.x[0]
        vel=pid_update(&velocity_controller,(kalman->x[0]));
        angular_vel=pid_update(&heading_controller,(kalman->x[0]));


        

    }



    return 1;
}