typedef struct 
{
    /* data */
    double kp;
    double ki;
    double kd;
    double last_output;
    double last_time;
    double last_input;
    double setpoint;
    double min_output;
    double max_output;
    double p_term;
    double i_term;
    double d_term;

} pid_type;


double pid_update(pid_type * pid, double current_input)
{   
    double now,dt,error,d_input,output;
    now=   clock();
    if (now - pid->last_time) dt = now -pid->last_time; else dt= 1e-16;
    if (dt<0) return 1;
    error=pid->setpoint - current_input;
    d_input=current_input - pid->last_input;
    pid->p_term=pid->kp*error;

    pid->i_term+= pid->ki * error* dt;

    pid->d_term= -1 * pid->kd * d_input / dt;

    output= pid->p_term + pid->i_term + pid->d_term;

    pid->last_output=output;
    pid->last_time=now;
    pid->last_input=current_input;


    if (output> pid->max_output) return pid->max_output;
    else if (output < pid->min_output) return pid-> min_output;
    else return output;
}

void pid_reset_i_term( pid_type * pid )
{
    pid->i_term=0;
}



pid_type init_pid(double kp, double ki, double kd, double setpoint)
{
    pid_type pid;
    pid.kp=kp;
    pid.ki=ki;
    pid.kd=kd;
    pid.setpoint=setpoint;

    return pid;

}


int main()
{

    pid_type pid;
    double r;
    r=pid_update(&pid, 10 );
    printf("%f",r);

    return 1;
}