public class CircularMotion {
    /*
        Centripetal Acceleration

        given a circle with speed, v
        and radius r
     */
    //given in meters
    static double getCycleDistance(double radius){
        return 2.0 * Math.PI * radius;
    }
    //time = time taken for one cycle to execute
    //given in m/s
    static double getVelocity(double area, double time){
        if(time != 0){
            return (area / time);
        } else {
            return 0;
        }
    }
    /*
        Period, T
            the amount of time it takes to do one cycle
    */
    static double getPeriod(double frequency){
        return (1 / frequency);
    }
    /*
        Angular speed
            change in angle with respect to time (in seconds)
            Number of radians traveled per second
     */
    //given in radians per second
    static double getAngularSpeed(double time){
        if(time != 0){
            return ((double)2.0 * Math.PI)  / time;
        } else {
            return 0;
        }
    }

    static double getVelocityFromAngularSpeed(double angularSpeed, double radius){
        //v = w * r
        return (angularSpeed * radius);
    }

    static double getAngularSpeedFromPeriod(double period){
        if(period != 0){
            return ((((double) 2.0) * Math.PI) / period);
        } else {
            return 0;
        }
    }

    static double getAngularSpeedFromFrequency(double frequency){
        return ((double)2.0 * Math.PI * frequency);
    }

    /*
        Frequency, F
            The number of cycles (full revolutions) per second
     */
    static double getFrequency(double angularSpeed){
        return (angularSpeed / (2.0 * Math.PI));
    }

    /*
        Get radians of a small movement about the center, dx
     */
    static double getRadsMovement(double dx, double radius){
        if(radius != 0){
            return (dx / radius);
        } else {
            return 0;
        }
    }
    static double getVelocityFromRadMovement(double dx, double dt){
        if(dt != 0){
            return dx / dt;
        } else {
            return 0;
        }
    }
}
