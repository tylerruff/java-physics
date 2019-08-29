public class ClassicalMechanics {
    /*
        Newton's Laws of Motion
        1. An object will continue in it's state of rest,
        or uniform motion, unless acted on by an external force.
        2. F = m * a
     */
    static double getForce(double mass, double acceleration){
        //1 kg m/sec^2 = 1 Newton
        return (mass * acceleration);
    }
	static double getCombinedForces(double[] masses, double[] accelerations){
		if(masses.length == accelerations.length){
			double[] forces = new double[masses.length];
			double totalForce = 0;
			for(int i = 0; i < masses.length; i++){
				forces[i] = ClassicalMechanics.getForce(masses[i], accelerations[i]);
				totalForce += forces[i]; //add to total list
			}
			return totalForce;
		} else {
			return 0;
		}
	}
    /*
        3. To every action, there is an equal and opposite reaction.
     */
    static double velocity(double distance, double time){
        //velocity is a vector, having both magnitude and direction
        if(time != 0) {
            //Multiple return statements,
            //works in Java, does NOT work in PHP
            //(PHP insists you put the return statement at the bottom)
            return (distance / time);
        } else {
            return 0;
        }
    }

    /*
     Acceleration = the rate of change of velocity
     ex units: m/s^2
     a = (v-u) / t
     */
    static double getAcceleration(double vel1, double vel2, double time){
        if(time != 0){
            return ((vel1 - vel2) / time);
        } else {
            return 0;
        }
    }
    /*
        End Velocity after accelerating at a acceleration for t time
        u = start velocity
        v = end velocity
        //v = u + at
     */
    static double calculateTravel(double initVelocity, double acceleration, double time){
        return (initVelocity + (acceleration * time));
    }

    /*
        Total Distance traveled
        s = ut + 1/2 (v - u) * t
		at = (v - u)
		(acceleration multiplyied by time equals the total change in velocity)
     */
    static double calculateDistance(double startVelocity, double endVelocity, double time){
        return (startVelocity * time) + ((0.5) * (endVelocity - startVelocity) * time);
    }

    static double calculateDistanceAlt(double startVelocity, double acceleration, double time){
        //s = ut + 1/2 a * (t) ^ 2
        return ((startVelocity * time) + (0.5) * acceleration * Math.pow(time, (double)2));
    }

    //K = 1/2 m * v ^2
    static double kineticEnergy(double mass, double velocity){
        return ((0.5) * mass * Math.pow(velocity, (double)2.0));
    }

    //F * s = 1/2(m * v ^ 2) - 1/2 (m * u ^ 2)
    //Total work done, the difference in kinetic energy
    //force times total distance traveled
    //expressed in Newton Meters or Joules
    static double workDone(double mass, double vel1, double vel2){
        return (ClassicalMechanics.kineticEnergy(mass, vel2) - ClassicalMechanics.kineticEnergy(mass, vel1));
    }

    /*
        Momentum in 1D
        p = m * v
        Expressed in kg  m/sec
        (Force is the rate of change of momentum)
     */
    static double momentum(double mass, double velocity){
        return (mass * velocity);
    }
    //Calculate force by the rate of change of start & end momentum
    static double calculateForcefromMomentum(double initMo, double finalMo, double time){
        if(time != 0){
            return ((finalMo - initMo) / time);
        } else {
            return 0;
        }
    }
    static double calculateForcefromVelocity(double initVel, double finalVel, double mass, double time) {
        if (time != 0) {
            return (mass * (finalVel - initVel)) / time;
        } else {
            return 0;
        }
    }

    /*
        Elastic and Inelastic Collisions
        Elastic = energy is conserved in the form of Kinetic Energy
        -------------------------------------------------------------
        Inelastic = energy is conserved, but not necessarily as kinetic energy,
        ex: energy lost as heat, sound, or light.
        -------------------------------------------------------------
     */

    /*
     When an object has potential to fall (or be thrown/or otherwise accelerated),
     it has potential energy.
     Potential Energy = m * g * h
     m = mass
     g = gravitational constant
     h = height from ground (ground is a position that would cause the body to rest, were it to fall)
     */
    static double potentialEnergy(double mass, double gracvConstant, double height){
        return mass * gracvConstant * height;
    }
    /*
        Once the fall, throw, or acceleration has been executed,
        that potential energy has been (and will be) converting that potential energy
        into kintetic energy,
        k = (1 / 2) m * (v ^ 2)
        until that body once again reaches a rest state.
        --------------------------------------------------
     */
    /*
    Impulse
        The change of momentum
        I = mv - mu
        I = F * t
        (force being the RATE of change of momentum)
        units: Newton Seconds
     */
    static double calculateImpulse(double mass, double velInit, double velEnd){
        //If an object, traveling at velocity v1,
        //hits a brick wall, then the impulse will show as follows:
        //I = 0 - mu, as the final velocity is 0 (assuming the wall stops the object)
        return ((velEnd * mass) - (velInit * mass));
    }

    /*
        Power = Energy or work per second
        Joules / Second, or Watts
     */
    static double getPower(double work, double time){
        if (time != 0) {
            return work / time;
        } else {
            return 0;
        }
    }
    /*
        Example of method overloading...
        Works in Java, does NOT work in PHP...
     */
    static double getPower(double force, double distance, double time){
        if(time != 0){
            return (force * distance) / time;
        } else {
            return 0;
        }
    }

    /*
        Moments
        M = Moment is the turning effect of a force
     */
    static double calculateMoment(double force, double distance){
        //distance = distance from pivot point
        return force * distance;
    }
    /*
        Seesaw example:
            m1               m2
            o                o
            |<--b--><---a--->|
        -------------------------
                   /\
            V                V
            F1               F2
        F1 = m1 * g
        F2 = m2 * g
            where,
            m1 = mass of person 1, F1 = force exerted by person 1
            m2 = mass of person 2, F2 = force exerted by person 2
            g = gravitational constant
            b = person 1 distance from pivot point
            a = person 2 distance from pivot point
       If, m1 * g * b = m2 * g * a
       Then, the seesaw will NOT move.
     */

    /*
    Torques
        t = F * d
       If a force is applied to either extant side to a pivot point,
       the torque can be found by multiplying the force by the distance
       of a single side from the pivot point.
       Pinwheel example
       A force is applied to either side of this pinwheel
              ************
            ****************
          ********************
     (F)<-----d---->*<----d----->(F)
          ********************
            ****************
             *************
               t = F * d
     */
    /*
        Projectile Motion
            When a projectile is fired at a certain velocity,
            and certain angle (theta), by which the projectile covers a certain trajectory.
            This trajectory often resembles an "arc"
     */
    //to calculate Range & max height
    //angle given in degrees
    static double projectileVerticalVelocity(double velocity, double angle){
        return velocity * Math.sin(ClassicalMechanics.degreeToRad(angle));
    }
    static double projectileVerticalTotalFallTime(double initVelocity, double gravity){
        if(gravity > 0){ gravity *= (-1); }
        //v = u + at
        //v = final velocity
        //u = initial velocity
        //a = acceleration (+ upwards, - downwards)
        //t = time
        //(v - u) / a =  t
        //if final velocity is zero
        return (((double)0) - initVelocity) / gravity;
    }
    //t = time it takes to get to position, P
    static double projectileHeightAtSpecificPoint(double initalVerticalVelocity, double gravity, double time){
        if (gravity > 0){
            gravity *= (-1);
        }
        return (initalVerticalVelocity * time) + (0.5) * gravity * Math.pow(time, 2);
    }
    //angle given in degrees
    static double projectileHorizontalVelocity(double velocity, double angle){
        return velocity * Math.cos(ClassicalMechanics.degreeToRad(angle));
    }
    //get vertical distance from origin at time, t
    static double projectileHorizontalDistance(double horizontalVelocity, double time){
        return horizontalVelocity * time;
    }
    //max height
    //v^2 = u^2 + 2 a s
    static double getMaxHeight(double initVelocity, double gravity){
        //0 = u ^ 2 + 2 * a * s
        //-u^2 / 2a = s
        if(gravity > 0){ gravity *= (-1); }
        return ((-1) * Math.pow(initVelocity, 2)) / (2.0 * gravity);
    }

    public static double degreeToRad(double deg){
        return deg * (Math.PI / 180.0);
    }
	
	/*
		Hooke's Law
			Getting the responsive force of a spring, given a specific spring constant.
	*/
	//k = F
	//    _
	//	  x
	public static double getSpringConstant(double forceReturned, double displacement){
		if(displacement != 0){
			return (-1) * (forceReturned / displacement);
		} else {
			return 0;
		}
	}
	public static double calculateElasticForce(double springConstant, double changeInX){
		//F = − kΔx
		return ((-1) * springConstant * changeInX);
	}
    public static void main(String[] args){
        //Example 1:
        // Projectile with inital vertical velocity 10 m/s^2
        // (The vertical velocity is calculated by
        //  multiplying the total init velocity (20 m/s^2) by sin(theta),
        //  where theta is the angle of the projectile's launch)
        // Gravitational constant simplified to 10
        // The maximum vertical height will be 5 meters
        System.out.println(ClassicalMechanics.getMaxHeight(ClassicalMechanics.projectileVerticalVelocity(20, 30), 10));
        //Output: 4.999999999999998 ~ 5.0
    }

}
