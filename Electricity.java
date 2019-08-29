public class Electricity {
    /*
    Two charges, Q1 & Q2
    they are a distance r apart
     */
    static final double permittivityConstant = (8.85 * Math.pow(10, -12));
    //P = V * I
	static double calculatePower(double voltage, double current){
		return (voltage * current);
	}
	//P = R * I^2
	static double calculatePowerAlt(double resistance, double current){
		return resistance * Math.pow(current, 2);
	}
	//E = V * I * t
	static double calculateEnergy(double voltage, double current, double time){
		return (voltage * current * time);
	}
	//E = V * Q
	static double calculateEnergy(double voltage, double charge){
		return (voltage * charge);
	}
	//W = V * Q
	static double calculateWork(double voltage, double charge){
		return (voltage * charge);
	}
	//P = E / t
	static double getPowerFromEnergy(double energy, double time){
		if(time != 0){
			return (energy / time);
		} else {
			return 0;
		}
	}
	static double getCharge(double current, double time){
		return (current * time);
	}
	//V = I * R
	static double calculateVoltage(double current, double resistance){
		return (current * resistance);
	}
	
	static double calculateForce(double q1, double q2, double r){
        return (q1 * q2) / (4.0 * Math.PI * permittivityConstant * Math.pow(r, 2));
    }
	
	/**
		Ohm's Law
			the current through a conductor between two points is directly proportional to the potential difference or voltage across the two points, and inversely proportional to the resistance between them
	*/
	//Current given in Amperes
	//Electrical Potential given in Volts
	//Resistance given in Ohms
	static double calculateCurrent(double potential, double resistance){
		if(resistance != 0){
			return (potential / resistance);
		} else {
			return 0;
		}
	}
	static double calculatePotential(double resistance, double current){
		return (resistance * current);
	}
	static double calculateResistance(double potential, double current){
		if(current != 0){
			return (potential / current);
		} else {
			return 0;
		}
	}
	
}
