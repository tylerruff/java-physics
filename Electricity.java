public class Electricity {
    /*
    Two charges, Q1 & Q2
    they are a distance r apart
     */
    static final double permittivityConstant = (8.85 * Math.pow(10, -12));
    static double calculateForce(double q1, double q2, double r){
        return (q1 * q2) / (4.0 * Math.PI * permittivityConstant * Math.pow(r, 2));
    }
}
