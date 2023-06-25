package Complex;

import NumberThoery.Numbers;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Ahmed
 */
/**
 * Represents a complex number with both real and imaginary parts.
 */
public class Complex {

    private final double real;
    private final double imaginary;
    private final double magnitude;
    private final double angle;

    /**
     * Gets the real part of the Complex number.
     *
     * @return The real part.
     */
    public double getReal() {
        return roundTo5Decimals(real);
    }

    /**
     * Gets the imaginary part of the Complex number.
     *
     * @return The imaginary part.
     */
    public double getImaginary() {
        return roundTo5Decimals(imaginary);
    }

    /**
     * Gets the magnitude (absolute value) of the Complex number.
     *
     * @return The magnitude.
     */
    private double getMagnitude() {
        return roundTo5Decimals(magnitude);
    }

    /**
     * Gets the phase (angle) of the Complex number.
     *
     * @return The angle in radians.
     */
    private double getAngle() {
        return roundTo5Decimals(angle);
    }

    /**
     * Private constructor to initialize the Complex object.
     *
     * @param real The real part of the complex number.
     * @param imaginary The imaginary part of the complex number.
     * @param magnitude The magnitude (absolute value) of the complex number.
     * @param angle The angle of the complex number.
     */
    private Complex(double real, double imaginary, double magnitude, double angle) {
        this.real = real;
        this.imaginary = imaginary;
        this.magnitude = magnitude;
        this.angle = angle;
    }

    /**
     * Creates a Complex number from Cartesian coordinates.
     *
     * @param real The real part of the complex number.
     * @param imaginary The imaginary part of the complex number.
     * @return A Complex object representing the Cartesian coordinates.
     */
    public static Complex createFromCartesian(double real, double imaginary) {
        double magnitude = calculateMagnitude(real, imaginary);
        double angle = calculateAngleInRadians(real, imaginary);
        return new Complex(real, imaginary, magnitude, angle);
    }

    /**
     * Creates a Complex number from polar coordinates.
     *
     * @param magnitude The magnitude (absolute value) of the complex number.
     * @param angle The angle of the complex number.
     * @return A Complex object representing the polar coordinates.
     */
    public static Complex createFromPolar(double magnitude, double angle) {
        double real = roundTo5Decimals(magnitude * Math.cos(angle));
        double imaginary = roundTo5Decimals(magnitude * Math.sin(angle));
        return new Complex(real, imaginary, magnitude, angle);
    }

    /**
     * Rounds a double value to five decimal places.
     *
     * @param value The double value to be rounded.
     * @return The rounded value with four decimal places.
     */
    private static double roundTo5Decimals(double value) {
        return Math.round(value * 100000) / 100000.0;
    }

    /**
     * Adds two Complex numbers and returns the sum.
     *
     * @param c1 The first Complex number.
     * @param c2 The second Complex number.
     * @return The sum of the two Complex numbers.
     */
    public static Complex sum(Complex c1, Complex c2) {
        double sumReal = roundTo5Decimals(c1.getReal() + c2.getReal());
        double sumImaginary = roundTo5Decimals(c1.getImaginary() + c2.getImaginary());
        return createFromCartesian(sumReal, sumImaginary);
    }

    /**
     * Subtracts the second Complex number from the first Complex number and
     * returns the result.
     *
     * @param c1 The first Complex number.
     * @param c2 The second Complex number.
     * @return The result of subtracting c2 from c1.
     */
    public static Complex subtract(Complex c1, Complex c2) {
        double resultReal = roundTo5Decimals(c1.getReal() - c2.getReal());
        double resultImaginary = roundTo5Decimals(c1.getImaginary() - c2.getImaginary());
        return createFromCartesian(resultReal, resultImaginary);
    }

    /**
     * Multiplies two Complex numbers and returns the product.
     *
     * @param c1 The first Complex number.
     * @param c2 The second Complex number.
     * @return The product of the two Complex numbers.
     */
    public static Complex multiply(Complex c1, Complex c2) {
        double productReal = roundTo5Decimals((c1.getReal() * c2.getReal()) - (c1.getImaginary() * c2.getImaginary()));
        double productImaginary = roundTo5Decimals((c1.getReal() * c2.getImaginary()) + (c2.getReal() * c1.getImaginary()));
        return createFromCartesian(productReal, productImaginary);
    }

    /**
     * Divides two Complex numbers and returns the result.
     *
     * @param dividend The dividend Complex number.
     * @param divisor The divisor Complex number.
     * @return The result of the division.
     */
    public static Complex divide(Complex dividend, Complex divisor) {
        double divisorSquared = divisor.getReal() * divisor.getReal() + divisor.getImaginary() * divisor.getImaginary();
        double realPart = ((dividend.getReal() * divisor.getReal()) + (dividend.getImaginary() * divisor.getImaginary())) / divisorSquared;
        double imaginaryPart = ((divisor.getReal() * dividend.getImaginary()) - (dividend.getReal() * divisor.getImaginary())) / divisorSquared;
        double roundedRealPart = roundTo5Decimals(realPart);
        double roundedImaginaryPart = roundTo5Decimals(imaginaryPart);
        return createFromCartesian(roundedRealPart, roundedImaginaryPart);
    }

    /**
     * Returns the complex number in rectangular form.
     *
     * @param c The Complex number.
     * @return The complex number in rectangular form as a string.
     */
    public static String toRectangularForm(Complex c) {
        double roundedReal = roundTo5Decimals(c.getReal());
        double roundedImaginary = roundTo5Decimals(c.getImaginary());
        return roundedReal + " + " + roundedImaginary + "i";
    }

    /**
     * Returns the polar form representation of the Complex number.
     *
     * @param complex The Complex number.
     * @return The polar form representation as a String.
     */
    public static String polarForm(Complex complex) {
        double magnitude = calculateMagnitude(complex);
        double angleInDegrees = calculateAngleInDegrees(complex);
        return magnitude + " * (cos(" + angleInDegrees + ") + i*sin(" + angleInDegrees + "))";
    }

    /**
     * Returns the polar form of a Complex number as a String with π
     * representation.
     *
     * @param complex The Complex number to convert to polar form.
     * @return The polar form of the Complex number as a String.
     */
    public static String getPolarFormWithPi(Complex complex) {
        double angleInDegrees = calculateAngleInDegrees(complex);
        int[] coeff = Numbers.convertDecimalToFraction(roundTo5Decimals(angleInDegrees / 180));
        String angleString = coeff[0] + "/" + coeff[1];

        String Polarform = calculateMagnitude(complex) + " * ( cos(PI" + angleString + ") + i sin(PI" + angleString + ") )";
        return Polarform;
    }

    /**
     * Creates a string representation of the complex number in angle notation.
     *
     * @param c The complex number.
     * @return The string representation in angle notation.
     */
    public static String angleNotation(Complex c) {
        String angleNotation = calculateMagnitude(c) + " " + "∠" + calculateAngleInDegrees(c);
        return angleNotation;
    }

    /**
     * Returns the complex number in exponential form.
     *
     * @param c The complex number.
     * @return The complex number in exponential form as a string.
     */
    public static String getExponentialForm(Complex c) {
        double magnitude = calculateMagnitude(c);
        double angleDegrees = calculateAngleInDegrees(c);
        return magnitude + " * e^(iπ" + angleDegrees + ")";
    }

    /**
     * Calculates the magnitude (absolute value) of a Complex number.
     *
     * @param complex The Complex number for which to calculate the magnitude.
     * @return The magnitude of the Complex number.
     */
    public static double calculateMagnitude(Complex complex) {
        return roundTo5Decimals(Math.sqrt((complex.getReal() * complex.getReal()) + (complex.getImaginary() * complex.getImaginary())));
    }

    /**
     * Calculates the magnitude (absolute value) of a Complex number given its
     * real and imaginary parts.
     *
     * @param real The real part of the Complex number.
     * @param imaginary The imaginary part of the Complex number.
     * @return The magnitude of the Complex number.
     */
    public static double calculateMagnitude(double real, double imaginary) {
        return roundTo5Decimals(Math.sqrt((real * real) + (imaginary * imaginary)));
    }

    /**
     * Calculates the angle in radians for a complex number.
     *
     * @param complex The complex number for which to calculate the angle.
     * @return The angle in radians.
     */
    public static double calculateAngleInRadians(Complex complex) {
        return roundTo5Decimals(Math.atan2(complex.getImaginary(), complex.getReal()));
    }

    /**
     * Calculates the angle (in radians) of a complex number given its real and
     * imaginary parts.
     *
     * @param real The real part of the complex number.
     * @param imaginary The imaginary part of the complex number.
     * @return The angle (in radians) of the complex number.
     */
    public static double calculateAngleInRadians(double real, double imaginary) {
        return roundTo5Decimals(Math.atan2(imaginary, real));
    }

    /**
     * Calculates the angle (in degrees) of a complex number given its real and
     * imaginary parts.
     *
     * @param real The real part of the complex number.
     * @param imaginary The imaginary part of the complex number.
     * @return The angle (in degrees) of the complex number.
     */
    public static double calculateAngleInDegrees(double real, double imaginary) {
        Complex complex = createFromCartesian(real, imaginary);
        double angleInRadians = complex.getAngle();
        double angleInDegrees = Math.toDegrees(angleInRadians);
        return roundTo5Decimals(angleInDegrees);
    }

    /**
     * Calculates the angle of a complex number in degrees.
     *
     * @param complex The complex number.
     * @return The angle of the complex number in degrees.
     */
    public static double calculateAngleInDegrees(Complex complex) {
        double angleInRadians = complex.getAngle();
        double angleInDegrees = Math.toDegrees(angleInRadians);
        return roundTo5Decimals(angleInDegrees);
    }

    /**
     * Checks if two complex numbers are equal.
     *
     * @param c1 The first complex number.
     * @param c2 The second complex number.
     * @return True if the complex numbers are equal, false otherwise.
     */
    public static boolean isEqual(Complex c1, Complex c2) {
        return Double.compare(c1.getReal(), c2.getReal()) == 0
                && Double.compare(c1.getImaginary(), c2.getImaginary()) == 0;
    }

    /**
     * Returns the conjugate of a given complex number.
     *
     * @param complex The complex number for which the conjugate is to be
     * calculated.
     * @return The conjugate of the complex number.
     */
    public static Complex calculateConjugate(Complex complex) {
        double realPart = roundTo5Decimals(complex.getReal());
        double imaginaryPart = roundTo5Decimals(-complex.getImaginary());
        return createFromCartesian(realPart, imaginaryPart);
    }

    /**
     * Converts polar coordinates to rectangular (Cartesian) coordinates.
     *
     * @param magnitude The magnitude (absolute value) of the complex number.
     * @param theta The angle (in degrees) of the complex number.
     * @return A Complex object representing the rectangular coordinates.
     */
    public static Complex toRectangular(double magnitude, double theta) {
        double real = roundTo5Decimals(magnitude * Math.cos(Math.toRadians(theta)));
        double imaginary = roundTo5Decimals(magnitude * Math.sin(Math.toRadians(theta)));
        return createFromCartesian(real, imaginary);
    }

    /**
     * Returns a string representation of the Complex number.
     *
     * @return The string representation of the Complex number.
     */
    @Override
    @SuppressWarnings("LocalVariableHidesMemberVariable")
    public String toString() {

        double real = getReal();
        double imaginary = getImaginary();

        if (imaginary < 0) {
            return real + " - " + (-imaginary) + "i";
        }
        if (imaginary == 0) {
            return String.valueOf(real);
        }
        if (real == 0) {
            return imaginary + "i";
        }

        return real + " + " + imaginary + "i";
    }

    /**
     * Calculates the complex power of a given complex number.
     *
     * @param c The complex number.
     * @param power The power to which the complex number is raised.
     * @return A list of complex numbers representing the result of the power
     * operation.
     */
    public static List<Complex> power(Complex c, double power) {
        List<Complex> result;
        int[] powerParts = Numbers.convertDecimalToFraction(power);
        double numerator = powerParts[0];
        int denominator = powerParts[1];

        if (powerParts[0] == 0 && powerParts[1] == 0) {
            // The power is irrational; the best we can do is get an approximation
            numerator = power;
            denominator = 1;
        }

        // Raise the number to the power of the numerator
        double absPow = Math.pow(c.getMagnitude(), numerator);
        double anglePow = c.getAngle() * numerator;
        double realPart = roundTo5Decimals(absPow * Math.cos(anglePow));
        double imgPart = roundTo5Decimals(absPow * Math.sin(anglePow));

        // Getting the denominator nth root
        result = nthRoot(createFromCartesian(realPart, imgPart), denominator);

        return result;
    }

    /**
     * Calculates the integer power of a given complex number.
     *
     * @param c The complex number.
     * @param power The power to which the complex number is raised.
     * @return A complex number representing the result of the power operation.
     */
    public static Complex power(Complex c, int power) {
        double absPow = Math.pow(c.getMagnitude(), power);
        double anglePow = c.getAngle() * power;
        double realPart = roundTo5Decimals(absPow * Math.cos(anglePow));
        double imgPart = roundTo5Decimals(absPow * Math.sin(anglePow));
        return createFromCartesian(realPart, imgPart);
    }

    /**
     * Calculates the nth roots of a given complex number.
     *
     * @param c The complex number.
     * @param n The degree of the root.
     * @return A list of complex numbers representing the nth roots.
     * @throws ArithmeticException If the 0th root is requested (undefined).
     */
    public static List<Complex> nthRoot(Complex c, int n) throws ArithmeticException {
        List<Complex> roots = new ArrayList<>();

        if (n == 1) {
            // No roots for n = 1
            roots.add(c);
        } else if (n <= 0) {
            throw new ArithmeticException("n must be greater than 0");
        } else {
            for (int k = 0; k < n; k++) {
                double angle = (c.getAngle() + (2 * Math.PI * k)) / n;
                if (angle > Math.PI) {
                    angle -= 2 * Math.PI;
                }
                double q = Math.pow(n, -1);
                Complex root = createFromPolar(Math.pow(c.getMagnitude(), q), angle);
                roots.add(root);
            }
        }
        return roots;
    }

    /**
     * Calculates the principal nth root of a given complex number.
     *
     * @param c The complex number.
     * @param n The degree of the root.
     * @return The principal nth root of the complex number.
     * @throws ArithmeticException If the 0th root is requested (undefined).
     */
    public static Complex principalRoot(Complex c, int n) throws ArithmeticException {
        List<Complex> roots = nthRoot(c, n);
        return roots.get(0);
    }
    
    
    /**
     * Calculates the complex power of a given complex number raised to another complex number.
     *
     * @param base   The base complex number.
     * @param exponent The exponent complex number.
     * @return The result of raising the base to the exponent.
     */
    public static Complex power(Complex base, Complex exponent) {
        double absBase = base.getMagnitude();
        double angleBase = base.getAngle();
        double realPart = exponent.getReal() * Math.log(absBase) - exponent.getImaginary() * angleBase;
        double imgPart = exponent.getReal() * angleBase + exponent.getImaginary() * Math.log(absBase);

        double resultReal = Math.exp(realPart) * Math.cos(imgPart);
        double resultImg = Math.exp(realPart) * Math.sin(imgPart);

        return createFromCartesian(resultReal, resultImg);
    }
    

    /**
     * Solves a quadratic equation and prints the results and can handle complex
     * roots.
     *
     * @param a The coefficient of x^2.
     * @param b The coefficient of x.
     * @param c The constant term.
     * @return An array of Complex objects representing the roots
     */
    public static Complex[] QuadraticEquationSolver(double a, double b, double c) {
        double discriminant = b * b - 4 * a * c;

        if (discriminant > 0) {
            double root1 = (-b + Math.sqrt(discriminant)) / (2 * a);
            double root2 = (-b - Math.sqrt(discriminant)) / (2 * a);
            return new Complex[]{createFromCartesian(root1, 0), createFromCartesian(root2, 0)};
        } else if (discriminant == 0) {
            double root = -b / (2 * a);
            return new Complex[]{createFromCartesian(root, 0)};
        } else {
            double realPart = -b / (2 * a);
            double imaginaryPart = Math.sqrt(Math.abs(discriminant)) / (2 * a);
            Complex root1 = createFromCartesian(realPart, imaginaryPart);
            Complex root2 = createFromCartesian(realPart, -imaginaryPart);
            return new Complex[]{root1, root2};
        }
    }

}
