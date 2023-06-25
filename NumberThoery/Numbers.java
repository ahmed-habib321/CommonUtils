/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package NumberThoery;

import java.text.DecimalFormat;

/**
 *
 * @author Ahmed
 */
public class Numbers {

    /**
     * Checks whether a given number is rational.
     *
     * @param number The input number to be evaluated. It can be a positive or
     * negative double value.
     * @return {@code true} if the number is determined to be rational,
     * {@code false} otherwise. The method returns {@code true} if the number is
     * close enough to an integer based on the calculations. The threshold for
     * determining closeness to an integer is set to 0.001 (1e-3). The method
     * performs a maximum of 20 iterations to make the determination.
     */
    public static boolean isRational(double number) {
        number = Math.abs(number);
        final double threshold = 1e-3;
        final int maxIterations = 20;

        for (int i = 0; i < maxIterations; i++) {
            double floorValue = Math.floor(number);
            if (number - floorValue < threshold) {
                return true;
            }
            number = 1 / (number - floorValue);
        }

        return false;
    }

    /**
     * Calculates the greatest common divisor (GCD) of two integers.
     *
     * @param a the first integer
     * @param b the second integer
     * @return the GCD of the two integers
     */
    public static int GCD(int a, int b) {
        a = Math.abs(a);
        b = Math.abs(b);
        int gcd = 1;
        for (int i = 1; i <= a && i <= b; i++) {
            if (a % i == 0 && b % i == 0) {
                gcd = i;
            }
        }
        return gcd;
    }

    /**
     * Calculates the factor of 10 needed to multiply a decimal number to make
     * it an integer.
     *
     * @param decimalNumber The decimal number to be multiplied.
     * @return The factor of 10 required to make the decimal number an integer.
     */
    public static int getDecimalToIntegerFactor(double decimalNumber) {
        int factor = 1;
        boolean isIntegerReached = isInteger(decimalNumber, factor);

        while (!isIntegerReached) {
            factor *= 10;
            isIntegerReached = isInteger(decimalNumber, factor);
        }

        return factor;
    }

    /**
     * Checks if the decimal number multiplied by the factor equals an integer.
     *
     * @param decimalNumber The decimal number to be multiplied.
     * @param factor The factor of 10 to be applied.
     * @return True if the multiplication results in an integer, false
     * otherwise.
     */
    private static boolean isInteger(double decimalNumber, int factor) {
        return (int) (decimalNumber * factor) == decimalNumber * factor;
    }

    /**
     * Converts a decimal number into a simplified fraction representation. The
     * fraction accuracy is up to 9 digits after the decimal point. If the
     * decimal number has more than 9 digits after the decimal point, the method
     * returns "0/0".
     *
     * @param decimalNumber The decimal number to be converted.
     * @return The simplified fraction representation of the decimal number.
     */
    public static int[] convertDecimalToFraction(double decimalNumber) {
        int denominator = getDecimalToIntegerFactor(decimalNumber);
        int numerator = (int) (decimalNumber * denominator);
        int gcd = GCD(numerator, denominator);
        
        int[] fraction = getNumeratorAndDenominator((numerator/gcd) + "/" + (denominator/gcd));
        return fraction;
    }

    /**
     * Parses a string representation of a fraction and returns the numerator
     * and denominator as an array.
     *
     * @param fractionString The string representation of the fraction in the
     * format "numerator/denominator".
     * @return An array containing the numerator and denominator.
     * @throws NumberFormatException If the input string is not in the correct
     * format or cannot be parsed.
     */
    public static int[] getNumeratorAndDenominator(String fractionString) throws NumberFormatException {
        String[] parts = fractionString.split("/");
        if (parts.length != 2) {
            throw new NumberFormatException("Invalid fraction format: " + fractionString);
        }
        int numerator = Integer.parseInt(parts[0]);
        int denominator = Integer.parseInt(parts[1]);
        return new int[]{numerator, denominator};
    }

    /**
     * Converts a decimal number to the specified base representation.
     *
     * @param decimalNumber The decimal number to be converted.
     * @param base The base to which the number should be converted (between 2
     * and 36).
     * @return The string representation of the decimal number in the specified
     * base.
     * @throws IllegalArgumentException If the specified base is not within the
     * valid range.
     */
    public static String convertToBase(int decimalNumber, int base) throws IllegalArgumentException {
        if (base < 2 || base > 36) {
            throw new IllegalArgumentException("Invalid base: " + base + ". Base must be between 2 and 36.");
        }

        StringBuilder result = new StringBuilder();
        int number = decimalNumber;

        while (number > 0) {
            int remainder = number % base;
            char digit = (remainder < 10) ? (char) (remainder + '0') : (char) (remainder - 10 + 'A');
            result.insert(0, digit);
            number /= base;
        }

        return result.toString();
    }

    /**
     * Converts a number from the specified base to its decimal representation.
     *
     * @param number The string representation of the number in the specified
     * base.
     * @param base The base of the number (between 2 and 36).
     * @return The decimal representation of the number.
     * @throws IllegalArgumentException If the specified base is not within the
     * valid range or the number contains invalid digits for the base.
     */
    public static int convertToDecimal(String number, int base) throws IllegalArgumentException {
        if (base < 2 || base > 36) {
            throw new IllegalArgumentException("Invalid base: " + base + ". Base must be between 2 and 36.");
        }

        int result = 0;

        for (int i = 0; i < number.length(); i++) {
            char digit = number.charAt(i);
            int value;

            if (Character.isDigit(digit)) {
                value = digit - '0';
            } else if (Character.isLetter(digit)) {
                value = Character.toUpperCase(digit) - 'A' + 10;
            } else {
                throw new IllegalArgumentException("Invalid digit in number: " + digit);
            }

            if (value >= base) {
                throw new IllegalArgumentException("Invalid digit in number: " + digit);
            }

            result = result * base + value;
        }

        return result;
    }

    /**
     * Checks if two integers are coprime, meaning they have no common positive
     * integer divisors other than 1.
     *
     * @param a The first integer.
     * @param b The second integer.
     * @return True if the integers are coprime, false otherwise.
     */
    public static boolean isCoprime(int a, int b) {
        return GCD(a, b) == 1;
    }

    /**
     * Calculates the Least Common Multiple (LCM) of two integers.
     *
     * @param a The first integer.
     * @param b The second integer.
     * @return The LCM of the two integers.
     */
    public static int calculateLCM(int a, int b) {
        return (a * b) / GCD(a, b);
    }

    /**
     * Calculates the factorial of a given number.
     *
     * @param number The number for which to calculate the factorial.
     * @return The factorial of the given number.
     * @throws IllegalArgumentException If the number is negative of decimal.
     */
    public static int factorial(int number) {
        if (number < 0 || !isInteger(number, 1)) {
            throw new IllegalArgumentException("Number cannot be negative or decimal.");
        }
        int factorial = 1;
        while (number > 0) {
            factorial *= number;
            number--;
        }
        return factorial;
    }

    /**
     * Solves a linear Diophantine equation of the form 'ax + by = c', where
     * 'a', 'b', and 'c' are integers.
     *
     * @param a The coefficient of 'x'.
     * @param b The coefficient of 'y'.
     * @param c The constant term.
     * @return An array containing the solution [x, y] if it exists; [0, 0]
     * otherwise.
     */
    public static int[] solveLinearDiophantine(int a, int b, int c) {
        int[] result = new int[2];

        // Find the greatest common divisor (GCD) of 'a' and 'b'
        int gcd = GCD(a, b);

        // Check if the equation has a solution
        if (c % gcd != 0) {
            result[0] = 0;  // No solution
            result[1] = 0;
            return result;
        }

        // Apply the extended Euclidean algorithm to find Bezout's identity
        int[] bezoutCoefficients = calculateBezoutCoefficients(a, b);
        int x0 = bezoutCoefficients[0] * (c / gcd);
        int y0 = bezoutCoefficients[1] * (c / gcd);

        result[0] = x0;
        result[1] = y0;
        return result;
    }

    /**
     * Calculates the coefficients of Bezout's identity for two numbers using
     * the extended Euclidean algorithm.
     *
     * @param a The first number.
     * @param b The second number.
     * @return An array containing the coefficients [x, y] of Bezout's identity
     * for 'a' and 'b'.
     */
    private static int[] calculateBezoutCoefficients(int a, int b) {
        if (b == 0) {
            return new int[]{1, 0};
        }
        int[] coefficients = calculateBezoutCoefficients(b, a % b);
        int x = coefficients[1];
        int y = coefficients[0] - (a / b) * coefficients[1];
        return new int[]{x, y};
    }

    /**
     * Checks if a number is a perfect number. A perfect number is a positive
     * integer that is equal to the sum of its proper divisors (excluding
     * itself).
     *
     * @param number The number to be checked.
     * @return True if the number is a perfect number, false otherwise.
     */
    public static boolean isPerfectNumber(int number) {
        if (number <= 1) {
            return false;
        }

        int sumOfDivisors = 1;  // Initialize the sum with 1 (as 1 is always a proper divisor)

        for (int divisor = 2; divisor <= Math.sqrt(number); divisor++) {
            if (number % divisor == 0) {
                sumOfDivisors += divisor;
                int otherDivisor = number / divisor;
                if (otherDivisor != divisor) {
                    sumOfDivisors += otherDivisor;
                }
            }
        }
        return sumOfDivisors == number;
    }

    /**
     * Checks if a given number is prime.
     *
     * @param number The number to be checked.
     * @return {@code true} if the number is prime, {@code false} otherwise.
     */
    public static boolean isPrime(int number) {
        if (number <= 1) {
            return false;
        }
        for (int i = 2; i <= Math.sqrt(number); i++) {
            if (number % i == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Rationalizes a decimal number by removing or truncating trailing zeros
     * after a specified number of decimal places.
     *
     * @param number The decimal number to be rationalized.
     * @param decimalPlaces The number of decimal places after which trailing
     * zeros should be removed.
     * @return The rationalized decimal number.
     */
    public static double rationalizeDecimal(double number, int decimalPlaces) {
        double factor = Math.pow(10, decimalPlaces);
        int multipliedNumber = (int) (number * factor);
        double rationalizedNumber = multipliedNumber / factor;

        return rationalizedNumber;
    }


/**
 * Formats a double number as a string using a decimal format.
 *
 * @param number the double number to be formatted
 * @return the formatted string representation of the number
 */
public static String formatDouble(double number) {
    DecimalFormat decimalFormat = new DecimalFormat("#.###################");
    decimalFormat.setDecimalSeparatorAlwaysShown(false);
    return decimalFormat.format(number);
}

}
