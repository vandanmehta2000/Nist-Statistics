public class Statistics {

    private static final double DEFAULT_EPSILON = 1.0E-14;

    private static final double[] DELTA = {
            0.833333333333333333333333333333E-01,
            -0.277777777777777777777777752282E-04,
            0.793650793650793650791732130419E-07,
            -0.595238095238095232389839236182E-09,
            0.841750841750832853294451671990E-11,
            -0.191752691751854612334149171243E-12,
            0.641025640510325475730918472625E-14,
            -0.295506514125338232839867823991E-15,
            0.179643716359402238723287696452E-16,
            -0.139228964661627791231203060395E-17,
            0.133802855014020915603275339093E-18,
            -0.154246009867966094273710216533E-19,
            0.197701992980957427278370133333E-20,
            -0.234065664793997056856992426667E-21,
            0.171348014966398575409015466667E-22
    };

    public static final double LANCZOS_G = 607.0 / 128.0;

    private static final double[] LANCZOS = {
            0.99999999999999709182,
            57.156235665862923517,
            -59.597960355475491248,
            14.136097974741747174,
            -0.49191381609762019978,
            0.33994649984811888699e-4,
            0.46523628927048575665e-4,
            -0.98374475304879564677e-4,
            0.15808870322491248884e-3,
            -0.21026444172410488319e-3,
            0.21743961811521264320e-3,
            -0.16431810653676389022e-3,
            0.84418223983852743293e-4,
            -0.26190838401581408670e-4,
            0.36899182659531622704e-5,
            };

    private static final double HALF_LOG_TWO_PI = 0.9189385332046727;

    private static final double HALF_LOG_2_PI = 0.5 * Math.log(2.0 * Math.PI);

    private static final double SQRT_TWO_PI = 2.506628274631000502;

    private static final double INV_GAMMA1P_M1_A0 = 0.611609510448141581788E-08;
    private static final double INV_GAMMA1P_M1_A1 = 0.624730830116465516210E-08;
    private static final double INV_GAMMA1P_M1_B1 = 0.203610414066806987300E+00;
    private static final double INV_GAMMA1P_M1_B2 = 0.266205348428949217746E-01;
    private static final double INV_GAMMA1P_M1_B3 = 0.493944979382446875238E-03;
    private static final double INV_GAMMA1P_M1_B4 = -0.851419432440314906588E-05;
    private static final double INV_GAMMA1P_M1_B5 = -0.643045481779353022248E-05;
    private static final double INV_GAMMA1P_M1_B6 = 0.992641840672773722196E-06;
    private static final double INV_GAMMA1P_M1_B7 = -0.607761895722825260739E-07;
    private static final double INV_GAMMA1P_M1_B8 = 0.195755836614639731882E-09;
    private static final double INV_GAMMA1P_M1_P0 = 0.6116095104481415817861E-08;
    private static final double INV_GAMMA1P_M1_P1 = 0.6871674113067198736152E-08;
    private static final double INV_GAMMA1P_M1_P2 = 0.6820161668496170657918E-09;
    private static final double INV_GAMMA1P_M1_P3 = 0.4686843322948848031080E-10;
    private static final double INV_GAMMA1P_M1_P4 = 0.1572833027710446286995E-11;
    private static final double INV_GAMMA1P_M1_P5 = -0.1249441572276366213222E-12;
    private static final double INV_GAMMA1P_M1_P6 = 0.4343529937408594255178E-14;
    private static final double INV_GAMMA1P_M1_Q1 = 0.3056961078365221025009E+00;
    private static final double INV_GAMMA1P_M1_Q2 = 0.5464213086042296536016E-01;
    private static final double INV_GAMMA1P_M1_Q3 = 0.4956830093825887312020E-02;
    private static final double INV_GAMMA1P_M1_Q4 = 0.2692369466186361192876E-03;
    private static final double INV_GAMMA1P_M1_C = -0.422784335098467139393487909917598E+00;
    private static final double INV_GAMMA1P_M1_C0 = 0.577215664901532860606512090082402E+00;
    private static final double INV_GAMMA1P_M1_C1 = -0.655878071520253881077019515145390E+00;
    private static final double INV_GAMMA1P_M1_C2 = -0.420026350340952355290039348754298E-01;
    private static final double INV_GAMMA1P_M1_C3 = 0.166538611382291489501700795102105E+00;
    private static final double INV_GAMMA1P_M1_C4 = -0.421977345555443367482083012891874E-01;
    private static final double INV_GAMMA1P_M1_C5 = -0.962197152787697356211492167234820E-02;
    private static final double INV_GAMMA1P_M1_C6 = 0.721894324666309954239501034044657E-02;
    private static final double INV_GAMMA1P_M1_C7 = -0.116516759185906511211397108401839E-02;
    private static final double INV_GAMMA1P_M1_C8 = -0.215241674114950972815729963053648E-03;
    private static final double INV_GAMMA1P_M1_C9 = 0.128050282388116186153198626328164E-03;
    private static final double INV_GAMMA1P_M1_C10 = -0.201348547807882386556893914210218E-04;
    private static final double INV_GAMMA1P_M1_C11 = -0.125049348214267065734535947383309E-05;
    private static final double INV_GAMMA1P_M1_C12 = 0.113302723198169588237412962033074E-05;
    private static final double INV_GAMMA1P_M1_C13 = -0.205633841697760710345015413002057E-06;

    public static double regularizedBeta(double x, double a, double b) {
        return regularizedBeta(x, a, b, DEFAULT_EPSILON, Integer.MAX_VALUE);
    }

    public static double regularizedBeta(double x,
                                         final double a, final double b,
                                         double epsilon, int maxIterations) {
        double ret;

        if (Double.isNaN(x) ||
                Double.isNaN(a) ||
                Double.isNaN(b) ||
                x < 0 ||
                x > 1 ||
                a <= 0 ||
                b <= 0) {
            ret = Double.NaN;
        }
        else if (x > (a + 1) / (2 + b + a) &&
                1 - x <= (b + 1) / (2 + b + a)) {
            ret = 1 - regularizedBeta(1 - x, b, a, epsilon, maxIterations);
        }
        else {
            ret = Math.exp((a * Math.log(x)) + (b * Math.log1p(-x)) -
                                   Math.log(a) - logBeta(a, b)) *
                    1.0 / evaluate(x, epsilon, maxIterations, a, b);
        }

        return ret;
    }

    public static double getA(int n, double x) {
        return 1.0;
    }

    public static double getB(int n, double x, double a, double b) {
        double ret;
        double m;
        if (n % 2 == 0) { // even
            m = n / 2.0;
            ret = (m * (b - m) * x) /
                    ((a + (2 * m) - 1) * (a + (2 * m)));
        }
        else {
            m = (n - 1.0) / 2.0;
            ret = -((a + m) * (a + b + m) * x) /
                    ((a + (2 * m)) * (a + (2 * m) + 1.0));
        }
        return ret;
    }

    public static double evaluate(double x, double epsilon, int maxIterations, double aold,
                                  double bold) {
        final double small = 1.0e-50;
        double hPrev = getA(0, x);

        // use the value of small as epsilon criteria for zero checks
        if (equality(hPrev, 0.0, small)) {
            hPrev = small;
        }

        int n = 1;
        double dPrev = 0.0;
        double cPrev = hPrev;
        double hN = hPrev;

        while (n < maxIterations) {
            final double a = getA(n, x);
            final double b = getB(n, x, aold, bold);

            double dN = a + b * dPrev;
            if (equality(dN, 0.0, small)) {
                dN = small;
            }
            double cN = a + b / cPrev;
            if (equality(cN, 0.0, small)) {
                cN = small;
            }

            dN = 1 / dN;
            final double deltaN = cN * dN;
            hN = hPrev * deltaN;

            if (Math.abs(deltaN - 1.0) < epsilon) {
                break;
            }

            dPrev = dN;
            cPrev = cN;
            hPrev = hN;
            n++;
        }

        return hN;
    }

    public static boolean equality(double x, double y, double eps) {
        return Math.abs(y - x) <= eps;
    }

    public static double logBeta(final double p, final double q) {
        if (Double.isNaN(p) || Double.isNaN(q) || (p <= 0.0) || (q <= 0.0)) {
            return Double.NaN;
        }

        final double a = Math.min(p, q);
        final double b = Math.max(p, q);
        if (a >= 10.0) {
            final double w = sumDeltaMinusDeltaSum(a, b);
            final double h = a / b;
            final double c = h / (1.0 + h);
            final double u = -(a - 0.5) * Math.log(c);
            final double v = b * Math.log1p(h);
            if (u <= v) {
                return (((-0.5 * Math.log(b) + HALF_LOG_TWO_PI) + w) - u) - v;

            }
            else {
                return (((-0.5 * Math.log(b) + HALF_LOG_TWO_PI) + w) - v) - u;

            }

        }
        else if (a > 2.0) {
            if (b > 1000.0) {
                final int n = (int) Math.floor(a - 1.0);
                double prod = 1.0;
                double ared = a;
                for (int i = 0; i < n; i++) {
                    ared -= 1.0;
                    prod *= ared / (1.0 + ared / b);

                }
                return (Math.log(prod) - n * Math.log(b)) +
                        (logGamma(ared) +
                                logGammaMinusLogGammaSum(ared, b));

            }
            else {
                double prod1 = 1.0;
                double ared = a;
                while (ared > 2.0) {
                    ared -= 1.0;
                    final double h = ared / b;
                    prod1 *= h / (1.0 + h);

                }
                if (b < 10.0) {
                    double prod2 = 1.0;
                    double bred = b;
                    while (bred > 2.0) {
                        bred -= 1.0;
                        prod2 *= bred / (ared + bred);

                    }
                    return Math.log(prod1) +
                            Math.log(prod2) +
                            (logGamma(ared) +
                                    (logGamma(bred) -
                                            logGammaSum(ared, bred)));

                }
                else {
                    return Math.log(prod1) +
                            logGamma(ared) +
                            logGammaMinusLogGammaSum(ared, b);

                }

            }

        }
        else if (a >= 1.0) {
            if (b > 2.0) {
                if (b < 10.0) {
                    double prod = 1.0;
                    double bred = b;
                    while (bred > 2.0) {
                        bred -= 1.0;
                        prod *= bred / (a + bred);

                    }
                    return Math.log(prod) +
                            (logGamma(a) +
                                    (logGamma(bred) -
                                            logGammaSum(a, bred)));

                }
                else {
                    return logGamma(a) +
                            logGammaMinusLogGammaSum(a, b);

                }

            }
            else {
                return logGamma(a) +
                        logGamma(b) -
                        logGammaSum(a, b);
            }

        }
        else {
            if (b >= 10.0) {
                return logGamma(a) +
                        logGammaMinusLogGammaSum(a, b);

            }
            else {
                return Math.log(gamma(a) * gamma(b) /
                                        gamma(a + b));
            }
        }
    }

    private static double sumDeltaMinusDeltaSum(final double p,
                                                final double q) {
        final double a = Math.min(p, q);
        final double b = Math.max(p, q);
        final double sqrtT = 10.0 / a;
        final double t = sqrtT * sqrtT;
        double z = DELTA[DELTA.length - 1];
        for (int i = DELTA.length - 2; i >= 0; i--) {
            z = t * z + DELTA[i];
        }
        return z / a + deltaMinusDeltaSum(a, b);
    }

    public static double logGamma(double x) {
        double ret;

        if (Double.isNaN(x) || (x <= 0.0)) {
            ret = Double.NaN;
        }
        else if (x < 0.5) {
            return logGamma1p(x) - Math.log(x);
        }
        else if (x <= 2.5) {
            return logGamma1p((x - 0.5) - 0.5);
        }
        else if (x <= 8.0) {
            final int n = (int) Math.floor(x - 1.5);
            double prod = 1.0;
            for (int i = 1; i <= n; i++) {
                prod *= x - i;
            }
            return logGamma1p(x - (n + 1)) + Math.log(prod);
        }
        else {
            double sum = lanczos(x);
            double tmp = x + LANCZOS_G + 0.5;
            ret = ((x + 0.5) * Math.log(tmp)) - tmp +
                    HALF_LOG_2_PI + Math.log(sum / x);
        }

        return ret;
    }

    public static double lanczos(final double x) {
        double sum = 0.0;
        for (int i = LANCZOS.length - 1; i > 0; --i) {
            sum += LANCZOS[i] / (x + i);
        }
        return sum + LANCZOS[0];
    }

    private static double logGammaSum(final double a, final double b) {

        final double x = (a - 1.0) + (b - 1.0);
        if (x <= 0.5) {
            return logGamma1p(1.0 + x);
        }
        else if (x <= 1.5) {
            return logGamma1p(x) + Math.log1p(x);
        }
        else {
            return logGamma1p(x - 1.0) + Math.log(x * (1.0 + x));
        }
    }

    public static double logGamma1p(final double x) {
        return -Math.log1p(invGamma1pm1(x));
    }

    public static double invGamma1pm1(final double x) {
        final double ret;
        final double t = x <= 0.5 ? x : (x - 0.5) - 0.5;
        if (t < 0.0) {
            final double a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
            double b = INV_GAMMA1P_M1_B8;
            b = INV_GAMMA1P_M1_B7 + t * b;
            b = INV_GAMMA1P_M1_B6 + t * b;
            b = INV_GAMMA1P_M1_B5 + t * b;
            b = INV_GAMMA1P_M1_B4 + t * b;
            b = INV_GAMMA1P_M1_B3 + t * b;
            b = INV_GAMMA1P_M1_B2 + t * b;
            b = INV_GAMMA1P_M1_B1 + t * b;
            b = 1.0 + t * b;

            double c = INV_GAMMA1P_M1_C13 + t * (a / b);
            c = INV_GAMMA1P_M1_C12 + t * c;
            c = INV_GAMMA1P_M1_C11 + t * c;
            c = INV_GAMMA1P_M1_C10 + t * c;
            c = INV_GAMMA1P_M1_C9 + t * c;
            c = INV_GAMMA1P_M1_C8 + t * c;
            c = INV_GAMMA1P_M1_C7 + t * c;
            c = INV_GAMMA1P_M1_C6 + t * c;
            c = INV_GAMMA1P_M1_C5 + t * c;
            c = INV_GAMMA1P_M1_C4 + t * c;
            c = INV_GAMMA1P_M1_C3 + t * c;
            c = INV_GAMMA1P_M1_C2 + t * c;
            c = INV_GAMMA1P_M1_C1 + t * c;
            c = INV_GAMMA1P_M1_C + t * c;
            if (x > 0.5) {
                ret = t * c / x;
            }
            else {
                ret = x * ((c + 0.5) + 0.5);
            }
        }
        else {
            double p = INV_GAMMA1P_M1_P6;
            p = INV_GAMMA1P_M1_P5 + t * p;
            p = INV_GAMMA1P_M1_P4 + t * p;
            p = INV_GAMMA1P_M1_P3 + t * p;
            p = INV_GAMMA1P_M1_P2 + t * p;
            p = INV_GAMMA1P_M1_P1 + t * p;
            p = INV_GAMMA1P_M1_P0 + t * p;

            double q = INV_GAMMA1P_M1_Q4;
            q = INV_GAMMA1P_M1_Q3 + t * q;
            q = INV_GAMMA1P_M1_Q2 + t * q;
            q = INV_GAMMA1P_M1_Q1 + t * q;
            q = 1.0 + t * q;

            double c = INV_GAMMA1P_M1_C13 + (p / q) * t;
            c = INV_GAMMA1P_M1_C12 + t * c;
            c = INV_GAMMA1P_M1_C11 + t * c;
            c = INV_GAMMA1P_M1_C10 + t * c;
            c = INV_GAMMA1P_M1_C9 + t * c;
            c = INV_GAMMA1P_M1_C8 + t * c;
            c = INV_GAMMA1P_M1_C7 + t * c;
            c = INV_GAMMA1P_M1_C6 + t * c;
            c = INV_GAMMA1P_M1_C5 + t * c;
            c = INV_GAMMA1P_M1_C4 + t * c;
            c = INV_GAMMA1P_M1_C3 + t * c;
            c = INV_GAMMA1P_M1_C2 + t * c;
            c = INV_GAMMA1P_M1_C1 + t * c;
            c = INV_GAMMA1P_M1_C0 + t * c;

            if (x > 0.5) {
                ret = (t / x) * ((c - 0.5) - 0.5);
            }
            else {
                ret = x * c;
            }
        }

        return ret;
    }

    public static double logGammaMinusLogGammaSum(final double a,
                                                  final double b) {
        final double d;
        final double w;
        if (a <= b) {
            d = b + (a - 0.5);
            w = deltaMinusDeltaSum(a, b);

        }
        else {
            d = a + (b - 0.5);
            w = deltaMinusDeltaSum(b, a);

        }

        final double u = d * Math.log1p(a / b);
        final double v = a * (Math.log(b) - 1.0);

        return u <= v ? (w - u) - v : (w - v) - u;

    }

    public static double gamma(final double x) {

        if ((x == Math.rint(x)) && (x <= 0.0)) {
            return Double.NaN;

        }

        final double ret;
        final double absX = Math.abs(x);
        if (absX <= 20.0) {
            if (x >= 1.0) {
                double prod = 1.0;
                double t = x;
                while (t > 2.5) {
                    t -= 1.0;
                    prod *= t;

                }
                ret = prod / (1.0 + invGamma1pm1(t - 1.0));

            }
            else {
                double prod = x;
                double t = x;
                while (t < -0.5) {
                    t += 1.0;
                    prod *= t;
                }
                ret = 1.0 / (prod * (1.0 + invGamma1pm1(t)));
            }
        }
        else {
            final double y = absX + LANCZOS_G + 0.5;
            final double gammaAbs = SQRT_TWO_PI / x *
                    Math.pow(y, absX + 0.5) *
                    Math.exp(-y) * lanczos(absX);
            if (x > 0.0) {
                ret = gammaAbs;
            }
            else {
                ret = -Math.PI /
                        (x * Math.sin(Math.PI * x) * gammaAbs);
            }
        }
        return ret;
    }

    public static double deltaMinusDeltaSum(final double a,
                                            final double b) {

        final double h = a / b;
        final double p = h / (1.0 + h);
        final double q = 1.0 / (1.0 + h);
        final double q2 = q * q;
        final double[] s = new double[DELTA.length];
        s[0] = 1.0;
        for (int i = 1; i < s.length; i++) {
            s[i] = 1.0 + (q + q2 * s[i - 1]);
        }
        final double sqrtT = 10.0 / b;
        final double t = sqrtT * sqrtT;
        double w = DELTA[DELTA.length - 1] * s[s.length - 1];
        for (int i = DELTA.length - 2; i >= 0; i--) {
            w = t * w + DELTA[i] * s[i];
        }
        return w * p / b;
    }

    public static double df(double v1, double v2, double n1, double n2) {
        return (((v1 / n1) + (v2 / n2)) * ((v1 / n1) + (v2 / n2))) /
                ((v1 * v1) / (n1 * n1 * (n1 - 1.0)) + (v2 * v2) /
                        (n2 * n2 * (n2 - 1.0)));
    }

    public static double tTest(EventSummary e1, EventSummary e2) {
        double pValue = tTest(e1.getMean(), e2.getMean(), e1.getVariance(), e2.getVariance(),
                              e1.getN(), e2.getN());
        return pValue;
    }

    public static double tTest(final double m1, final double m2,
                               final double v1, final double v2,
                               final double n1, final double n2) {

        final double t = Math.abs(t(m1, m2, v1, v2, n1, n2));
        final double degreesOfFreedom = df(v1, v2, n1, n2);
        double pValue = 2.0 * cumulativeProbability(-t, degreesOfFreedom);
        return pValue;
    }

    public static double t(final double m1, final double m2,
                           final double v1, final double v2,
                           final double n1, final double n2) {
        return (m1 - m2) / Math.sqrt((v1 / n1) + (v2 / n2));
    }

    public static double cumulativeProbability(double x, double degreesOfFreedom) {
        double ret;
        if (x == 0) {
            ret = 0.5;
        }
        else {
            double t =
                    regularizedBeta(
                            degreesOfFreedom / (degreesOfFreedom + (x * x)),
                            0.5 * degreesOfFreedom,
                            0.5);
            if (x < 0.0) {
                ret = 0.5 * t;
            }
            else {
                ret = 1.0 - 0.5 * t;
            }
        }

        return ret;
    }
}
