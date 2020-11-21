import static java.lang.Math.*;

/**
 * @see http://docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite2/JHS197_liite2.html
 * @see http://docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite3/JHS197_liite3.html
 */
public class EUREFFIN {

    private static final double a = 6378137.0;
    private static final double f = 1.0 / 298.257222101;
    private static final double e = sqrt(2 * f - (f * f));

    private static final double k0 = 0.9996;
    private static final double l0 = toRadians(27.0);
    private static final double E0 = 500000.0;

    private static final double n = f / (2.0 - f);
    private static final double n2 = n * n;
    private static final double n3 = n * n * n;
    private static final double n4 = n * n * n * n;

    private static final double A1 = (a / (1.0 + n)) * (1.0 + (n2 / 4.0) + (n4 / 64.0));

    private static final double h1 = (n / 2.0) - (2.0 * n2 / 3.0) + (37.0 * n3 / 96.0) - (n4 / 360.0);
    private static final double h2 = (n2 / 48.0) + (n3 / 15.0) - (437.0 * n4 / 1440.0);
    private static final double h3 = (17.0 * n3 / 480.0) - (37.0 * n4 / 840.0);
    private static final double h4 = (4397.0 * n4 / 161280.0);

    private static final double h1_ = (n / 2.0) - (2.0 * n2 / 3.0) + (5.0 * n3 / 16.0) + (41.0 * n4 / 180.0);
    private static final double h2_ = (13.0 * n2 / 48.0) - (3.0 * n3 / 5.0) + (557.0 * n4 / 1440.0);
    private static final double h3_ = (61.0 * n3 / 240.0) - (103.0 * n4 / 140.0);
    private static final double h4_ = (49561.0 * n4 / 161280.0);

    private static final double epsilon = 1e-12;

    public static void toETRSTM35FIN(double lon, double lat, double[] out, int off) {
        lat = toRadians(lat);
        lon = toRadians(lon);

        double Q1 = asinh(tan(lat));
        double Q2 = atanh(e * sin(lat));
        double Q = Q1 - (e * Q2);

        double l = lon - l0;
        double B = atan(sinh(Q));
        double n_ = atanh(cos(B) * sin(l));

        double ks_ = asin(sin(B) / sech(n_));

        double ks1 = h1_ * sin(2.0 * ks_) * cosh(2.0 * n_);
        double ks2 = h2_ * sin(4.0 * ks_) * cosh(4.0 * n_);
        double ks3 = h3_ * sin(6.0 * ks_) * cosh(6.0 * n_);
        double ks4 = h4_ * sin(8.0 * ks_) * cosh(8.0 * n_);

        double n1_ = h1_ * cos(2.0 * ks_) * sinh(2.0 * n_);
        double n2_ = h2_ * cos(4.0 * ks_) * sinh(4.0 * n_);
        double n3_ = h3_ * cos(6.0 * ks_) * sinh(6.0 * n_);
        double n4_ = h4_ * cos(8.0 * ks_) * sinh(8.0 * n_);

        double ks = ks_ + ks1 + ks2 + ks3 + ks4;
        double nn = n_ + n1_ + n2_ + n3_ + n4_;

        double N = A1 * ks * k0;
        double E = A1 * nn * k0 + E0;

        out[off + 0] = E;
        out[off + 1] = N;
    }

    public static void toETRSGeodetic(double E, double N, double[] out, int off) {
        double ks = N / (A1 * k0);
        double nn = (E - E0) / (A1 * k0);

        double ks1_ = h1 * sin(2.0 * ks) * cosh(2.0 * nn);
        double ks2_ = h2 * sin(4.0 * ks) * cosh(4.0 * nn);
        double ks3_ = h3 * sin(6.0 * ks) * cosh(6.0 * nn);
        double ks4_ = h4 * sin(8.0 * ks) * cosh(8.0 * nn);

        double nn1_ = h1 * cos(2.0 * ks) * sinh(2.0 * nn);
        double nn2_ = h2 * cos(4.0 * ks) * sinh(4.0 * nn);
        double nn3_ = h3 * cos(6.0 * ks) * sinh(6.0 * nn);
        double nn4_ = h4 * cos(8.0 * ks) * sinh(8.0 * nn);

        double ks_ = ks - ks1_ - ks2_ - ks3_ - ks4_;
        double nn_ = nn - nn1_ - nn2_ - nn3_ - nn4_;

        double B = asin(sech(nn_) * sin(ks_));
        double l = asin(tanh(nn_) / cos(B));

        double Q = asinh(tan(B));
        double Q1 = Q + e * atanh(e * tanh(Q));

        double delta;
        do {
            double Q2 = Q + e * atanh(e * tanh(Q1));
            delta = Q2 - Q1;
            Q1 = Q2;
        } while (Math.abs(delta) > epsilon);

        double lat = atan(sinh(Q1));
        double lon = l0 + l;

        out[off + 0] = toDegrees(lon);
        out[off + 1] = toDegrees(lat);
    }

    private static double asinh(double x) {
        return Math.log(x + Math.sqrt(x*x + 1.0));
    }

    private static double sech(double x) {
        return 1.0 / Math.cosh(x);
    }

    private static double atanh(double x) {
        return 0.5 * Math.log((1.0 + x) / (1.0 - x));
    }

    public static void main(String[] args) {
        double[] tmp = new double[2];
        double E = 106256.35958;
        double N = 6715706.37705;
        toETRSGeodetic(E, N, tmp, 0);
        System.out.printf("%.7f %.7f%n", E, N);
        double lon = tmp[0];
        double lat = tmp[1];
        toETRSTM35FIN(lon, lat, tmp, 0);
        System.out.printf("%.7f %.7f%n", lon, lat);
        double E_ = tmp[0];
        double N_ = tmp[1];
        System.out.printf("%.7f %.7f%n", E_, N_);
    }

}
