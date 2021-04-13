public class EventSummary {
    private double mean;
    private double max;
    private double min;
    private long numberOfElements;
    private double sum;
    private double stdDeviation;
    private double distance;
    private String configuration;

    public EventSummary() {

    }

    public EventSummary(double max, double mean, double min, long numberOfElements,
                        double stdDeviation, double sum, double distance, String configuration) {
        this.max = max;
        this.min = min;
        this.mean = mean;
        this.numberOfElements = numberOfElements;
        this.stdDeviation = stdDeviation;
        this.sum = sum;
        this.distance = distance;
        this.configuration = configuration;
    }

    public String getConfiguration() {
        return configuration;
    }

    public double getDistance() {
        return distance;
    }

    public double getMax() {
        return max;
    }

    public double getMean() {
        return mean;
    }

    public double getMin() {
        return min;
    }

    public long getN() {
        return numberOfElements;
    }

    public double getStandardDeviation() {
        return stdDeviation;
    }

    public double getVariance() {
        return stdDeviation * stdDeviation;
    }

    public double getSum() {
        return sum;
    }
}
