import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdStats;

public class FileWriter {

    public static boolean shouldCopy(String line, int i) {
        if (i < 7) {
            return true;
        }
        else {
            return line.contains("Bluetooth");
        }
    }

    public static boolean isSimple(String fileName) {
        In input = new In("test\\" + fileName);
        String file = input.readAll();
        return (!file.contains("walking"));
    }

    public static void writeStats(String fileName) {
        int[] data = OrientationDataStats.getRSSIArray(fileName);
        double mean = StdStats.mean(data);
        double stdDev = StdStats.stddev(data);
        double max = StdStats.max(data);
        double min = StdStats.min(data);
        long numberOfElements = data.length;
        double sum = mean * numberOfElements;
        double distance = OrientationDataStats.getDistance(fileName);
        String configuration = OrientationDataStats.getConfig(fileName);
        StdOut.println(
                fileName + "\t" + distance + "\t" + max + "\t" + min + "\t" + mean + "\t" + stdDev
                        + "\t" + numberOfElements + "\t" + sum + "\t" + configuration);
    }

    public static int[] merge(int[] a, int[] b) {
        int[] merged = new int[a.length + b.length];
        for (int i = 0; i < merged.length; i++) {
            if (i < a.length) {
                merged[i] = a[i];
            }
            else {
                merged[i] = b[i - a.length];
            }
        }
        return merged;
    }

    public static void writeClusterSummary(String[] cluster) {
        int weight = cluster.length;
        double distance = OrientationDataStats.getDistance(cluster[0]);
        String config = OrientationDataStats.getConfig(cluster[0]);
        int[] data = OrientationDataStats.getRSSIArray(cluster[0]);
        for (int i = 1; i < weight; i++) {
            data = merge(data, OrientationDataStats.getRSSIArray(cluster[i]));
        }
        double mean = StdStats.mean(data);
        double stdDev = StdStats.stddev(data);
        double max = StdStats.max(data);
        double min = StdStats.min(data);
        long numberOfElements = data.length;
        double sum = mean * numberOfElements;
        StdOut.println(distance + "\t" + max + "\t" + min + "\t" + mean + "\t" + stdDev
                               + "\t" + numberOfElements + "\t" + sum + "\t" + weight + "\t"
                               + config);
    }

    public static void main(String[] args) {
        /*HashMap<String, ArrayList<String>> distributions = new HashMap<String, ArrayList<String>>();
        StdIn.readLine();
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            StdIn.readDouble();
            String key = OrientationDataStats.getConfig(word) + " " + OrientationDataStats
                    .getDistance(word);
            if (distributions.containsKey(key)) {
                distributions.get(key).add(word);
            }
            else {
                ArrayList<String> files = new ArrayList<String>();
                files.add(word);
                distributions.put(key, files);
            }
        }
        for (String sets : distributions.keySet()) {
            ArrayList<String> currentSet = distributions.get(sets);
            int[] data = OrientationDataStats.getRSSIArray(currentSet.get(0));
            double distance = OrientationDataStats.getDistance(currentSet.get(0));
            String configuration = OrientationDataStats.getConfig(currentSet.get(0));
            for (int i = 1; i < currentSet.size(); i++) {
                data = merge(data, OrientationDataStats.getRSSIArray(currentSet.get(i)));
            }
            double mean = StdStats.mean(data);
            double stdDev = StdStats.stddev(data);
            double max = StdStats.max(data);
            double min = StdStats.min(data);
            long numberOfElements = data.length;
            double sum = mean * numberOfElements;
            StdOut.println(
                    distance + "\t" + max + "\t" + min + "\t" + mean + "\t" + stdDev
                            + "\t" + numberOfElements + "\t" + sum + "\t" + configuration);
        }*/
        while (!StdIn.isEmpty()) {
            StdIn.readLine();
            String line = StdIn.readLine();
            line = line.trim();
            String[] cluster = line.split(" ");
            /*for (int i = 0; i < cluster.length; i++) {
                StdOut.println(cluster[i]);
            }*/
            writeClusterSummary(cluster);
        }

    }
}
