import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;
import edu.princeton.cs.algs4.Stopwatch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class Clusters {

    private class Roots implements Comparable<Roots> {
        public String fileName;
        public int numberOfChildren;

        public Roots(String fileName, int numberOfChildren) {
            this.fileName = fileName;
            this.numberOfChildren = numberOfChildren;
        }

        @Override
        public int compareTo(Roots other) {
            return other.numberOfChildren - this.numberOfChildren;
        }
    }

    public HashMap<String, EventSummary> events;
    private HashMap<String, HashMap<String, Double>> records;
    private Roots[] setOfRoots;
    private HashMap<String, Integer> isAvailable;
    public HashMap<String, Integer> weights;
    private boolean isSecondConstructor;

    public Clusters(String fileName) {
        isSecondConstructor = false;
        events = new HashMap<String, EventSummary>();
        In input = new In(fileName);
        while (!input.isEmpty()) {
            String eventName = input.readString();
            double distance = input.readDouble();
            double max = input.readDouble();
            double min = input.readDouble();
            double mean = input.readDouble();
            double stdDev = input.readDouble();
            long numberOfElements = input.readLong();
            double sum = input.readDouble();
            String configuration = input.readLine();
            EventSummary currentEvent = new EventSummary(max, mean, min, numberOfElements, stdDev,
                                                         sum, distance, configuration);
            events.put(eventName, currentEvent);
        }
        this.records = new HashMap<String, HashMap<String, Double>>();
    }

    public Clusters(String fileName, boolean override) {
        events = new HashMap<String, EventSummary>();
        weights = new HashMap<String, Integer>();
        In input = new In(fileName);
        this.isSecondConstructor = override;
        int i = 0;
        while (!input.isEmpty()) {
            String eventName = Integer.toString(i);
            double distance = input.readDouble();
            double max = input.readDouble();
            double min = input.readDouble();
            double mean = input.readDouble();
            double stdDev = input.readDouble();
            long numberOfElements = input.readLong();
            double sum = input.readDouble();
            int weight = input.readInt();
            String configuration = input.readLine();
            EventSummary currentEvent = new EventSummary(max, mean, min, numberOfElements, stdDev,
                                                         sum, distance, configuration);
            events.put(eventName, currentEvent);
            weights.put(eventName, weight);
            i++;
        }
        this.records = new HashMap<String, HashMap<String, Double>>();
    }

    public double constructRecords(double threshold) {
        int length = events.size() - 1;
        double sum = 0;
        for (String file1 : events.keySet()) {
            HashMap<String, Double> accepted = new HashMap<String, Double>();
            EventSummary event1 = events.get(file1);
            for (String file2 : events.keySet()) {
                /*if (file1.equals(file2)) {
                    continue;
                }*/
                EventSummary event2 = events.get(file2);
                if (event1.getConfiguration().equals(event2.getConfiguration())
                        && event1.getDistance() == event2.getDistance()) {
                    double currentPValue = Statistics.tTest(event1, event2);
                    if (currentPValue > threshold) {
                        accepted.put(file2, currentPValue);
                    }
                }
            }
            if (accepted.size() > 0) {
                //StdOut.println();
                //StdOut.println("Cluster Size of " + file1 + " = " + accepted.size());
                //StdOut.println();
                records.put(file1, accepted);
                sum = sum + accepted.size();
            }
        }
        double averageAcceptedSize = sum / records.size();
        setOfRoots = new Roots[records.size()];
        int j = 0;
        for (String key : records.keySet()) {
            setOfRoots[j] = new Roots(key, records.get(key).size());
            j++;
        }
        Arrays.sort(setOfRoots);
        return averageAcceptedSize;
    }

    public int getWeight(ArrayList<String> summaries) {
        int weight = 0;
        for (int i = 0; i < summaries.size(); i++) {
            weight = weight + weights.get(summaries.get(i));
        }
        return weight;
    }

    public String[] constructCluster(String fileName, int trials) {
        HashMap<String, Double> currentRoot = records.get(fileName);
        String[] files = new String[currentRoot.size()];
        String[] maxCluster = new String[1];
        int maxWeight = 0;
        int j = 0;
        for (String key : currentRoot.keySet()) {
            files[j] = key;
            j++;
        }
        for (int i = 0; i < trials; i++) {
            ArrayList<String> currentCluster = new ArrayList<String>();
            currentCluster.add(fileName);
            int chosen = 0;
            for (int k = 0; k < files.length; k++) {
                int index = StdRandom.uniform(0, files.length - chosen);
                if (!isAvailable.containsKey(files[index])) {
                    chosen++;
                    String tmp = files[index];
                    files[index] = files[files.length - chosen];
                    files[files.length - chosen] = tmp;
                    continue;
                }
                HashMap<String, Double> currentNode = records.get(files[index]);
                boolean isClique = true;
                for (int l = 0; l < currentCluster.size(); l++) {
                    if (!currentNode.containsKey(currentCluster.get(l))) {
                        isClique = false;
                    }
                }
                if (isClique) {
                    currentCluster.add(files[index]);
                }
                chosen++;
                String tmp = files[index];
                files[index] = files[files.length - chosen];
                files[files.length - chosen] = tmp;
            }
            if (isSecondConstructor) {
                int currentWeight = getWeight(currentCluster);
                if (currentWeight > maxWeight) {
                    maxWeight = currentWeight;
                    maxCluster = new String[currentCluster.size()];
                    for (int l = 0; l < maxCluster.length; l++) {
                        maxCluster[l] = currentCluster.get(l);
                    }
                }
            }
            else {
                if (currentCluster.size() >= maxCluster.length) {
                    maxCluster = new String[currentCluster.size()];
                    for (int l = 0; l < maxCluster.length; l++) {
                        maxCluster[l] = currentCluster.get(l);
                    }
                }
            }
        }
        for (int i = 0; i < maxCluster.length; i++) {
            isAvailable.remove(maxCluster[i]);
        }
        return maxCluster;
    }

    public String[][] constructClusters(int trials) {
        isAvailable = new HashMap<String, Integer>();
        for (int i = 0; i < setOfRoots.length; i++) {
            isAvailable.put(setOfRoots[i].fileName, 0);
        }
        ArrayList<String[]> clusters = new ArrayList<String[]>();
        for (int i = 0; i < setOfRoots.length; i++) {
            if (isAvailable.containsKey(setOfRoots[i].fileName)) {
                clusters.add(constructCluster(setOfRoots[i].fileName, trials));
            }
        }
        String[][] allClusters = new String[clusters.size()][];
        for (int i = 0; i < allClusters.length; i++) {
            allClusters[i] = clusters.get(i);
        }
        return allClusters;
    }

    public static void main(String[] args) {
        Stopwatch stopwatch = new Stopwatch();
        Clusters clusters = new Clusters(args[0]);
        //StdOut.println("Constructed");
        double averageAcceptance = clusters.constructRecords(0.2);
        //StdOut.println(averageAcceptance);
        int trials = Integer.parseInt(args[1]);
        String[][] allClusters = clusters.constructClusters(trials);
        double sum = 0.0;
        int[] clusterStats = new int[5];
        for (int i = 0; i < allClusters.length; i++) {
            sum = sum + allClusters[i].length;
            StdOut.println(allClusters[i].length);
            for (int k = 0; k < clusterStats.length; k++) {
                if (allClusters[i].length - 1 > (25 * (k + 1))) {
                    clusterStats[k] = clusterStats[k] + 1;
                }
            }
            for (int j = 0; j < allClusters[i].length; j++) {
                StdOut.print(allClusters[i][j] + " ");
            }
            StdOut.println();
        }
        double averageClusterSize = sum / allClusters.length;
        /*StdOut.println("Average Cluster Size = " + averageClusterSize);
        StdOut.println("Number of Cluster = " + allClusters.length);
        for (int i = 0; i < clusterStats.length; i++) {
            int number = 25 * (i + 1);
            StdOut.println("The number of clusters with size more than " + number + " = "
                                   + clusterStats[i]);
        }
        StdOut.println("Time taken = " + stopwatch.elapsedTime());*/
    }
}
