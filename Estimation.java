import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;
import edu.princeton.cs.algs4.Stopwatch;

import java.util.ArrayList;
import java.util.HashMap;

public class Estimation {

    private HashMap<String, EventSummary> events;
    private HashMap<String, HashMap<EventSummary, Integer>> clusters;

    public Estimation(String fileName) {
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
        input.close();
    }

    public Estimation(String allSummary, String clusterSummary) {
        events = new HashMap<String, EventSummary>();
        In input = new In(allSummary);
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
            configuration = configuration.trim();
            EventSummary currentEvent = new EventSummary(max, mean, min, numberOfElements, stdDev,
                                                         sum, distance, configuration);
            events.put(eventName, currentEvent);
        }
        input.close();
        input = new In(clusterSummary);
        clusters = new HashMap<String, HashMap<EventSummary, Integer>>();
        while (!input.isEmpty()) {
            double distance = input.readDouble();
            double max = input.readDouble();
            double min = input.readDouble();
            double mean = input.readDouble();
            double stdDev = input.readDouble();
            long numberOfElements = input.readLong();
            double sum = input.readDouble();
            int weight = input.readInt();
            String configuration = input.readLine();
            configuration = configuration.trim();
            EventSummary currentEvent = new EventSummary(max, mean, min, numberOfElements, stdDev,
                                                         sum, distance, configuration);
            String key = "" + distance + configuration;
            if (clusters.containsKey(key)) {
                clusters.get(key).put(currentEvent, weight);
            }
            else {
                HashMap<EventSummary, Integer> currentList = new HashMap<EventSummary, Integer>();
                currentList.put(currentEvent, weight);
                clusters.put(key, currentList);
            }
        }
        input.close();
    }

    public Estimation(String[] files, Estimation superSet) {
        events = new HashMap<String, EventSummary>();
        for (int i = 0; i < files.length; i++) {
            events.put(files[i], superSet.events.get(files[i]));
        }
    }

    public double meanPValue(EventSummary givenEvent, ArrayList<EventSummary> dataSet) {
        double sum = 0.0;
        for (int i = 0; i < dataSet.size(); i++) {
            EventSummary currentEvent = dataSet.get(i);
            double currentPValue = Statistics
                    .tTest(givenEvent.getMean(), currentEvent.getMean(), givenEvent.getVariance(),
                           currentEvent.getVariance(), givenEvent.getN(), currentEvent.getN());
            sum = sum + currentPValue;
        }
        return sum / dataSet.size();
    }

    public double guessDistance(String fileName) {
        EventSummary givenEvent = events.get(fileName);
        String configuration = givenEvent.getConfiguration();
        ArrayList<EventSummary> distance4 = new ArrayList<EventSummary>();
        ArrayList<EventSummary> distance6 = new ArrayList<EventSummary>();
        ArrayList<EventSummary> distance10 = new ArrayList<EventSummary>();
        ArrayList<EventSummary> distance15 = new ArrayList<EventSummary>();
        for (String file : events.keySet()) {
            EventSummary currentEvent = events.get(file);
            double currentDistance = currentEvent.getDistance();
            String currentConfiguration = currentEvent.getConfiguration();
            if (configuration.equals("-1") || configuration.equals(currentConfiguration)) {
                if (currentDistance == 4) {
                    distance4.add(currentEvent);
                }
                if (currentDistance == 6) {
                    distance6.add(currentEvent);
                }
                if (currentDistance == 10) {
                    distance10.add(currentEvent);
                }
                if (currentDistance == 15) {
                    distance15.add(currentEvent);
                }
            }
        }
        double[] distances = { 4, 6, 10, 15 };
        double[] meanPValues = {
                meanPValue(givenEvent, distance4), meanPValue(givenEvent, distance6),
                meanPValue(givenEvent, distance10), meanPValue(givenEvent, distance15),
                };
        int currentMax = 0;
        for (int i = 0; i < 4; i++) {
            /*if (meanPValues[i] > 0.25) {
                StdOut.println(meanPValues[i]);
            }*/
            if (meanPValues[i] > meanPValues[currentMax]) {
                currentMax = i;
            }
        }
        return distances[currentMax];
    }

    public int[] testGuesses() {
        int[] results = { 0, 0, 0, 0 };
        for (String fileName : events.keySet()) {
            double currentDistance = events.get(fileName).getDistance();
            double guessedDistance = guessDistanceClusterSummary(fileName);
            if (currentDistance == guessedDistance) {
                results[0] = results[0] + 1;
            }
            if (currentDistance < 6.1 && guessedDistance < 6.1) {
                results[1] = results[1] + 1;
            }
            if (currentDistance > 6.1 && guessedDistance > 6.1) {
                results[1] = results[1] + 1;
            }
            if (currentDistance < 6.1 && guessedDistance > 6.1) {
                results[3] = results[3] + 1;
            }
            if (currentDistance > 6.1 && guessedDistance < 6.1) {
                results[2] = results[2] + 1;
            }
        }
        return results;
    }

    public double guessDistance(String fileName, String[] trainSet) {
        EventSummary givenEvent = events.get(fileName);
        String configuration = givenEvent.getConfiguration();
        ArrayList<EventSummary> distance4 = new ArrayList<EventSummary>();
        ArrayList<EventSummary> distance6 = new ArrayList<EventSummary>();
        ArrayList<EventSummary> distance10 = new ArrayList<EventSummary>();
        ArrayList<EventSummary> distance15 = new ArrayList<EventSummary>();
        for (String file : trainSet) {
            EventSummary currentEvent = events.get(file);
            double currentDistance = currentEvent.getDistance();
            String currentConfiguration = currentEvent.getConfiguration();
            if (configuration.equals("-1") || configuration.equals(currentConfiguration)) {
                if (currentDistance == 4) {
                    distance4.add(currentEvent);
                }
                if (currentDistance == 6) {
                    distance6.add(currentEvent);
                }
                if (currentDistance == 10) {
                    distance10.add(currentEvent);
                }
                if (currentDistance == 15) {
                    distance15.add(currentEvent);
                }
            }
        }
        double[] distances = { 4, 6, 10, 15 };
        double[] meanPValues = {
                meanPValue(givenEvent, distance4), meanPValue(givenEvent, distance6),
                meanPValue(givenEvent, distance10), meanPValue(givenEvent, distance15),
                };
        int currentMax = 0;
        for (int i = 0; i < 4; i++) {
            /*if (meanPValues[i] > 0.25) {
                StdOut.println(meanPValues[i]);
            }*/
            if (meanPValues[i] > meanPValues[currentMax]) {
                currentMax = i;
            }
        }
        return distances[currentMax];
    }

    public int[] testRandomDivision() {
        String[] allFiles = new String[events.size()];
        int small = (2 * allFiles.length) / 10;
        int j = 0;
        for (String fileName : events.keySet()) {
            allFiles[j] = fileName;
            j++;
        }
        StdRandom.shuffle(allFiles);
        String[] testFiles = new String[small];
        String[] trainFiles = new String[allFiles.length - small];

        for (int i = 0; i < allFiles.length; i++) {
            if (i < small) {
                testFiles[i] = allFiles[i];
            }
            else {
                trainFiles[i - small] = allFiles[i];
            }
        }

        int[] results = { 0, 0, 0, 0 };
        for (String fileName : testFiles) {
            double currentDistance = events.get(fileName).getDistance();
            double guessedDistance = guessDistance(fileName, trainFiles);
            if (currentDistance == guessedDistance) {
                results[0] = results[0] + 1;
            }
            if (currentDistance < 6.1 && guessedDistance < 6.1) {
                results[1] = results[1] + 1;
            }
            if (currentDistance > 6.1 && guessedDistance > 6.1) {
                results[1] = results[1] + 1;
            }
            if (currentDistance < 6.1 && guessedDistance > 6.1) {
                results[3] = results[3] + 1;
            }
            if (currentDistance > 6.1 && guessedDistance < 6.1) {
                results[2] = results[2] + 1;
            }
        }
        return results;
    }

    public String[] getKeys(HashMap<String, EventSummary> sets, String config) {
        String[] keys = new String[4];
        for (String key : sets.keySet()) {
            //StdOut.println("\n" + key + " Hello\n");
            if (Integer.parseInt(key.substring(0, 1)) == 4 && key.contains(config)) {
                keys[0] = key;
                //StdOut.println("\n" + key + "\n");
            }
            if (Integer.parseInt(key.substring(0, 1)) == 6 && key.contains(config)) {
                keys[1] = key;
                //StdOut.println("\n" + key + "\n");
            }
            if (key.charAt(0) == '1' && key.charAt(1) == '0' && key.contains(config)) {
                keys[2] = key;
                //StdOut.println("\n" + key + "\n");
            }
            if (key.charAt(0) == '1' && key.charAt(1) == '5' && key.contains(config)) {
                keys[3] = key;
                //StdOut.println("\n" + key + "\n");
            }
        }
        return keys;
    }

    public double guessDistance(String fileName, boolean overRide) {
        In input = new In("stationarySummary.tsv");
        HashMap<String, EventSummary> sets = new HashMap<String, EventSummary>();
        while (!input.isEmpty()) {
            double distance = input.readDouble();
            double max = input.readDouble();
            double min = input.readDouble();
            double mean = input.readDouble();
            double stdDev = input.readDouble();
            long numberOfElements = input.readLong();
            double sum = input.readDouble();
            String configuration = input.readLine();
            configuration = configuration.trim();
            EventSummary currentEvent = new EventSummary(max, mean, min, numberOfElements, stdDev,
                                                         sum, distance, configuration);
            String key = distance + configuration;
            //StdOut.println(key);
            sets.put(key, currentEvent);
        }
        /*String givenConfiguration = OrientationDataStats.getConfig(fileName);
        int[] givenData = OrientationDataStats.getRSSIArray(fileName);
        double distance = OrientationDataStats.getDistance(fileName);
        double max = StdStats.max(givenData);
        double min = StdStats.min(givenData);
        double mean = StdStats.mean(givenData);
        double stdDev = StdStats.stddev(givenData);
        long numberOfElements = givenData.length;
        double sum = mean * numberOfElements;
        EventSummary givenEvent = new EventSummary(max, mean, min, numberOfElements, stdDev,
                                                   sum, distance, givenConfiguration);*/
        EventSummary givenEvent = events.get(fileName);
        String givenConfiguration = givenEvent.getConfiguration();
        givenConfiguration = givenConfiguration.trim();
        givenConfiguration = " " + givenConfiguration;
        givenConfiguration = givenConfiguration.trim();
        String[] keys = getKeys(sets, givenConfiguration);
        EventSummary[] events = {
                sets.get(keys[0]), sets.get(keys[1]), sets.get(keys[2]), sets.get(keys[3])
        };
        double[] distances = { 4, 6, 10, 15 };
        double[] pValues = {
                Statistics.tTest(givenEvent, events[0]), Statistics.tTest(givenEvent, events[1]),
                Statistics.tTest(givenEvent, events[2]), Statistics.tTest(givenEvent, events[3]),
                };
        int maxIndex = 0;
        for (int i = 1; i < 4; i++) {
            if (pValues[i] > pValues[maxIndex]) {
                maxIndex = i;
            }
        }
        return distances[maxIndex];
    }

    public double meanPValue(EventSummary givenEvent, HashMap<EventSummary, Integer> cluster) {
        double sum = 0;
        double numberOfElements = 0;
        for (EventSummary e2 : cluster.keySet()) {
            int weight = cluster.get(e2);
            numberOfElements = numberOfElements + weight;
            double currentPValue = Statistics.tTest(givenEvent, e2);
            sum = sum + (weight * currentPValue);
        }
        double meanP = sum / numberOfElements;
        return meanP;
    }

    public double guessDistanceClusterSummary(String fileName) {
        EventSummary givenEvent = events.get(fileName);
        String configuration = givenEvent.getConfiguration();

        double[] distances = { 4.0, 6.0, 10.0, 15.0 };
        double[] meanPValues = new double[4];
        String[] keys = {
                "" + distances[0] + configuration, "" + distances[1] + configuration,
                "" + distances[2] + configuration, "" + distances[3] + configuration,
                };
        for (int i = 0; i < keys.length; i++) {
            meanPValues[i] = meanPValue(givenEvent, clusters.get(keys[i]));
        }

        int currentMax = 0;
        for (int i = 0; i < 4; i++) {
            if (meanPValues[i] > meanPValues[currentMax]) {
                currentMax = i;
            }
        }
        return distances[currentMax];
    }

    public static void main(String[] args) {
        Stopwatch stopwatch = new Stopwatch();
        Estimation estimation;
        if (args.length > 1) {
            estimation = new Estimation(args[0], args[1]);
        }
        else {
            estimation = new Estimation(args[0]);
        }
        StdOut.println("Time taken to construct = " + stopwatch.elapsedTime());
        int[] results = estimation.testGuesses();
        StdOut.println("Number of correct distance guesses = " + results[0]);
        StdOut.println("Number of correct event guesses = " + results[1]);
        StdOut.println("Number of false positives = " + results[2]);
        StdOut.println("Number of false negatives = " + results[3]);
        double den = results[1] + results[2] + results[3];
        double accuracy = results[1] / den;
        StdOut.println("Accuracy = " + accuracy);
        StdOut.println("Time taken for the first test = " + stopwatch.elapsedTime() + "\n");

        /*int[] randomDivisionResults = estimation.testRandomDivision();
        StdOut.println("Number of correct distance guesses = " + randomDivisionResults[0]);
        StdOut.println("Number of correct event guesses = " + randomDivisionResults[1]);
        StdOut.println("Number of false positives = " + randomDivisionResults[2]);
        StdOut.println("Number of false negatives = " + randomDivisionResults[3]);
        den = randomDivisionResults[1] + randomDivisionResults[2] + randomDivisionResults[3];
        accuracy = randomDivisionResults[1] / den;
        StdOut.println("Accuracy = " + accuracy);
        StdOut.println("Time taken for the first test = " + stopwatch.elapsedTime());*/
    }
}
