/* *****************************************************************************
 *  Name:    Vandan Mehta
 *
 *  Description:  Performs statisticcal operations on NIST dataset for TC4TL
 *                challange.
 *
 **************************************************************************** */

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdStats;
import org.apache.commons.math3.stat.inference.TTest;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class DataStats {

    /**
     * The below device setting is the most frequent (321/935 files) in the training data and hence
     * it is set as a criteria when we want to calculate data while keeping the device
     * configurations same.
     */
    public static final String SELECTED_DEVICES
            = " TXDevice,iPhone11 TXPower,8 RXDevice,iPhone11_Pro_max";

    /**
     * A helper function to check if the device configuration is the same as 'SELECTED_DEVICES'
     *
     * @param fileName is the name of the file
     * @return True if the device configuration in @param fileName is same as 'SELECTED_DEVICES',
     * false otherwise.
     */
    public static boolean isSelected(String fileName) {
        In input = new In("dataFiles\\" + fileName);
        input.readLine();
        String deviceData = "";
        for (int i = 0; i < 3; i++) {
            String line = input.readString();
            deviceData = deviceData + " " + line;
        }
        input.close();
        return deviceData.equals(SELECTED_DEVICES);
    }

    /**
     * Returns the Std Deviation of the RSSI values in the given file.
     *
     * @param filename is the name of the file.
     * @return The Std Deviation of the RSSI values in the file 'fileName'.
     */
    public static double stdDev(String filename) {
        int[] rssiArray = getRSSIArray(filename);
        return StdStats.stddev(rssiArray);
    }

    /**
     * Returns the Std Deviation of the RSSI values in the given file.
     *
     * @param filename is the name of the file.
     * @return The Std Deviation of the RSSI values in the file 'fileName'.
     */
    public static double mean(String filename) {
        int[] rssiArray = getRSSIArray(filename);
        return StdStats.mean(rssiArray);
    }

    /**
     * Returns the RSSI values stored in the given file in and int array.
     *
     * @param filename: Name of the file.
     * @return The array containing the RSSI values stored in the file 'fileName'.
     */
    public static int[] getRSSIArray(String filename) {
        In input = new In("dataFiles\\" + filename);
        input.readString();
        ArrayList<Integer> rssiValues = new ArrayList<Integer>();

        input.readString();
        input.readString();
        input.readString();
        while (input.hasNextLine()) {
            String line = input.readLine();
            String[] parts = line.trim().split(",");
            for (int i = 0; i < parts.length; i++) {
                if (i == 2) {
                    rssiValues.add(Integer.parseInt(parts[i]));
                    //StdOut.println(parts[i]);
                }
            }
        }
        int[] rssiArray = new int[rssiValues.size()];
        for (int i = 0; i < rssiArray.length; i++) {
            rssiArray[i] = rssiValues.get(i);
        }
        return rssiArray;
    }

    /**
     * The file EventsAndDistances.tsv which contains the names of all data files followed by their
     * respective distances in feet is redirected to Std Input in order to run this function.
     * <p>
     * This function opens every data file and calculates their mean and prints the Highest Mean
     * found among all the files and the distance between the 2 devices in the event considered in
     * that file.
     */
    public static void highestMean() {
        double currentHighest = Double.NEGATIVE_INFINITY;
        double currentDistance = 0;
        StdIn.readLine();
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            double distance = StdIn.readDouble();
            if (isSelected(word) && distance == 15) {
                double currentMean = mean(word);
                if (currentMean > currentHighest) {
                    currentDistance = distance;
                    currentHighest = currentMean;
                }
            }
        }

        StdOut.println("HighestMean = " + currentHighest + " And Distance = " + currentDistance);
    }

    /**
     * This file calculates the mean RSSI values of all files where the devices are kept at
     * 'desiredDistance'. It prints the mean RSSI value of all these means.
     *
     * @param desiredDistance The distance at which we want to calculate the mean of Means.
     */
    public static void meanOfMeans(double desiredDistance) {
        double currentSum = 0;
        double numberOfItems = 0;
        StdIn.readLine();
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            double distance = StdIn.readDouble();
            if (isSelected(word) && distance == desiredDistance) {
                currentSum = currentSum + mean(word);
                numberOfItems++;
            }
        }
        double average = currentSum / numberOfItems;

        StdOut.println("Mean of Means = " + average + " And Distance = " + desiredDistance);
    }

    /**
     * A helper function used to print the different device configurations in  all the files along
     * with their frequencies.
     */
    public static void deviceFreq() {
        HashMap<String, Integer> devices = new HashMap<String, Integer>();
        StdIn.readLine();
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            double distance = StdIn.readDouble();
            // StdOut.println(word + "\t" + distance);
            In input = new In("dataFiles\\" + word);
            input.readLine();
            String deviceData = "";
            for (int i = 0; i < 3; i++) {
                String line = input.readString();
                deviceData = deviceData + " " + line;
            }
            input.close();
            if (devices.containsKey(deviceData)) {
                int currentFreq = devices.get(deviceData);
                currentFreq++;
                devices.put(deviceData, currentFreq);
            }
            else {
                devices.put(deviceData, 1);
            }

        }
        int maxValue = 0;
        for (Map.Entry<String, Integer> mapElement : devices.entrySet()) {
            String key = mapElement.getKey();
            int value = mapElement.getValue();
            if (value > maxValue) {
                StdOut.println(key + "\t" + value);
                maxValue = value;
            }
        }
    }

    /**
     * Given 2 files, it performs a T-test where the sets being comppared are the sets of RSSi
     * values of each file. It returns the probaility whether these 2 data sets belong to the same
     * distribution or not.
     *
     * @param file1 Name of the first file containing the RSSI values.
     * @param file2 Name of the second file containing the RSSI values.
     * @return The probability of the event that both the sets belong to the same distribution.
     */
    public static double tTest(String file1, String file2) {
        int[] rssiCalculated = getRSSIArray(file1);
        int[] rssiExpected = getRSSIArray(file2);
        double[] set1 = new double[rssiCalculated.length];
        double[] set2 = new double[rssiExpected.length];
        for (int i = 0; i < set1.length; i++) {
            set1[i] = rssiCalculated[i];
        }

        for (int i = 0; i < rssiExpected.length; i++) {
            set2[i] = rssiExpected[i];
        }
        TTest tTest = new TTest();
        double pValue = tTest.tTest(set1, set2);
        return pValue;
    }

    /**
     * Returns the set of all file names where the distance between the devices is @param distance.
     *
     * @param distance The desired distance between the 2 devices.
     * @return The set of file names where the distance between the devices is @param distance.
     */
    public static String[] sameDistanceFiles(double distance) {
        ArrayList<String> fileList = new ArrayList<String>();
        In input = new In("EventsAndDistances.tsv");
        input.readLine();
        int j = 0;
        while (j < 935 && input.hasNextChar()) {
            String fileName = input.readString();
            double fileDistance = input.readDouble();
            if (fileDistance == distance && isSelected(fileName)) {
                fileList.add(fileName);
            }
            j++;
        }
        input.close();
        String[] files = new String[fileList.size()];
        for (int i = 0; i < files.length; i++) {
            files[i] = fileList.get(i);
        }
        return files;
    }

    /**
     * After collecting all the files where the distance between the 2 devices is @param distance,
     * it considers all pairs of files (one pair at a time) and computes the p-values for the T-test
     * between their RSSI data sets. Finally it returns the mean of these choose(n,2) p-values where
     * n is the number of files where the distance between the devices is @param distance.
     *
     * @param distance The desired distance between the 2 devices.
     * @return The average p-value after performing t-tests on each pair of files where the distance
     * between the devices is @param distance.
     */
    public static double meanPValues(double distance) {
        String[] fileNames = sameDistanceFiles(distance);
        double pValueSum = 0;
        double numberOfEntries = 0;
        double maxPValue = 0;
        for (int i = 0; i < fileNames.length; i++) {
            for (int j = i + 1; j < fileNames.length; j++) {

                double currentPValue = tTest(fileNames[i], fileNames[j]);
                pValueSum = pValueSum + currentPValue;
                numberOfEntries++;
                /* if (numberOfEntries % 10000 == 0) {
                    StdOut.println("Number of Entries = " + numberOfEntries + " current sum = "
                                           + pValueSum);
                }*/
                if (currentPValue > maxPValue) {
                    maxPValue = currentPValue;
                }

            }
        }
        StdOut.println("Number of files with distance " + distance + " = " + fileNames.length);
        StdOut.println("Maximum PValue = " + maxPValue);
        double meanPvalue = pValueSum / numberOfEntries;
        return meanPvalue;
    }

    /**
     * Here we collect 2 groups of files. Group I consists of files with distance = @param distance1
     * and group II consists of files with distance = @param distance2. For each file in group I, we
     * perform a t-test and compute the p-value with every file in group II. Finally we return the
     * mean of these m*n p-values where 'm' is the number of files in group I and 'n' is the number
     * of files in group II.
     *
     * @param distance1 The desired distance between devices in the first group of files.
     * @param distance2 The desired distance between devices in the second group of files.
     * @return The mean of the p-values computed in the m*n t-tests.
     */
    public static double differentDistancesMeanPValues(double distance1, double distance2) {
        String[] fileNames1 = sameDistanceFiles(distance1);
        String[] fileNames2 = sameDistanceFiles(distance2);
        double pValueSum = 0;
        double numberOfEntries = 0;
        double maxPValue = 0;
        String currentMaxFile1 = "";
        String currentMaxFile2 = "";
        for (int i = 0; i < fileNames1.length; i++) {
            for (int j = 0; j < fileNames2.length; j++) {
                double currentPValue = tTest(fileNames1[i], fileNames2[j]);
                pValueSum = pValueSum + currentPValue;
                numberOfEntries++;
                /* if (numberOfEntries % 20000 == 0) {
                    StdOut.println("Number of Entries = " + numberOfEntries + " current sum = "
                                           + pValueSum);
                }*/
                if (currentPValue > maxPValue) {
                    maxPValue = currentPValue;
                    currentMaxFile1 = fileNames1[i];
                    currentMaxFile2 = fileNames2[j];
                }
            }
        }
        StdOut.println("Number of files with distances " + distance1 + " and " + distance2 + " = "
                               + fileNames1.length + " and " + fileNames2.length + " respectively");
        StdOut.println("Maximum PValue = " + maxPValue + " Files: " + currentMaxFile1 + " and "
                               + currentMaxFile2);
        double meanPvalue = pValueSum / numberOfEntries;
        return meanPvalue;
    }

    public static void main(String[] args) {
        double distance = Double.parseDouble(args[0]);
        double averagePValue = meanPValues(distance);
        StdOut.println("Average P-Value for distance " + distance + " = " + averagePValue);
        double distance2 = Double.parseDouble(args[1]);
        double averagePValue2 = meanPValues(distance2);
        StdOut.println("Average P-Value for distance " + distance2 + " = " + averagePValue2);
        double diffAvg = differentDistancesMeanPValues(distance, distance2);
        StdOut.println("Average P-Value for distances " + distance + " and " + distance2 +
                               " = " + diffAvg);
    }
}
