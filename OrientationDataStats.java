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

public class OrientationDataStats {

    public static int weirdCount = 0;

    /**
     * The below device setting is the most frequent (321/935 files) in the training data and hence
     * it is set as a criteria when we want to calculate data while keeping the device
     * configurations same.
     */
    public static final String[] CONFIGURATIONS = {
            " TXDevice,iPhone11 TXPower,8 RXDevice,iPhone11_Pro_max",
            " TXDevice,iPhone6s TXPower,12 RXDevice,iPhone6s",
            " TXDevice,iPhoneXR TXPower,12 RXDevice,iPhoneXR",
            " TXDevice,iPhoneXS TXPower,12 RXDevice,iPhone8",
            " TXDevice,iPhoneXR TXPower,12 RXDevice,iPhone11_Pro",
            " TXDevice,iPhone11_Pro TXPower,7 RXDevice,iPhoneXR"
    };

    public static final String[] ORIENTATIONS = {
            " TXCarry,unknown RXCarry,unknown", " TXCarry,pocket RXCarry,hand",
            " TXCarry,hand RXCarry,pocket", " TXCarry,pocket RXCarry,pocket",
            " TXCarry,hand RXCarry,hand"
    };

    /**
     * A helper function to check if the device configuration is the same as 'SELECTED_DEVICES'
     *
     * @param fileName is the name of the file
     * @return True if the device configuration in @param fileName is same as 'SELECTED_DEVICES',
     * false otherwise.
     */
    public static boolean isSelected(String fileName, int config, int orientation) {
        boolean selection = true;
        String deviceData = getConfig(fileName);
        String deviceOrientation = getOrientation(fileName);
        if (orientation >= 0) {
            selection = deviceOrientation.equals(ORIENTATIONS[orientation]);
        }
        if (config >= 0) {
            selection = selection && deviceData.equals(CONFIGURATIONS[config]);
        }
        return selection;
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

    public static double getDistance(String filename) {
        In input = new In("refinedTestData\\" + filename);
        String line = input.readString();
        String[] parts = line.split(",");
        double distance = Double.parseDouble(parts[1]);
        input.close();
        return distance;
    }

    public static String getConfig(String filename) {
        In input = new In("refinedTestData\\" + filename);
        input.readLine();
        /* */
        input.readLine();
        String deviceConfig = input.readLine();
        /* */
        //String deviceConfig = "";

        /*for (int i = 0; i < 3; i++) {
            deviceConfig = deviceConfig + " " + input.readLine();
        }*/
        input.close();
        return deviceConfig;
    }

    public static String getOrientation(String filename) {
        In input = new In("refinedTestData\\" + filename);
        for (int i = 0; i < 4; i++) {
            input.readLine();
        }
        String orientation = "";
        for (int i = 0; i < 2; i++) {
            orientation = orientation + " " + input.readLine();
        }
        input.close();
        return orientation;
    }

    public static int[] getData(String fileName) {
        int[] data = new int[3];
        data[0] = (int) (getDistance(fileName));
        data[1] = -1;
        data[2] = -1;
        String config = getConfig(fileName);
        for (int i = 0; i < CONFIGURATIONS.length; i++) {
            if (config.equals(CONFIGURATIONS[i])) {
                data[1] = i;
                break;
            }
        }
        String orientation = getOrientation(fileName);
        for (int i = 0; i < ORIENTATIONS.length; i++) {
            if (orientation.equals(ORIENTATIONS[i])) {
                data[2] = i;
                break;
            }
        }
        return data;
    }

    /**
     * Returns the RSSI values stored in the given file in and int array.
     *
     * @param filename: Name of the file.
     * @return The array containing the RSSI values stored in the file 'fileName'.
     */
    public static int[] getRSSIArray(String filename) {
        In input = new In("refinedTestData\\" + filename);

        ArrayList<Integer> rssiValues = new ArrayList<Integer>();

        for (int i = 0; i < 8; i++) {
            input.readLine();
        }
        boolean isWeird = false;
        while (input.hasNextLine()) {
            String line = input.readLine();
            String[] parts = line.trim().split(",");
            for (int i = 0; i < parts.length; i++) {
                if (i == 2) {
                    int rssi = Integer.parseInt(parts[i]);
                    if (rssi < 30) {
                        rssiValues.add(Integer.parseInt(parts[i]));
                        //StdOut.println(parts[i]);
                    }
                    else {
                        isWeird = true;
                    }
                }
            }
        }
        /*if (isWeird) {
            StdOut.println("File " + filename + " has some positive rssi values.");
            weirdCount++;
        }*/
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
    public static void highestMean(int config, int orientations) {
        double currentHighest = Double.NEGATIVE_INFINITY;
        double currentDistance = 0;
        String highestFile = "";
        StdIn.readLine();
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            double distance = StdIn.readDouble();
            if (isSelected(word, config, orientations)) {
                double currentMean = mean(word);
                if (currentMean > currentHighest) {
                    currentDistance = distance;
                    currentHighest = currentMean;
                    highestFile = word;
                }
            }
        }

        StdOut.println("Number of files with positive rssi values = " + weirdCount);
        StdOut.println("HighestMean = " + currentHighest + " And Distance = " + currentDistance
                               + " And filename = " + highestFile);
    }

    /**
     * This file calculates the mean RSSI values of all files where the devices are kept at
     * 'desiredDistance'. It prints the mean RSSI value of all these means.
     *
     * @param desiredDistance The distance at which we want to calculate the mean of Means.
     */
    public static void meanOfMeans(double desiredDistance, int config, int orientation) {
        double currentSum = 0;
        double numberOfItems = 0;
        StdIn.readLine();
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            double distance = StdIn.readDouble();
            if (isSelected(word, config, orientation) && distance == desiredDistance) {
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
        int count = 0;
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            double distance = StdIn.readDouble();
            // StdOut.println(word + "\t" + distance);
            String deviceData = getConfig(word);
            if (devices.containsKey(deviceData)) {
                int currentFreq = devices.get(deviceData);
                currentFreq++;
                devices.put(deviceData, currentFreq);
            }
            else {
                devices.put(deviceData, 1);
            }
            count++;

        }
        int maxValue = 0;
        for (Map.Entry<String, Integer> mapElement : devices.entrySet()) {
            String key = mapElement.getKey();
            int value = mapElement.getValue();
            if (value > maxValue) {
                maxValue = value;
            }
            StdOut.println(key + "\t" + value);
        }
        StdOut.println("Number of files = " + count);
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
    public static String[] sameDistanceFiles(double distance, int config, int orientation) {
        ArrayList<String> fileList = new ArrayList<String>();
        In input = new In("EventsAndDistances.tsv");
        input.readLine();
        int j = 0;
        while (j < 935 && input.hasNextChar()) {
            String fileName = input.readString();
            double fileDistance = input.readDouble();
            if (fileDistance == distance && isSelected(fileName, config, orientation)) {
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
    public static double meanPValues(double distance, int config, int orientation) {
        String[] fileNames = sameDistanceFiles(distance, config, orientation);
        double pValueSum = 0;
        double numberOfEntries = 0;
        double maxPValue = 0;
        for (int i = 0; i < fileNames.length; i++) {
            for (int j = i + 1; j < fileNames.length; j++) {

                double currentPValue = tTest(fileNames[i], fileNames[j]);
                pValueSum = pValueSum + currentPValue;
                numberOfEntries++;
                if (numberOfEntries % 5000 == 0) {
                    StdOut.println("Number of Entries = " + numberOfEntries + " current sum = "
                                           + pValueSum);
                }
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
    public static double differentDistancesMeanPValues(double distance1, double distance2,
                                                       int config, int orientation) {
        String[] fileNames1 = sameDistanceFiles(distance1, config, orientation);
        String[] fileNames2 = sameDistanceFiles(distance2, config, orientation);
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
                if (numberOfEntries % 5000 == 0) {
                    StdOut.println("Number of Entries = " + numberOfEntries + " current sum = "
                                           + pValueSum);
                }
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

    public static double guessDistance(String fileName) {
        double actualDistance = getDistance(fileName);
        In input = new In("EventsAndDistances.tsv");
        input.readLine();
        int j = 0;
        String currentMaxFile = "";
        double currentMaxPValue = 0;
        while (j < 934 && input.hasNextChar()) {
            String currentFileName = input.readString();
            double fileDistance = input.readDouble();
            if (fileName.equals(currentFileName)) {
                continue;
            }
            double currentPValue = tTest(fileName, currentFileName);
            if (currentPValue > currentMaxPValue) {
                currentMaxPValue = currentPValue;
                currentMaxFile = currentFileName;
            }
            j++;
        }
        input.close();
        double calculatedDistance = getDistance(currentMaxFile);
        StdOut.println("Actual Distance = " + actualDistance + " Calculated Distance = "
                               + calculatedDistance + " Input Filename = " + fileName
                               + " Matched Filename = " + currentMaxFile);
        return calculatedDistance;
    }

    public static int testGuesses() {
        String[] allFiles = new String[934];
        In input = new In("EventsAndDistances.tsv");
        input.readLine();
        int j = 0;
        while (j < 934 && input.hasNextChar()) {
            String currentFileName = input.readString();
            double fileDistance = input.readDouble();
            allFiles[j] = currentFileName;
            j++;
        }
        input.close();
        int correctDistanceGuesses = 0;
        int correctEventGuesses = 0;
        for (int i = 0; i < allFiles.length; i++) {
            double guessedDistance = guessDistance(allFiles[i]);
            double actualDistance = getDistance(allFiles[i]);
            if (actualDistance == guessedDistance) {
                correctDistanceGuesses++;
            }
            if (actualDistance > 6.1 && guessedDistance > 6.1) {
                correctEventGuesses++;
            }
            if (actualDistance < 6.1 && guessedDistance < 6.1) {
                correctEventGuesses++;
            }
        }
        StdOut.println("The number of correct distance guesses = " + correctDistanceGuesses);
        StdOut.println("The number of correct event guesses = " + correctEventGuesses);
        return correctDistanceGuesses;
    }


    public static void main(String[] args) {
        /*double distance = Double.parseDouble(args[0]);
        double distance2 = Double.parseDouble(args[1]);
        int config = Integer.parseInt(args[2]);
        int orientation = Integer.parseInt(args[3]);*/
        /*double averagePValue = meanPValues(distance, config, orientation);
        StdOut.println("Average P-Value for distance " + distance + " = " + averagePValue);
        double averagePValue2 = meanPValues(distance2, config, orientation);
        StdOut.println("Average P-Value for distance " + distance2 + " = " + averagePValue2);*/
        /*double diffAvg = differentDistancesMeanPValues(distance, distance2, config, orientation);
        StdOut.println("Average P-Value for distances " + distance + " and " + distance2 +
                               " = " + diffAvg);*/
        deviceFreq();

        //testGuesses();
    }
}
