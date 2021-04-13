import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;
import edu.princeton.cs.algs4.StdStats;
import edu.princeton.cs.algs4.Stopwatch;
import org.apache.commons.math3.stat.inference.TTest;

import java.util.ArrayList;
import java.util.HashMap;

public class Analysis {

    private class FileData {
        private String fileName;
        private double distance;
        private int configuration;
        private int orientation;
        private int[] rssiData;

        public FileData(String fileName) {
            this.fileName = fileName;
            this.distance = OrientationDataStats.getDistance(fileName);
            int[] data = OrientationDataStats.getData(fileName);
            this.configuration = data[1];
            this.orientation = data[2];
            this.rssiData = OrientationDataStats.getRSSIArray(fileName);
        }

    }

    private HashMap<String, FileData> fileMap;

    /**
     * Opens the file containing the names of all files in the data set and constructs a set of
     * TC4TL events from each file in the data set.
     *
     * @param setOfFiles is the file containing the names of all the files in the data set.
     */
    public Analysis(String setOfFiles) {
        fileMap = new HashMap<String, FileData>();
        In input = new In(setOfFiles);
        input.readLine();
        int i = 0;
        while (!input.isEmpty()) {
            String fileName = input.readString();
            input.readDouble();
            FileData currentFile = new FileData(fileName);
            fileMap.put(fileName, currentFile);
            /*if (i % 500 == 0) {
                StdOut.println("Constructing file number " + i + " file name " + fileName);
            }
            i++;*/
        }
        input.close();
    }

    public Analysis(Analysis fullSet, String[] fileNames) {
        fileMap = new HashMap<String, FileData>();
        for (int i = 0; i < fileNames.length; i++) {
            FileData currentData = fullSet.fileMap.get(fileNames[i]);
            fileMap.put(fileNames[i], currentData);
        }
    }

    /**
     * Finds the mean rssi value in a given event.
     *
     * @param fileName The name of the file containing the data of the event.
     * @return The mean rssi value in the given event.
     */
    public double mean(String fileName) {
        int[] rssi = fileMap.get(fileName).rssiData;
        return StdStats.mean(rssi);
    }

    public double stddev(String fileName) {
        int[] rssi = fileMap.get(fileName).rssiData;
        return StdStats.stddev(rssi);
    }

    /**
     * Checks if the event whose data is stored in the file @param fileName falls under the category
     * mentioned using integers @param config (the Tx and Rx devices and TxPower) and @param
     * orientation (Whether the phone is handheld or in pocket). If config < 0, we ignore the device
     * details of this event. Similarly if orientation < 0, we ignore the orientation details of
     * this event.
     *
     * @param fileName    The name of the file containing the data of the given event.
     * @param config      Encodes The Tx and Rx devices and Tx power
     * @param orientation Encodes whether the devices are handheld or in pocket.
     * @return True if the given event has the same config (if config >= 0) and same orientation if
     * (orientation >= 0), false otherwise.
     */
    public boolean isSelected(String fileName, int config, int orientation) {
        boolean selection = true;
        int configurations = fileMap.get(fileName).configuration;
        if (config >= 0) {
            selection = (config == configurations);
        }
        int currentOrientation = fileMap.get(fileName).orientation;
        if (orientation >= 0) {
            selection = selection && (orientation == currentOrientation);
        }
        return selection;
    }

    /**
     * A program which parses through all the events which exhibit the given configuration and
     * orientation and calculates each of their mean rssi values. It returns the highest mean rssi
     * value and prints the details of the event which exhibits the highest mean rssi values.
     *
     * @param config      The given configuration (Tx and Rx devices and TxPower).
     * @param orientation The given orientation (Whether the phones are handheld or in pocket).
     * @return The highest mean rssi value among the mean rssi values of events which exhibit the
     * given configuration and orientation.
     */
    public double highestMean(int config, int orientation) {
        double currentHighest = Double.NEGATIVE_INFINITY;
        double currentDistance = 0;
        String highestFile = "";
        for (String key : fileMap.keySet()) {
            if (isSelected(key, config, orientation)) {
                FileData currentFile = fileMap.get(key);
                double distance = currentFile.distance;
                double currentMean = mean(key);
                if (currentMean > currentHighest) {
                    currentHighest = currentMean;
                    currentDistance = distance;
                    highestFile = key;
                }
            }
        }
        StdOut.println("Highest Mean = " + currentHighest + " Distance = " + currentDistance
                               + " file name = " + highestFile);
        return currentHighest;
    }

    /**
     * This function takes 2 events whose data is stored in files 'file1' and 'file2' and performs a
     * T-Test between the set of their respective rssi values. It returns the probability that both
     * these sets belong to the same distribution.
     *
     * @param file1 The name of the file containing the data of the first event.
     * @param file2 The name of the file containing the data of the second event.
     * @return The probability that both the set of rssi values belong to the same distribution.
     */
    public double Ttest(String file1, String file2) {
        int[] rssi1 = fileMap.get(file1).rssiData;
        int[] rssi2 = fileMap.get(file2).rssiData;
        double[] set1 = new double[rssi1.length];
        double[] set2 = new double[rssi2.length];
        for (int i = 0; i < set1.length; i++) {
            set1[i] = rssi1[i];
        }

        for (int i = 0; i < rssi2.length; i++) {
            set2[i] = rssi2[i];
        }
        TTest tTest = new TTest();
        double pValue = tTest.tTest(set1, set2);
        return pValue;
    }

    /**
     * This function takes as input another data set (test), an event from another data set (stored
     * in file named file1) and an event from the current data set (stored in file2). It collects
     * the set of rssi values for both of these events and performs a T-Test between these sets. It
     * returns the probability that these sets belong to the same distribution.
     *
     * @param test  The other data set.
     * @param file1 Name of the file containing the data of an event from the data set 'test'.
     * @param file2 Name of the file containing the data of an event from the current data set.
     * @return The probability that the sets of rssi values from 'file1' and 'file2' both belong to
     * the same distribution.
     */
    public double Ttest(Analysis test, String file1, String file2) {
        int[] rssi1 = test.fileMap.get(file1).rssiData;
        int[] rssi2 = fileMap.get(file2).rssiData;
        double[] set1 = new double[rssi1.length];
        double[] set2 = new double[rssi2.length];
        for (int i = 0; i < set1.length; i++) {
            set1[i] = rssi1[i];
        }

        for (int i = 0; i < rssi2.length; i++) {
            set2[i] = rssi2[i];
        }
        TTest tTest = new TTest();
        double pValue = tTest.tTest(set1, set2);
        return pValue;
    }

    /**
     * This function returns the set of files which encode events where the distance between the 2
     * devices is @param distance, the configuration is @param config and the orientation of the
     * devices is @param orientation.
     *
     * @param distance    The desired distance between the devices.
     * @param config      The desired configuration (Tx and Rx device and TxPower).
     * @param orientation The desired orientation (Whether phones are handheld or in pocket).
     * @return The names of all the files which exhibit the given 3 properties(@param config, @param
     * distance and @param orientation).
     */
    public String[] sameDistanceFiles(double distance, int config, int orientation) {
        ArrayList<String> fileList = new ArrayList<String>();
        for (String key : fileMap.keySet()) {
            FileData currentFile = fileMap.get(key);
            double currentDistance = currentFile.distance;
            if ((currentDistance == distance) && isSelected(key, config, orientation)) {
                fileList.add(key);
            }
        }
        String[] files = new String[fileList.size()];
        for (int i = 0; i < files.length; i++) {
            files[i] = fileList.get(i);
        }
        return files;
    }

    /**
     * This function groups all the events which exhibit the 3 given properties i.e. the distance,
     * configuration and orientation and stores their corresponding file names in an array. It now
     * takes each unordered pair of events from this group and performs a T-Test and computes the
     * P-Value (Probability that the rssi values in both the events belong to the same set). After
     * computing n choose 2 such P-Values (where n = number of events in the group), it computes the
     * mean P-Values and returns it as well as prints some statistical data about the group.
     *
     * @param distance    The desired distance between the devices in the group of events.
     * @param config      The configuration (Tx and Rx device and TxPower) of devices in the group
     *                    of events.
     * @param orientation The orientation of devices in the group of events (handheld or in
     *                    pocket).
     * @return Takes each unordered pair of events and performs a T-Test between their rssi value
     * sets. Computes the probability that these 2 sets belong to the same distribution for C(n,2)
     * such pairs and returns the mean of these C(n,2) P-Values.
     */
    public double meanPValues(double distance, int config, int orientation) {
        String[] fileNames = sameDistanceFiles(distance, config, orientation);
        double pValueSum = 0;
        double numberOfEntries = 0;
        double maxPValue = 0;
        for (int i = 0; i < fileNames.length; i++) {
            for (int j = i + 1; j < fileNames.length; j++) {

                double currentPValue = Ttest(fileNames[i], fileNames[j]);
                pValueSum = pValueSum + currentPValue;
                numberOfEntries++;
                /*if (numberOfEntries % 50000 == 0) {
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
        StdOut.println("Mean P-Value for distance " + distance + "ft = " + meanPvalue);
        return meanPvalue;
    }

    public double differentDistancesMeanPValue(double distance1, double distance2, int config,
                                               int orientation) {
        String[] fileNames1 = sameDistanceFiles(distance1, config, orientation);
        String[] fileNames2 = sameDistanceFiles(distance2, config, orientation);
        double pValueSum = 0;
        double numberOfEntries = 0;
        double maxPValue = 0;
        String currentMaxFile1 = "";
        String currentMaxFile2 = "";
        for (int i = 0; i < fileNames1.length; i++) {
            for (int j = 0; j < fileNames2.length; j++) {
                double currentPValue = Ttest(fileNames1[i], fileNames2[j]);
                pValueSum = pValueSum + currentPValue;
                numberOfEntries++;
                /*if (numberOfEntries % 50000 == 0) {
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
        StdOut.println("Mean P-Value for distances " + distance1 + " and " + distance2 + " = "
                               + meanPvalue);
        return meanPvalue;
    }

    public double meanPValueFile(String fileName, String[] files) {
        double currentSum = 0;
        for (int i = 0; i < files.length; i++) {
            currentSum = currentSum + Ttest(fileName, files[i]);
        }
        double meanPValue = currentSum / files.length;
        return meanPValue;
    }

    public double meanPValueFile(Analysis test, String fileName, String[] files) {
        double currentSum = 0;
        for (int i = 0; i < files.length; i++) {
            currentSum = currentSum + Ttest(test, fileName, files[i]);
        }
        double meanPValue = currentSum / files.length;
        return meanPValue;
    }

    public double guessDistance(String fileName) {
        String currentMaxFile = "";
        double currentMaxPValue = 0;
        FileData givenFile = fileMap.get(fileName);
        int givenOrientation = givenFile.orientation;
        int givenConfiguration = givenFile.configuration;
        for (String key : fileMap.keySet()) {
            if (key.equals(fileName)) {
                continue;
            }
            if (!isSelected(key, givenConfiguration, givenOrientation)) {
                continue;
            }
            double currentPValue = Ttest(fileName, key);
            if (currentPValue > currentMaxPValue) {
                currentMaxFile = key;
                currentMaxPValue = currentPValue;
            }
        }
        double calculatedDistance = fileMap.get(currentMaxFile).distance;
        double actualDistance = givenFile.distance;
        StdOut.println("Actual Distance = " + actualDistance + " Calculated Distance = "
                               + calculatedDistance + " Input Filename = " + fileName
                               + " Matched Filename = " + currentMaxFile);
        return calculatedDistance;
    }

    public double guessDistanceApproach2(String fileName) {
        FileData givenFile = fileMap.get(fileName);
        int givenOrientation = givenFile.orientation;
        int givenConfiguration = givenFile.configuration;
        String[] dist4 = sameDistanceFiles(4, givenConfiguration, givenOrientation);
        String[] dist6 = sameDistanceFiles(6, givenConfiguration, givenOrientation);
        String[] dist10 = sameDistanceFiles(10, givenConfiguration, givenOrientation);
        String[] dist15 = sameDistanceFiles(15, givenConfiguration, givenOrientation);
        double[] pValueMeans = new double[4];
        pValueMeans[0] = meanPValueFile(fileName, dist4);
        pValueMeans[1] = meanPValueFile(fileName, dist6);
        pValueMeans[2] = meanPValueFile(fileName, dist10);
        pValueMeans[3] = meanPValueFile(fileName, dist15);
        double[] distances = { 4, 6, 10, 15 };
        double currentMaxMean = pValueMeans[0];
        double distance = 4;
        for (int i = 1; i < 4; i++) {
            if (pValueMeans[i] > currentMaxMean) {
                currentMaxMean = pValueMeans[i];
                distance = distances[i];
            }
        }
        return distance;
    }

    public int testGuess() {
        int correctDistanceGuesses = 0;
        int correctEventGuesses = 0;
        int falsePositives = 0;
        int falseNegatives = 0;
        HashMap<String, Integer> descriptions = new HashMap<String, Integer>();
        int i = 0;
        for (String key : fileMap.keySet()) {
            double guessedDistance = guessDistanceApproach2(key);
            FileData currentFile = fileMap.get(key);
            double actualDistance = currentFile.distance;
            if (actualDistance == guessedDistance) {
                correctDistanceGuesses++;
            }
            if (actualDistance > 6.1 && guessedDistance > 6.1) {
                correctEventGuesses++;
            }
            if (actualDistance < 6.1 && guessedDistance < 6.1) {
                correctEventGuesses++;
            }
            if (actualDistance < 6.1 && guessedDistance > 6.1) {
                String currentDes = "" + currentFile.configuration + currentFile.orientation;
                if (descriptions.containsKey(currentDes)) {
                    int newValue = descriptions.get(currentDes) + 1;
                    descriptions.put(currentDes, newValue);
                }
                else {
                    descriptions.put(currentDes, 1);
                }
                falseNegatives++;
            }
            if (actualDistance > 6.1 && guessedDistance < 6.1) {
                String currentDes = "" + currentFile.configuration + currentFile.orientation;
                if (descriptions.containsKey(currentDes)) {
                    int newValue = descriptions.get(currentDes) + 1;
                    descriptions.put(currentDes, newValue);
                }
                else {
                    descriptions.put(currentDes, 1);
                }
                falsePositives++;
            }
            /*if (i % 500 == 0) {
                StdOut.println("Estimating distance of file " + key + " which is filenumber " + i);
            }
            i++;*/

        }
        StdOut.println("The number of correct distance guesses = " + correctDistanceGuesses);
        StdOut.println("The number of correct event guesses = " + correctEventGuesses);
        StdOut.println("The number of False Positives = " + falsePositives);
        StdOut.println("The number of False Negatives = " + falseNegatives);
        for (String key : descriptions.keySet()) {
            String config = "" + key.charAt(0);
            String orient = "" + key.charAt(1);
            int configuration = Integer.parseInt(config);
            int orientation = Integer.parseInt(orient);
            config = OrientationDataStats.CONFIGURATIONS[configuration];
            orient = OrientationDataStats.ORIENTATIONS[orientation];
            String eventDescription = config + orient;
            //StdOut.println(eventDescription + "\t" + descriptions.get(key));
        }
        return correctDistanceGuesses;
    }

    public double guessDistance(String fileName, Analysis test) {
        FileData givenFile = test.fileMap.get(fileName);
        int givenOrientation = givenFile.orientation;
        int givenConfiguration = givenFile.configuration;

        String[] dist4 = sameDistanceFiles(4, givenConfiguration, givenOrientation);
        String[] dist6 = sameDistanceFiles(6, givenConfiguration, givenOrientation);
        String[] dist10 = sameDistanceFiles(10, givenConfiguration, givenOrientation);
        String[] dist15 = sameDistanceFiles(15, givenConfiguration, givenOrientation);

        int totalNumberOfFiles = dist4.length + dist6.length + dist10.length + dist15.length;
        //StdOut.println(totalNumberOfFiles);

        double[] pValueMeans = new double[4];
        pValueMeans[0] = meanPValueFile(test, fileName, dist4);
        pValueMeans[1] = meanPValueFile(test, fileName, dist6);
        pValueMeans[2] = meanPValueFile(test, fileName, dist10);
        pValueMeans[3] = meanPValueFile(test, fileName, dist15);
        double[] distances = { 4, 6, 10, 15 };
        double currentMaxMean = pValueMeans[0];
        double distance = 4;
        for (int i = 1; i < 4; i++) {
            if (pValueMeans[i] > currentMaxMean) {
                currentMaxMean = pValueMeans[i];
                distance = distances[i];
            }
        }
        return distance;
    }

    public int testGuess(Analysis test) {
        int correctDistanceGuesses = 0;
        int correctEventGuesses = 0;
        int falsePositives = 0;
        int falseNegatives = 0;
        //HashMap<String, Integer> descriptions = new HashMap<String, Integer>();
        for (String key : test.fileMap.keySet()) {
            double guessedDistance = guessDistance(key, test);
            FileData currentFile = test.fileMap.get(key);
            double actualDistance = currentFile.distance;
            if (actualDistance == guessedDistance) {
                correctDistanceGuesses++;
            }
            if (actualDistance > 6.1 && guessedDistance > 6.1) {
                correctEventGuesses++;
            }
            if (actualDistance < 6.1 && guessedDistance < 6.1) {
                correctEventGuesses++;
            }
            if (actualDistance < 6.1 && guessedDistance > 6.1) {
                falseNegatives++;
            }
            if (actualDistance > 6.1 && guessedDistance < 6.1) {
                falsePositives++;
            }
        }
        //StdOut.println("The number of correct distance guesses = " + correctDistanceGuesses);
        //StdOut.println("The number of correct event guesses = " + correctEventGuesses);
        //StdOut.println("The number of False Positives = " + falsePositives);
        //StdOut.println("The number of False Negatives = " + falseNegatives);
        /*for (String key : descriptions.keySet()) {
            String config = "" + key.charAt(0);
            String orient = "" + key.charAt(1);
            int configuration = Integer.parseInt(config);
            int orientation = Integer.parseInt(orient);
            config = OrientationDataStats.CONFIGURATIONS[configuration];
            orient = OrientationDataStats.ORIENTATIONS[orientation];
            String eventDescription = config + orient;
            //StdOut.println(eventDescription + "\t" + descriptions.get(key));
        }*/
        return correctEventGuesses;
    }

    public static int randomTestGuessWithin(Analysis fullSet) {
        String[] fileNamesArray = new String[fullSet.fileMap.size()];
        int j = 0;
        for (String key : fullSet.fileMap.keySet()) {
            fileNamesArray[j] = key;
            j++;
        }
        StdRandom.shuffle(fileNamesArray);
        int large = (8 * fileNamesArray.length) / 10;
        int small = fileNamesArray.length - large;
        String[] trainSet = new String[large];
        String[] testSet = new String[small];
        for (int i = 0; i < fileNamesArray.length; i++) {
            if (i < large) {
                trainSet[i] = fileNamesArray[i];
            }
            else {
                testSet[i - large] = fileNamesArray[i];
            }
        }
        Analysis train = new Analysis(fullSet, trainSet);
        Analysis test = new Analysis(fullSet, testSet);
        return train.testGuess(test);
    }

    public static int randomTestGuessSecond(Analysis trainSuperSet, Analysis testSuperSet) {
        String[] fileNamesArray = new String[testSuperSet.fileMap.size()];
        int j = 0;
        for (String key : testSuperSet.fileMap.keySet()) {
            fileNamesArray[j] = key;
            j++;
        }
        StdRandom.shuffle(fileNamesArray);
        int large = (8 * fileNamesArray.length) / 10;
        int small = fileNamesArray.length - large;
        //StdOut.println("Number of files = " + small);
        String[] testSet = new String[small];
        for (int i = large; i < fileNamesArray.length; i++) {
            testSet[i - large] = fileNamesArray[i];
        }
        Analysis test = new Analysis(testSuperSet, testSet);
        String[] trainFileNames = new String[trainSuperSet.fileMap.size() - small];
        j = 0;
        for (String key : trainSuperSet.fileMap.keySet()) {
            boolean contains = false;
            for (int k = 0; k < testSet.length; k++) {
                if (testSet[k].equals(key)) {
                    contains = true;
                    break;
                }
            }
            if (j >= trainFileNames.length) {
                break;
            }
            if (!contains) {
                trainFileNames[j] = key;
                j++;
            }
        }
        Analysis train = new Analysis(trainSuperSet, trainFileNames);
        return train.testGuess(test);
    }

    public static void main(String[] args) {
        /*Analysis train = new Analysis("trainEventsData.tsv");
        StdOut.println();
        Analysis test = new Analysis("EventsAndDistances.tsv");
        StdOut.println();
        train.testGuess(test);
        StdOut.println();
        test.testGuess();
        test.testGuess(train);
        StdOut.println();
        analysis.highestMean(-1, -1);
        StdOut.println();
        analysis.meanPValues(6, -1, -1);
        StdOut.println();
        analysis.differentDistancesMeanPValue(4, 6, -1, -1);
        StdOut.println();*/

        Analysis fullSet = new Analysis("trainEventsData.tsv");
        Stopwatch stopwatch = new Stopwatch();
        double sum = 0;
        int size = Math.round((fullSet.fileMap.size() * 2) / 10);
        int trials = Integer.parseInt(args[0]);
        StdOut.println("Total number of events = " + size);
        for (int i = 0; i < trials; i++) {
            sum = sum + randomTestGuessWithin(fullSet);
            // StdOut.println("Test number " + i + " time spent = " + stopwatch.elapsedTime());
            // StdOut.println();
        }
        double meanCorrectGuesses = sum / trials;
        StdOut.println("Mean correct Event guesses = " + meanCorrectGuesses);
        double percentageCorrect = (meanCorrectGuesses / size) * 100;
        StdOut.println("Mean correct Event guesses = " + meanCorrectGuesses);
        StdOut.println(
                "When a single dataset is segregated such that testing is 20% and training is 80%, accuracy = "
                        + percentageCorrect + "%");
        /*Analysis trainSuperSet = new Analysis("trainEventsData.tsv");
        Analysis testSuperSet = new Analysis("simpleDataSet.tsv");
        int size = Math.round((trainSuperSet.fileMap.size() * 2) / 10);
        Stopwatch stopwatch = new Stopwatch();
        double sum = 0;
        int trials = Integer.parseInt(args[0]);
        StdOut.println("Total number of events = " + size);
        for (int i = 0; i < trials; i++) {
            sum = sum + randomTestGuessSecond(testSuperSet, trainSuperSet);
            // StdOut.println("Test number " + i + " time spent = " + stopwatch.elapsedTime());
            // StdOut.println();
        }
        double meanCorrectGuesses = sum / trials;
        double percentageCorrect = (meanCorrectGuesses / size) * 100;
        StdOut.println("Mean correct Event guesses = " + meanCorrectGuesses);
        StdOut.println(
                "When the testing set is more general and the training set is simple, accuracy = "
                        + percentageCorrect);*/
        //trainSuperSet.testGuess(testSuperSet);

    }

}
