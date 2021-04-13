# Distance-Estimation2

dev: Dirctory containing .csv files given in the NIST data set.

dataFiles: Directory containing edited files from dev which contain only Distance in ft, TxDevice, TxPower, RxDevice and RSSI readings.

DataStats.java: Program used for statistically analyzing the files in dataFiles.

EventsAndDistances.tsv: A tab seperated file where on each line we have a file name (of a file in dataFiles) and the distance between the two devices in that file in feet.

FileWriter.java: A small program used to modify the given NIST files into the files in dataFiles.

tc4tl_dev_system_output.tsv: A tab seperated file where on each line we have a file name (of a file in dataFiles) and the distance between the two devices in that file in meters.

Orientation.java: Program used for statistically analyzing the files in dataFiles while considering parameters such as their orientation, devices, TxPower, etc.

Analysis.java: Has almost the same functionality as Orientation.java except Orientation.java is not very efficient since it uses a lot of file reading while Analysis.java uses it only once in the beginning and hence is much more efficient.
