import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.Out;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

public class FileWriter {

    public static boolean shouldCopy(String line, int i) {
        if (i < 3) {
            return true;
        }
        else if (i < 7) {
            return false;
        }
        else {
            return line.contains("Bluetooth");
        }
    }

    public static void main(String[] args) {
        StdIn.readLine();
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            double distance = StdIn.readDouble();
            StdOut.println(word + "\t" + distance);
            In input = new In("dev\\" + word);
            StringBuilder allData = new StringBuilder();
            int i = 0;
            while (input.hasNextLine()) {
                String line = input.readLine();
                if (shouldCopy(line, i)) {
                    if (i > 0) {
                        allData.append("\n");
                    }
                    allData.append(line);
                }
                i++;
            }
            Out output = new Out("dataFiles\\" + word);
            distance = distance / 0.30;
            output.println("Distance," + distance);
            output.print(allData.toString());
            input.close();
            output.close();
        }

    }
}
