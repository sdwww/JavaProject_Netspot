import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Scanner;

public class ComputeJaccard {
    HashSet<DynaGraph.EdgeTimetick> ground_truth;
    int count = 0;
    int overlap = 0;


    public void loadGroundTruth(String fFileName) throws IOException {
        ground_truth = new HashSet<DynaGraph.EdgeTimetick>();

        int count = 0;
        int overlap = 0;
        Scanner scanner = new Scanner(new FileInputStream(fFileName), Config.fEncoding);
        try {
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                count++;
                if (line.trim().length() == 0 || line.trim().charAt(0) == '#')
                    continue;
                Scanner scannerLine = new Scanner(line);
                scannerLine.useDelimiter(",");
                int s, d, t;
                float w;
                if (scannerLine.hasNext()) {
                    String token = scannerLine.next().trim();
                    s = new Integer(token);

                    if (!scannerLine.hasNext())
                        throw new IOException("Incomplete row # " + count);
                    token = scannerLine.next().trim();
                    d = new Integer(token);

                    if (!scannerLine.hasNext())
                        throw new IOException("Incomplete row # " + count);
                    token = scannerLine.next().trim();
                    t = new Integer(token);

                    ground_truth.add(new DynaGraph.EdgeTimetick(new Edge(s, d), t));
                }
            }
        } finally {
            scanner.close();
        }
    }

    public void add(Subgraph patGraphNew, Subsequence patSeqNew) {
        for (int t : patSeqNew) {
            for (Edge e : patGraphNew) {
                if (ground_truth.contains(new DynaGraph.EdgeTimetick(e, t)))
                    overlap++;   // TODO: there is an error here: save in a hashSet to manage overlap
                count++;
            }
        }
    }

    public double getJaccard() {
        if (count + ground_truth.size() - overlap == 0) {
            assert overlap == 0;
            return 0;
        }
        return overlap / (count + ground_truth.size() - overlap);  // TODO: something is not working here
    }
}
