import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import static java.lang.System.exit;

public class SigSpot {
    static public String fEncoding = "UTF-8";
    static String path = "./GenGraph/";

    private static class Parameters {
        public static String network_file;
        public static double threshold;
        public static int num_stop;
        public static String output_file;
        // optional
        public static double mu;
        public static NetAlgorithms.Type method;
        public static String statistic_file;
        public static String ground_truth_file;

        private static void print_usage_and_exit() {
            System.out.println("usage: SigSpot <network_file> <threshold> <num_stop> <output_file> " +
                    "[-method hsmss(default)|neigh|con|multcon|multneigh|baseline|maxagg|max|rand|maxed] " +
                    "[-mu <mu>(default 0.01)] [-stat <statistics_file>] [-ground <ground_truth_file>]");
            exit(1);
        }

        public static void parseCommandLine(String[] argv) {
            if (argv.length < 4) {
                print_usage_and_exit();
            }
            int i = 0;

            network_file = argv[i++];
            threshold = Double.parseDouble(argv[i++]);
            num_stop = Integer.parseInt(argv[i++]);
            output_file = argv[i++];

            //optional
            mu = 0.01;
            method = NetAlgorithms.Type.hsmss2;
            statistic_file = null;
            ground_truth_file = null;
            while (i < argv.length) {
                String opt = argv[i++];
                if (i >= argv.length) {
                    print_usage_and_exit();
                }
                if (opt.compareTo("-method") == 0) {
                    String strMethod = argv[i++];
                    if (strMethod.equals("neigh"))
                        method = NetAlgorithms.Type.neigh;
                    else if (strMethod.equals("con"))
                        method = NetAlgorithms.Type.con;
                    else if (strMethod.equals("multcon"))
                        method = NetAlgorithms.Type.multcon;
                    else if (strMethod.equals("multneigh"))
                        method = NetAlgorithms.Type.multneigh;
                    else if (strMethod.equals("hsmss"))
                        method = NetAlgorithms.Type.hsmss;
                    else if (strMethod.equals("baseline"))
                        method = NetAlgorithms.Type.baseline;
                    else if (strMethod.equals("baseline_opt"))
                        method = NetAlgorithms.Type.baseline_opt;
                    else if (strMethod.equals("maxagg"))
                        method = NetAlgorithms.Type.maxagg;
                    else if (strMethod.equals("max"))
                        method = NetAlgorithms.Type.max;
                    else if (strMethod.equals("rand"))
                        method = NetAlgorithms.Type.rand;
                    else if (strMethod.equals("max2"))
                        method = NetAlgorithms.Type.max2;
                    else if (strMethod.equals("multcon2"))
                        method = NetAlgorithms.Type.multcon2;
                    else if (strMethod.equals("multcon3"))
                        method = NetAlgorithms.Type.multcon3;
                    else if (strMethod.equals("mf"))
                        method = NetAlgorithms.Type.mf;
                    else if (strMethod.equals("hsmss2"))
                        method = NetAlgorithms.Type.hsmss2;
                    else if (strMethod.equals("hsmss3"))
                        method = NetAlgorithms.Type.hsmss3;
                    else if (strMethod.equals("maxed"))
                        method = NetAlgorithms.Type.maxed;
                    else {
                        print_usage_and_exit();
                    }
                } else if (opt.compareTo("-ground") == 0) {
                    ground_truth_file = argv[i++];
                } else if (opt.compareTo("-stat") == 0) {
                    statistic_file = argv[i++];
                } else if (opt.compareTo("-mu") == 0) {
                    mu = Double.parseDouble(argv[i++]);
                } else {
                    print_usage_and_exit();
                }
            }
        }
    }

    public static void main(String[] argv) throws Exception {
        //wikipedia990.quadruples 10 10 wikipedia990.out
        Parameters.parseCommandLine(argv);

        System.out.println("Loading " + Parameters.network_file);
        ElapsedTime et = ElapsedTime.Start();
        DynaGraph network = new DynaGraph();
        network.read(Parameters.network_file, Parameters.mu);

        System.out.println("Loading_time: " + ElapsedTime.Stop(et));

        List<Subgraph> patGraphOut = new LinkedList<Subgraph>();
        List<Subsequence> patSeqOut = new LinkedList<Subsequence>();
        ArrayList<Double> times = new ArrayList<Double>();
        ArrayList<Double> scores = new ArrayList<Double>();
        ArrayList<Double> jaccards = new ArrayList<Double>();

        DynaGraph networkClone = network.clone();

        System.out.println("Running...");
        et = ElapsedTime.Start();

        if (Parameters.method == NetAlgorithms.Type.baseline) {
            NetAlgorithms.TopKMeden(network, patGraphOut, patSeqOut, networkClone, (float) Parameters.threshold,
                    false, times, scores, Integer.MAX_VALUE);
        } else if (Parameters.method == NetAlgorithms.Type.baseline_opt) {
            NetAlgorithms.TopKMeden(network, patGraphOut, patSeqOut, networkClone, (float) Parameters.threshold,
                    true, times, scores, Integer.MAX_VALUE);
        } else {
            ComputeJaccard compJaccard = null;
            if (Parameters.ground_truth_file != null) {
                compJaccard = new ComputeJaccard();
                try {
                    compJaccard.loadGroundTruth(Parameters.ground_truth_file);
                } catch (Exception e) {
                    e.printStackTrace();
                    System.out.println("Error in loading the ground truth file. " +
                            "The Jaccard coefficient will not be computed");
                    Parameters.ground_truth_file = null;
                }
            }

            NetAlgorithms.RefineExperiment(network, patGraphOut, patSeqOut,
                    networkClone, new TopDown(), Parameters.num_stop, (float) Parameters.threshold,
                    Integer.MAX_VALUE, Parameters.method, compJaccard, times, scores, jaccards);
        }

        float elapsedTime = ElapsedTime.Stop(et);
        System.out.println("Running time: " + elapsedTime);

        save(Parameters.output_file, network, patGraphOut, patSeqOut, scores);
        if (Parameters.statistic_file != null)
            saveStat(Parameters.statistic_file, Parameters.threshold, patGraphOut, patSeqOut, scores, times, jaccards);
    }

    private static void save(String fileName, DynaGraph network, List<Subgraph> patGraphOut,
                             List<Subsequence> patSeqOut, ArrayList<Double> scores) throws IOException {
        File file = new File(fileName);
        Writer out = new BufferedWriter(new FileWriter(file));

        Iterator<Subgraph> it = patGraphOut.iterator();
        Iterator<Subsequence> itS = patSeqOut.iterator();
        Iterator<Double> itSc = scores.iterator();

        int i = 0;
        while (it.hasNext() && itS.hasNext() && itSc.hasNext()) {
            Subgraph s = it.next();
            Subsequence ss = itS.next();
            Double sc = itSc.next();
            out.write("\n# Region: " + i + " Score: " + sc + "\n");
            for (int j = 0; j < ss.size(); j++) {
                for (Edge e : s)
                    out.write("" + network.MapNames.getName(e.source()) + ", " + network.MapNames.getName(e.dest())
                            + ", " + ss.get(j) + "\n");
            }
            i++;
        }
        out.close();
    }

    private static void saveStat(String fileName, double threshold, List<Subgraph> patGraphOut,
                                 List<Subsequence> patSeqOut, ArrayList<Double> scores, ArrayList<Double> times,
                                 ArrayList<Double> jaccards) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
        double sum = 0.0;
        Iterator<Subsequence> it = patSeqOut.iterator();
        Iterator<Subgraph> it2 = patGraphOut.iterator();
        for (int i = 0; i < times.size(); i++) {
            if (!it.hasNext())
                throw new AssertionError("Unexpected end of list patterns");
            if (!it2.hasNext())
                throw new AssertionError("Unexpected end of list patterns");
            if (scores.get(i) >= threshold) sum += scores.get(i);
            Subsequence seq = it.next();
            Subgraph subgraph = it2.next();
            double jaccard;
            if (jaccards == null || jaccards.size() == 0)
                jaccard = 0.0;
            else
                jaccard = jaccards.get(i);
            bw.write(times.get(i) + "\t" + scores.get(i) + "\t" + sum + "\t" + jaccard + "\t" + seq.getPrintShort()
                    + "\t" + subgraph.getPrintShort() + "\n");
        }
        bw.close();
    }
}
