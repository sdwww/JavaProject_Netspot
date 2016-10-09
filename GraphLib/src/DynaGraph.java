import java.io.*;
import java.util.*;


public class DynaGraph extends Vector<Graph> {

    private HashMap<Edge, Sequence> allFastSubsequences = null;

    public class _MapNames {
        private HashMap<String, Integer> mapName2Index = new HashMap<String, Integer>();
        private HashMap<Integer, String> mapIndex2Name = new HashMap<Integer, String>();

        int count = 0;

        public String getName(int index) {
            String name = mapIndex2Name.get(index);
            assert name != null;
            return name;
        }

        public int getIndex(String name) {
            Integer i = mapName2Index.get(name);
            if (i == null) {
                mapName2Index.put(name, count);
                mapIndex2Name.put(count, name);
                i = count;
                count++;
            }
            return i;
        }
    }

    _MapNames  MapNames = new _MapNames();

    HashMap<Edge, ArrayList<Edge>> edgeLists = null;

    public void print(Subgraph patGraph, Subsequence patSeq) {
        System.out.print("Edges\t");
        for (Edge e : patGraph) {
            System.out.print("\t(" + MapNames.getName(e.source()) + "," + MapNames.getName(e.dest()) + ")");
        }
        System.out.println("");
        System.out.print("Scores\t");
        for (Edge e : patGraph) {
            Subgraph sg = new Subgraph();
            sg.add(e);
            System.out.print("\t" + GetScore(sg, patSeq));
        }
        System.out.println("");
        System.out.print("MSS\t");
        for (Edge e : patGraph) {
            Subgraph sg = new Subgraph();
            sg.add(e);
            float max = Float.NEGATIVE_INFINITY;
            for (int t : patSeq) {
                Subsequence ss = GetSubsequence(e).MaximumScoreSubsequence(t);
                float sc = GetScore(sg, ss);
                if (sc > max)
                    max = sc;
            }
            System.out.print("\t" + max);
        }
        System.out.println("");
        for (int t : patSeq) {
            Subsequence seq = new Subsequence();
            seq.add(t);
            System.out.print("" + t + "\t" + GetScore(patGraph, seq));
            for (Edge e : patGraph) {
                System.out.print("\t" + this.get(t).get(e));
            }
            System.out.println("");
        }
    }

    private Sequence GetSubsequence(Edge e) {
        Sequence seq = new Sequence();
        seq.setSize(size());
        for (int i = 0; i < size(); i++) {
            Graph g = get(i);
            seq.set(i, g.get(e));
        }

        return seq;
    }

    class Entry {
        public Entry(int ls, int ld, int lt, double lw) {
            s = ls;
            d = ld;
            t = lt;
            w = lw;
        }

        public int s, d, t;
        public double w;
    }

    // usefull if the graph structure is unchanging, so there's no point in recalculating it for each timeslice
    public void getEdgeLists() {
        edgeLists = new HashMap<Edge, ArrayList<Edge>>();
        Graph g = get(0);
        for (Edge e : g.keySet()) {
            ArrayList<Edge> cur = new ArrayList<Edge>();
            for (Edge n : g.keySet()) {
                if (n.dest() == e.source() || n.source() == e.dest() || n.source() == e.source() || n.dest() == e.dest())
                    cur.add(n);
            }
            edgeLists.put(e, cur);
        }

    }

    public float GetScore(Subgraph subgraph, Subsequence subsequence) {
        float score = 0;
        for (int t : subsequence) {
            for (Edge e : subgraph) {
                score = score + this.get(t).get(e);
            }
        }
        return score;
    }

    public float GetEnergy(Subgraph subgraph, Subsequence subsequence) {
        float score = 0;
        for (int t : subsequence) {
            for (Edge e : subgraph) {
                float w = this.get(t).get(e);
                if (w > 0)
                    score = score + w;
            }
        }
        return score;
    }

    public Graph GetInducedSubgraph(Subsequence subsequence) {
        Graph g = new Graph();

        for (int t : subsequence) {
            g.add(this.get(t));
        }

        return g;
    }

    public Sequence GetInducedSubsequence(Subgraph subgraph) {
        Sequence seq = new Sequence();
        seq.setSize(size());
        for (int i = 0; i < size(); i++) {
            Graph g = get(i);
            Graph sg = g.GetSubgraph(subgraph);
            seq.set(i, sg.GetTotalWeight());
        }

        return seq;
    }

    public Sequence GetFastInducedSubsequence(Subgraph subgraph) {
        Sequence seq = new Sequence();
        seq.setSize(size());
        for (int i = 0; i < size(); i++) {
            Graph g = get(i);
            seq.set(i, g.GetSubgraphWeight(subgraph));
        }

        return seq;
    }

    public HashMap<Edge, Sequence> getAllSubsequences() {
        HashMap<Edge, Sequence> res = new HashMap<Edge, Sequence>();
        // initialize
        for (Edge e : get(0).keySet()) {
            res.put(e, new Sequence());
        }

        for (int i = 0; i < this.size(); i++) {
            for (Edge e : this.get(i).keySet()) {
                res.get(e).add(this.get(i).get(e));
            }
        }

        return res;
    }

    public void updateAllFastSubsequences(Edge ed, int t) {
        if (allFastSubsequences != null) {
            allFastSubsequences.get(ed).setElementAt(0f, t);
        }
    }

    public HashMap<Edge, Sequence> getAllFastSubsequences() {
        if (allFastSubsequences == null)
            allFastSubsequences = getAllSubsequences();

        return allFastSubsequences;
    }


    // TODO(razvan,misael): this is whaere the seed generator and refinement need to hook
    // returns the hashmap from edge to the
    // scores of the connected positive component including the edge
    // for every time. the vallue is an array of the size of the times
    public HashMap<Edge, float[]> GetMatrixPositiveConnectedComponent() {
        HashMap<Edge, float[]> res = new HashMap<Edge, float[]>();
        // initialize
        for (Edge e : get(0).keySet()) {
            res.put(e, new float[size()]);
        }
        // fill
        HashMap<Edge, Float> scores = new HashMap<Edge, Float>();
        for (int t = 0; t < size(); t++) {
            //System.out.println("components at: " + t);
            get(t).getEdgePositiveComponentScores(scores);

            for (Edge e : get(t).keySet()) {
                //System.out.println(t + ": " + MapNames.getName(e.source()) + "->" + MapNames.getName(e.dest()) + " = " + scores.get(e) + " ");
                res.get(e)[t] = scores.get(e);
            }
        }
        return res;
    }

    public HashMap<Edge, float[]> GetMatrixPositiveNeighbourScores(HashMap<Integer, ArrayList<Edge>> changes) {

        if (edgeLists == null)
            getEdgeLists();

        HashMap<Edge, float[]> res = new HashMap<Edge, float[]>();
        // initialize
        for (Edge e : get(0).keySet()) {
            res.put(e, new float[size()]);
        }
        // fill
        HashMap<Edge, Float> scores = new HashMap<Edge, Float>();
        for (int t = 0; t < size(); t++) {
            if (changes.containsKey(t))
                scores = get(t).getEdgePositiveNeighbourScores(changes.get(t), edgeLists);
            else
                scores = get(t).getEdgePositiveNeighbourScores(null, edgeLists);
            for (Edge e : get(t).keySet()) {
                res.get(e)[t] = scores.get(e);
            }
        }
        return res;
    }

    public HashMap<Edge, float[]> GetMatrixPositiveNeighbourScores() {
        return GetMatrixPositiveNeighbourScores(new HashMap<Integer, ArrayList<Edge>>());
    }

    public HashMap<Integer, ArrayList<Edge>> Erase(Subgraph patGraph, Subsequence patSeq) {
        HashMap<Integer, ArrayList<Edge>> rez = new HashMap<Integer, ArrayList<Edge>>();

        for (Iterator<Integer> itT = patSeq.iterator(); itT.hasNext(); ) {
            int t = itT.next();
            rez.put(t, new ArrayList<Edge>());
        }

        for (Iterator<Edge> itEd = patGraph.iterator(); itEd.hasNext(); ) {
            Edge ed = itEd.next();
            for (Iterator<Integer> itT = patSeq.iterator(); itT.hasNext(); ) {
                int t = itT.next();
                updateAllFastSubsequences(ed, t);
                Graph cur = get(t);
                if (cur.get(ed) > 0) {
                    cur.put(ed, 0f);
                    rez.get(t).add(ed);
                }
            }
        }


        return rez;
    }

    @Override
    public DynaGraph clone() {
        DynaGraph net = new DynaGraph();
        net.setSize(this.size());
        for (int i = 0; i < this.size(); i++) {
            net.set(i, this.get(i).clone());
        }
        return net;
    }

    public void write(String fFileName) throws IOException {
        double miu = 0.01;
        write(fFileName, miu);
    }

    public void write(String fFileName, double miu) throws IOException {
        Writer out = new OutputStreamWriter(new FileOutputStream(fFileName), Config.fEncoding);
        try {
            for (int i = 0; i < size(); i++) {
                Graph g = this.get(i);
                if (g == null) {
                    continue;
                }

                for (Edge e : g.keySet()) {
                    double w = miu / Math.pow(2, g.get(e));   // inverse of score conversion in read
                    out.write(MapNames.getName(e.source()) + "," + MapNames.getName(e.dest()) + "," + String.valueOf(i) + "," + String.valueOf(w) + "\n");
                }
            }
        } finally {
            out.close();
        }
    }

    /**
     * Read the contents of the given file.
     */
    public void read(String fFileName) throws IOException {
        double miu = 0.01;
        read(fFileName, miu);
    }

    public void read(String fFileName, double miu) throws IOException {
        Scanner scanner = new Scanner(new FileInputStream(fFileName), Config.fEncoding);
        LinkedList<Entry> l = new LinkedList<Entry>();
        int node_count = -1;
        int time_count = -1;
        try {
            while (scanner.hasNextLine()) {
                Entry e = processLine(scanner.nextLine());
                if (e == null)
                    continue;
                l.add(e);

                node_count = Math.max(node_count, Math.max(e.s, e.d));
                time_count = Math.max(time_count, e.t);
            }
            node_count++;
            time_count++;

            this.setSize(time_count);
            for (Entry e : l) {
                Graph g = this.get(e.t);
                if (g == null) {
                    g = new Graph();
                    this.set(e.t, g);
                }

                if (e.s > e.d) {
                    int dest = e.s;
                    e.s = e.d;
                    e.d = dest;
                } //swap
                float score;
                if (miu != 0)
                    score = (float) (-Math.log(e.w / miu) / Math.log(2));
                else
                    score = (float) e.w;

                g.put(new Edge(e.s, e.d), score);
            }
        } finally {
            scanner.close();
        }
    }

    protected Entry processLine(String aLine) throws IOException {
        //use a second Scanner to parse the content of each line
        Scanner scanner = new Scanner(aLine);
        scanner.useDelimiter(",");
        int s, d, t;
        float w;
        if (scanner.hasNext()) {
            String token = scanner.next().trim();
            s = MapNames.getIndex(token);

            if (!scanner.hasNext())
                throw new IOException("Incomplete row");
            token = scanner.next().trim();
            d = MapNames.getIndex(token);

            if (!scanner.hasNext())
                throw new IOException("Incomplete row");
            token = scanner.next().trim();
            t = new Integer(token);

            if (!scanner.hasNext())
                throw new IOException("Incomplete row");
            token = scanner.next().trim();
            w = new Float(token);
            return new Entry(s, d, t, w);
        }
        return null;
    }

    // this is a unit test for positive component scores
    private static boolean testPosCompMat() {
        // create a graph
        HashMap<Edge, Float> actpcs = new HashMap<Edge, Float>();
        Graph g = new Graph();
        Edge e;
        e = new Edge(0, 1);
        g.put(e, (float) 2.0);
        actpcs.put(e, (float) 2);
        e = new Edge(1, 2);
        g.put(e, (float) -1.0);
        actpcs.put(e, (float) -1);
        e = new Edge(1, 3);
        g.put(e, (float) -1.0);
        actpcs.put(e, (float) -1);
        e = new Edge(1, 4);
        g.put(e, (float) -8.0);
        actpcs.put(e, (float) -8);
        e = new Edge(2, 3);
        g.put(e, (float) -2.0);
        actpcs.put(e, (float) -2);
        e = new Edge(3, 4);
        g.put(e, (float) 2.0);
        actpcs.put(e, (float) 6);
        e = new Edge(3, 5);
        g.put(e, (float) 1.0);
        actpcs.put(e, (float) 6);
        e = new Edge(4, 5);
        g.put(e, (float) 3.0);
        actpcs.put(e, (float) 6);

        DynaGraph dg = new DynaGraph();
        dg.add(g);

        HashMap<Edge, float[]> mat = dg.GetMatrixPositiveConnectedComponent();

        // test if actpcs is the same as poscompres
        assert (mat.size() == actpcs.size()) : "Wrong number of positive components";
        for (Edge ed : g.keySet()) {
            assert (mat.containsKey(ed)) :
                    "Edge " + ed.toString() + " missing from positive components\n";
            assert (mat.get(ed).length == 1) : "wrong number of time stamps";
            assert (((Float) mat.get(ed)[0]).equals(actpcs.get(ed))) :
                    "Wrong pos comp score for edge " + ed.toString() + "= " + mat.get(ed)[0] +
                            ". Should be " + actpcs.get(ed) + " instead.\n";
        }

        return true;
    }

    public static void main(String[] args) throws IOException {
        if (testPosCompMat()) {
            System.out.print("Tested GraphLib.DynaGraph successfully.\n");
        }
        ;
    }

    public static class EdgeTimetick {
        EdgeTimetick(Edge _ed, int _t) {
            ed = _ed;
            t = _t;
        }

        @Override
        public boolean equals(Object obj) {
            EdgeTimetick other = (EdgeTimetick) obj;
            return (ed.equals(other.ed)) && (t == other.t);    //To change body of overridden methods use File | Settings | File Templates.
        }

        @Override
        public int hashCode() {
            return 253523342 * ed.hashCode() + t;
        }

        @Override
        public String toString() {
            return ed.toString() + "->" + t;
        }

        public Edge ed;
        public int t;
    }
}
