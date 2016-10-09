import java.util.*;

public class GraphFast {
    public class ComparableCluster implements Comparable<ComparableCluster> {
        public Double s, n;

        public ComparableCluster(double _s, double _n) {
            s = _s;
            n = _n;
        }

        public int compareTo(ComparableCluster o) {
            return s.compareTo(o.s);
        }

        public String toString() {
            return "(" + s + "," + n + ")";
        }
    }

    public static class ComparableEdge implements Comparable<ComparableEdge> {
        public Double w;
        public Integer il;
        public Integer ih;
        public Integer v1;
        public Integer v2;

        public ComparableEdge(double _w, int _il, int _ih, int _v1, int _v2) {
            w = _w;
            il = _il;
            ih = _ih;
            v1 = _v1;
            v2 = _v2;
        }

        public int compareTo(ComparableEdge o) {
            if (w.compareTo(o.w) != 0) return w.compareTo(o.w);
            else return il.compareTo(o.il);
        }

        @Override
        public String toString() {
            return v1 + "-" + v2 + "[" + w + "]";
        }

    }

    static Random r = new Random();

    // -------------------- Fields

    // the second vertex id on an edge
    public int[] endv;

    // size n + 1: edges of vertex i are between endv[v[i]] and endv[v[i+1]]
    public int[] ind;

    // size m edge weights
    public double[] we;

    // size n node weights
    public double[] wn;

    // size n nodes names
    public String[] names;
    public HashMap<String, Integer> n2i = new HashMap<String, Integer>();

    public void buildN2I() {
        if (n2i.size() == 0) {
            for (int i = 0; i < getn(); i++)
                n2i.put(names[i], i);
        }
    }


    // Dual graph in which
    public GraphFast agg = null;
    // this is a mapping between super nodes and
    // their constituent nodes. positive edges can be then
    // retrieved in the original graph
    public ArrayList<String> aggnode_map = null;
    // this is a mapping between remaining non-positive edges
    // and the original edges.
    public ArrayList<Integer> aggedge_emap = null;

    // op score if read from the input
    public double opt_nw = -1;


    // ------------------------ Constructors --------------------------

    public GraphFast(String name, double value) {
        wn = new double[1];
        wn[0] = value;
        names = new String[1];
        names[0] = name;
    }

    public GraphFast(int n, int m) {
        ind = new int[n + 1];
        Arrays.fill(ind, 0);
        endv = new int[m];
        Arrays.fill(endv, 0);
        we = new double[m];
        Arrays.fill(we, Double.NEGATIVE_INFINITY);
        wn = new double[n];
        Arrays.fill(wn, Double.NEGATIVE_INFINITY);
        names = new String[n];
    }

    //shallow copy
    public GraphFast(GraphFast g) {
        ind = g.ind;
        endv = g.endv;
        we = g.we;
        wn = g.wn;
        names = g.names;
    }

    // Builds a graph from a list of edges. the edges should be in both directions if undirected
    // edges is in the format
    // name of edge1.nodes1, name of edge1.node2, name of edge2.node1 ....
    // node names is just the sequence of node names
    public GraphFast(String[] edges,
                     String[] node_names) {
        this(node_names.length, edges.length / 2);
        assert (edges.length % 2 == 0);

        // collect names of nodes
        HashMap<String, Integer> n2i = new HashMap<String, Integer>();
        for (int i = 0; i < node_names.length; i++) {
            n2i.put(node_names[i], i);
            names[i] = node_names[i];
        }
        // collect a list of edges for each node
        HashMap<Integer, ArrayList<Integer>> nbrs = new HashMap<Integer, ArrayList<Integer>>();
        int n1, n2;
        for (int i = 0; i < edges.length; i += 2) {
            n1 = n2i.get(edges[i]);
            n2 = n2i.get(edges[i + 1]);
            if (!nbrs.containsKey(n1)) nbrs.put(n1, new ArrayList<Integer>());
            nbrs.get(n1).add(n2);
        }
        // put the collected edges in the graph
        for (int i = 0; i < names.length; i++) {
            ind[i + 1] = ind[i] + nbrs.get(i).size();
            for (int j = 0; j < nbrs.get(i).size(); j++) {
                endv[ind[i] + j] = nbrs.get(i).get(j);
            }
        }
    }


    // Builds a graph from a list of edges
    public GraphFast(String[] edges,
                     double[] weights_e,
                     String[] node_names,
                     double[] weights_n) {
        this(node_names.length, weights_e.length);
        assert (edges.length % 2 == 0 &&
                weights_e.length == edges.length / 2 &&
                node_names.length == weights_n.length);

        // collect names of nodes
        HashMap<String, Integer> n2i = new HashMap<String, Integer>();
        for (int i = 0; i < node_names.length; i++) {
            n2i.put(node_names[i], i);
            names[i] = node_names[i];
            wn[i] = weights_n[i];
        }
        // collect a list of edges for each node
        HashMap<Integer, ArrayList<Integer>> nbrs = new HashMap<Integer, ArrayList<Integer>>();
        HashMap<Integer, ArrayList<Double>> nbrsw = new HashMap<Integer, ArrayList<Double>>();
        int n1, n2;
        for (int i = 0; i < edges.length; i += 2) {
            n1 = n2i.get(edges[i]);
            n2 = n2i.get(edges[i + 1]);
            if (!nbrs.containsKey(n1)) nbrs.put(n1, new ArrayList<Integer>());
            if (!nbrsw.containsKey(n1)) nbrsw.put(n1, new ArrayList<Double>());
            nbrs.get(n1).add(n2);
            nbrsw.get(n1).add(weights_e[i / 2]);
        }
        // put the collected edges in the graph
        for (int i = 0; i < names.length; i++) {
            ind[i + 1] = ind[i] + nbrs.get(i).size();
            for (int j = 0; j < nbrs.get(i).size(); j++) {
                endv[ind[i] + j] = nbrs.get(i).get(j);
                we[ind[i] + j] = nbrsw.get(i).get(j);
            }
        }
    }

    // ------------------------ Subgraph --------------------------

    public GraphFast getStructure() {
        GraphFast g = new GraphFast(getn(), getm());
        g.ind = Arrays.copyOf(ind, ind.length);
        g.endv = Arrays.copyOf(endv, endv.length);
        g.names = Arrays.copyOf(names, names.length);
        Arrays.fill(g.we, 0.0);
        Arrays.fill(g.wn, 0.0);
        return g;
    }

    public GraphFast subgraph(ArrayList<Integer> edges_ind) {
        int[] ei = new int[edges_ind.size()];
        for (int i = 0; i < edges_ind.size(); i++) ei[i] = edges_ind.get(i);
        return subgraph(ei);
    }

    // returns the indices of edges in this that are also in g
    // returns null if g is not a subgraph of this
    public ArrayList<Integer> getSubgraphEdgesStructureOnly(GraphFast g) {
        ArrayList<Integer> res = new ArrayList<Integer>();
        if (g.getn() == 1) return res;
        HashSet<String> gedges = new HashSet<String>();
        for (int i = 0; i < g.getn(); i++) {
            for (int j = g.ind[i]; j < g.ind[i + 1]; j++) {
                gedges.add(g.names[i] + "#" + g.names[g.endv[j]]);
            }
        }
        for (int i = 0; i < getn(); i++) {
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                if (gedges.contains(names[i] + "#" + names[endv[j]])) {
                    res.add(j);
                    gedges.remove(names[i] + "#" + names[endv[j]]);
                }
            }
        }
        if (gedges.size() > 0) return null;
        return res;
    }

    public GraphFast subgraph(int[] edge_indices) {
        String[] edges = new String[edge_indices.length * 2];
        double[] weights_e = new double[edge_indices.length];
        TreeSet<Integer> nodes = new TreeSet<Integer>();
        // compose the graph, corresponding to the induced tree
        for (int i = 0; i < edge_indices.length; i++) {
            nodes.add(endv[edge_indices[i]]);
            nodes.add(getOrigin(edge_indices[i]));
            edges[2 * i + 1] = names[endv[edge_indices[i]]];
            edges[2 * i] = names[getOrigin(edge_indices[i])];
            weights_e[i] = we[edge_indices[i]];
        }

        String[] node_names = new String[nodes.size()];
        double[] weights_n = new double[nodes.size()];
        int pos = 0;
        for (int n : nodes) {
            node_names[pos] = names[n];
            weights_n[pos] = wn[n];
            pos++;
        }
        return new GraphFast(edges, weights_e, node_names, weights_n);
    }

    public GraphFast subgraph(TreeSet<ComparableEdge> eds) {
        String[] edges = new String[eds.size() * 4];
        double[] weights_e = new double[eds.size() * 2];
        HashSet<Integer> nodes = new HashSet<Integer>();
        // compose the graph, corresponding to the induced tree
        int i = 0;
        for (ComparableEdge ce : eds) {
            nodes.add(ce.v1);
            nodes.add(ce.v2);
            edges[2 * i + 1] = names[ce.v1];
            edges[2 * i] = names[ce.v2];
            weights_e[i] = ce.w;
            i++;
            edges[2 * i + 1] = names[ce.v2];
            edges[2 * i] = names[ce.v1];
            weights_e[i] = ce.w;
            i++;
        }

        String[] node_names = new String[nodes.size()];
        double[] weights_n = new double[nodes.size()];
        int pos = 0;
        for (int n : nodes) {
            node_names[pos] = names[n];
            weights_n[pos] = wn[n];
            pos++;
        }
        return new GraphFast(edges, weights_e, node_names, weights_n);
    }

    public GraphFast inducedSubgraph(ArrayList<Integer> nodes) {
        // compute the number of edges
        int n = nodes.size();
        int m = 0;
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = ind[nodes.get(i)]; j < ind[nodes.get(i) + 1]; j++) {
                if (nodes.indexOf(endv[j]) >= 0) {
                    m++;
                }
            }
        }
        GraphFast g = new GraphFast(n, m);
        for (int i = 0; i < nodes.size(); i++) {
            g.names[i] = names[nodes.get(i)];
            g.wn[i] = wn[nodes.get(i)];
        }

        int eidx = 0;
        g.ind[0] = 0;
        for (int i = 0; i < nodes.size(); i++) {
            g.ind[i + 1] = g.ind[i];
            for (int j = ind[nodes.get(i)]; j < ind[nodes.get(i) + 1]; j++) {
                if (nodes.indexOf(endv[j]) >= 0) {
                    g.endv[eidx] = nodes.indexOf(endv[j]);
                    g.we[eidx] = we[j];
                    g.ind[i + 1]++;
                    eidx++;
                }
            }
        }
        return g;
    }

    public GraphFast getLargestConnectedComponent() {
        int[] comp = getConnectedComponents();
        HashMap<Integer, Integer> ccoutns = new HashMap<Integer, Integer>();
        for (int i = 0; i < getn(); i++) {
            if (!ccoutns.containsKey(comp[i])) ccoutns.put(comp[i], 1);
            else ccoutns.put(comp[i], 1 + ccoutns.get(comp[i]));
        }
        int lccid = -1, lccn = -1;
        for (Integer c : ccoutns.keySet()) {
            if (lccn < ccoutns.get(c)) {
                lccid = c;
                lccn = ccoutns.get(c);
            }
        }
        ArrayList<Integer> nodes = new ArrayList<Integer>();
        for (int i = 0; i < getn(); i++)
            if (comp[i] == lccid)
                nodes.add(i);
        return inducedSubgraph(nodes);
    }

    // ------------------------ Scores --------------------------
    public double getScore() {
        double res = 0;
        for (int i = 0; i < getn(); i++) {
            res += wn[i];
            if (getm() != 0) {
                for (int j = ind[i]; j < ind[i + 1]; j++) {
                    if (i < endv[j]) { // Count only once
                        res -= we[j];
                    }
                }
            }
        }
        return res;
    }

    public double getP() {
        double sum = 0;
        for (int i = 0; i < getn(); i++)
            sum += wn[i];
        return sum;
    }

    public double getN() {
        double res = 0;
        if (getm() != 0) {
            for (int i = 0; i < getn(); i++) {
                for (int j = ind[i]; j < ind[i + 1]; j++) {
                    if (i < endv[j]) {
                        res += we[j];
                    }
                }
            }
        }
        return res;
    }

    private double getSumOfPositiveEdges(String nodes) {
        double sum = 0;
        ArrayList<Integer> edges = getPositiveEdgesAmong(nodes);
        for (int ed : edges) {
            assert (we[ed] > 0);
            sum += we[ed];
        }
        return sum;
    }

    // ------------------------ Transform --------------------------

    // pass the names of a set of nodes and their collective new score
    // and get a graph with these nodes combined in a super node and
    // edges to another node: min(of all original edges)
    public GraphFast collapseCluster(ArrayList<String> cluster, double nw, String name) {
        // first collect the new edges and rename where needed
        String n1, n2, key;
        HashMap<String, Double> nedges = new HashMap<String, Double>();
        for (int i = 0; i < getn(); i++) {
            n1 = names[i];
            if (cluster.contains(n1)) n1 = name;
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                n2 = names[endv[j]];
                if (cluster.contains(n2)) n2 = name;
                if (n1.equals(n2)) continue;
                key = n1 + "#" + n2;
                if (nedges.containsKey(key)) {
                    if (we[j] < nedges.get(key)) {
                        nedges.put(key, we[j]);
                    }
                } else {
                    nedges.put(key, we[j]);
                }
            }
        }

        // set the names and wn
        int npos = 0, epos = 0;
        GraphFast g = new GraphFast(getn() - cluster.size() + 1, nedges.size());
        for (int i = 0; i < getn(); i++) {
            n1 = names[i];
            if (cluster.contains(n1)) continue;
            g.names[npos] = names[i];
            g.wn[npos] = wn[i];
            npos++;
        }
        g.names[npos] = name;
        g.wn[npos] = nw;
        // Now set the ind, endv and we
        for (int i = 0; i < g.getn(); i++) {
            n1 = g.names[i];
            for (String k : nedges.keySet()) {
                if (n1.equals(k.split("#")[0])) {
                    g.endv[epos] = g.getNodeIndex(k.split("#")[1]);
                    g.we[epos] = nedges.get(k);
                    epos++;
                    g.ind[i + 1] = epos;
                }
            }
        }
        return g;
    }


    // returns a mapping from edges to components
    public int[] getPositiveEdgeComponents() {
        int[] res = new int[getm()];

        // these are the takesn edges
        boolean[] te = new boolean[getm()];
        Arrays.fill(te, false);
        int[] clusters = new int[getn()];
        Arrays.fill(clusters, 0);
        //TODO
        int sind = 0;
        int cid = 0;
        double score;
        Integer current;


        Queue<Integer> q = new PriorityQueue<Integer>();
        while (true) {
            assert (q.isEmpty());
            // initialize a partition
            while (sind < getm() && (te[sind] || we[sind] < 0.0)) sind++;
            if (sind == getm()) break;
            cid++;
            score = 0.0;
            q.add(Math.min(getOrigin(sind), endv[sind]));
            while (null != (current = q.poll())) {
                clusters[current] = cid;
                for (int i = ind[current]; i < ind[current + 1]; i++) {
                    if (te[i]) continue;
                    if (we[i] >= 0.0) {
                        // put the other node tin q if not there yet and not considered
                        if (clusters[endv[i]] == 0 && !q.contains(endv[i]))
                            q.add(endv[i]);
                        if (current < endv[i]) {
                            te[i] = true;
                            te[getReverseEdgeIndex(i)] = true;
                            score += we[i];
                        }
                    } else { // negative edge
                        if (clusters[endv[i]] == clusters[current]) {
                            te[i] = true;
                            te[getReverseEdgeIndex(i)] = true;
                        }
                    }
                }
            }
//			if (!p.containsKey(score)) p.put(score, 1);
//			else p.put(score, p.get(score) + 1);
        }


        return res;
    }


    //TODO: optimize, use getPositiveConnectedComponents
    // Takes a general +/- edge-weighted graph and
    // aggregates connected positive components into
    // positively weighted super-nodes with a corresponding
    // score
    public GraphFast getAggNodeWeighted() {

        int[] map = new int[getn()];
        for (int i = 0; i < getn(); i++) map[i] = i;
        //TODO improve inefficient
        // map nodes to the smallest node they see through
        // a positive edge
        boolean change = true;
        while (change) {
            change = false;
            for (int i = 0; i < ind.length - 1; i++) {
                for (int j = ind[i]; j < ind[i + 1]; j++) {
                    if ((we[j] >= 0) &&
                            (map[endv[j]] != map[i])
                            ) {
                        map[i] = Math.min(map[endv[j]], map[i]);
                        change = true;
                    }
                }
            }
        }

        // Map these components to nodes
        // create the nodes and their weights
        HashMap<Integer, String> map_multi = new HashMap<Integer, String>();
        for (int i = 0; i < map.length; i++) {
            if (!map_multi.containsKey(map[i]))
                map_multi.put(map[i], "" + i);
            else
                map_multi.put(map[i], map_multi.get(map[i]) + "#" + i);
        }
        aggnode_map = new ArrayList<String>();
        for (String s : map_multi.values()) {
            aggnode_map.add(s);
        }

        HashMap<String, Integer> newnodes2oldedgeidx = new HashMap<String, Integer>();
        // Build edges among new nodes, map new negative edges to original edges
        int edge = -1;
        for (int i = 0; i < aggnode_map.size(); i++) {
            for (int j = i + 1; j < aggnode_map.size(); j++) {
                edge = getBestEdgeBetween(aggnode_map.get(i), aggnode_map.get(j));
                if (edge != -1) {
                    newnodes2oldedgeidx.put(i + "#" + j, edge);
                    newnodes2oldedgeidx.put(j + "#" + i, edge);
                }
            }
        }
        agg = new GraphFast(aggnode_map.size(), newnodes2oldedgeidx.size());
        aggedge_emap = new ArrayList<Integer>();
        int cnt = 0;
        int ed = -1;
        agg.ind[0] = 0;
        // build up the graph iterrating through the edgeset
        for (int i = 0; i < aggnode_map.size(); i++) {
            agg.wn[i] = getSumOfPositiveEdges(aggnode_map.get(i));
            agg.names[i] = "" + i;
            agg.ind[i + 1] = cnt;
            for (int j = 0; j < aggnode_map.size(); j++) {
                if (i == j) continue;
                if (newnodes2oldedgeidx.containsKey(i + "#" + j)) {
                    ed = newnodes2oldedgeidx.get(i + "#" + j);
                    assert (we[ed] < 0) : "Trying to add edge >=0 in the node-weighted graph";
                    agg.endv[cnt] = j;
                    agg.ind[i + 1]++;
                    agg.we[cnt] = (-1.0) * we[ed];
                    // finally update mapping and a count
                    aggedge_emap.add(ed);
                    cnt++;
                }
            }
        }
        return agg;
    }

    // Computes positive components, their scores and highest score negative edge adjacent to them
    public void getPosNegClusters(ArrayList<ComparableCluster> clusts) {
        // track if an edges is taken
        boolean[] processed = new boolean[getm()];
        Arrays.fill(processed, false);

        // the node to clusters mapping
        int[] clusters = new int[getn()];
        Arrays.fill(clusters, 0);

        int sind = 0;
        int cid = 0;
        double score, max_neg;
        Integer current;

        Queue<Integer> q = new PriorityQueue<Integer>();
        while (true) {
            assert (q.isEmpty());
            while (sind < getm()) {
                if (!processed[sind] && we[sind] >= 0.0) break;
                sind++;
            }
            if (sind == getm()) break;
            // new cluster init the score
            cid++;
            score = 0.0;
            q.add(Math.min(getOrigin(sind), endv[sind]));
            HashSet<Integer> adjacent = new HashSet<Integer>();
            while (null != (current = q.poll())) {
                clusters[current] = cid;
                adjacent.clear();
                for (int i = ind[current]; i < ind[current + 1]; i++) {
                    if (processed[i]) continue;
                    if (we[i] >= 0.0) {
                        processed[i] = true;
                        // put the other node tin q if not there yet and not considered
                        if (clusters[endv[i]] == 0 && !q.contains(endv[i]))
                            q.add(endv[i]);
                        // Take each edge only once. The score is counted only once
                        if (current < endv[i]) {
                            processed[getReverseEdgeIndex(i)] = true;
                            score += we[i];
                        }
                    } else { // negative edge
                        adjacent.add(i);
                    }
                }
            }
            max_neg = Double.NEGATIVE_INFINITY;
            for (Integer i : adjacent) if (we[i] > max_neg) max_neg = we[i];
            if (score + max_neg > 0)
                clusts.add(new ComparableCluster(score, max_neg));
        }
    }

    // Computes positive component score frequencies and negative connection weight frequencies
    // Counts the score of undirected edges only once (does not double count for the two edges)
    public void getPosNegDists(TreeMap<Double, Integer> p, TreeMap<Double, Integer> n) {
        // track if an edges is taken
        boolean[] processed = new boolean[getm()];
        Arrays.fill(processed, false);

        // the node to clusters mapping
        int[] clusters = new int[getn()];
        Arrays.fill(clusters, 0);

        int sind = 0;
        int cid = 0;
        double score;
        Integer current;

        Queue<Integer> q = new PriorityQueue<Integer>();
        while (true) {
            assert (q.isEmpty());
            // initialize a partition. scroll to the next positive not taken edge
            while (sind < getm()) {
                if (!processed[sind] && we[sind] >= 0.0) break;
                sind++;
            }
            if (sind == getm()) break;
            // new cluster init the score
            cid++;
            score = 0.0;

            q.add(Math.min(getOrigin(sind), endv[sind]));
            while (null != (current = q.poll())) {
                clusters[current] = cid;
                for (int i = ind[current]; i < ind[current + 1]; i++) {
                    if (processed[i]) continue; // skip taken
                    if (we[i] >= 0.0) {
                        processed[i] = true;
                        // put the adjacent node in q if not there yet and not considered
                        if (clusters[endv[i]] == 0 && !q.contains(endv[i]))
                            q.add(endv[i]);
                        // Take each edge only once. The score is counted only once
                        if (current < endv[i]) {
                            processed[getReverseEdgeIndex(i)] = true;
                            score += we[i];
                        }
                    } else { // negative edge. mark it processed if intracluster
                        if (clusters[endv[i]] == clusters[current]) {
                            processed[i] = true;
                            processed[getReverseEdgeIndex(i)] = true;
                        }
                    }
                }
            }
            if (!p.containsKey(score))
                p.put(score, 1);
            else
                p.put(score, p.get(score) + 1);
        }

        // collect not taken negative edges
        for (int i = 0; i < getm(); i++) {
            if (!processed[i]) {
                if (!n.containsKey(we[i])) n.put(we[i], 1);
                else n.put(we[i], n.get(we[i]) + 1);
            }
        }
    }

    // ------------------------ Search --------------------------

    public ArrayList<Integer> getPositiveEdgesAmong(String nodes) {
        ArrayList<Integer> res = new ArrayList<Integer>();
        int ed;
        if (nodes.split("#").length == 1) return res;
        HashSet<Integer> nodeidx = new HashSet<Integer>();
        for (String s : nodes.split("#")) nodeidx.add(Integer.parseInt(s));
        for (Integer n1 : nodeidx) {
            for (Integer n2 : nodeidx) {
                if (n2 <= n1) continue;
                if ((ed = getEdgeIndex(n1, n2)) != -1) {
                    if (we[ed] > 0) {
                        res.add(ed);
                    }
                }
            }
        }
        return res;
    }

    private int getBestEdgeBetween(String nodes1, String nodes2) {
        HashSet<Integer> nodesidx1 = new HashSet<Integer>();
        HashSet<Integer> nodesidx2 = new HashSet<Integer>();
        for (String s : nodes1.split("#")) nodesidx1.add(Integer.parseInt(s));
        for (String s : nodes2.split("#")) nodesidx2.add(Integer.parseInt(s));
        int max_eidx = -1, eidx = -1;
        double maxe = Double.NEGATIVE_INFINITY;
        for (int i : nodesidx1) {
            for (int j : nodesidx2) {
                eidx = getEdgeIndex(i, j);
                if (eidx >= 0) {
                    if (we[eidx] > maxe) {
                        maxe = we[eidx];
                        max_eidx = eidx;
                    }
                }
            }
        }
        return max_eidx;
    }

    // performs a search from n1 to n2 using only edges
    public boolean sees(int n1, int n2, ArrayList<Integer> edges) {
        ArrayList<Integer> nodes = new ArrayList<Integer>();
        nodes.add(n1);
        int idx = 0, curr;
        while (idx < nodes.size()) {
            curr = nodes.get(idx);
            for (int i = ind[curr]; i < ind[curr + 1]; i++) {
                if (!edges.contains(i) &&
                        !edges.contains(getReverseEdgeIndex(i))) continue;
                if (endv[i] == n2) return true;
                if (!nodes.contains(endv[i])) nodes.add(endv[i]);

            }
            idx++;
        }
        return false;
    }


    // ------------------------ Accessors --------------------------

    public void setEdgeWeights(double[] ew) {
        for (int i = 0; i < getm(); i++) we[i] = ew[i];
    }

    public boolean isUndirected() {
        if (null == ind) return true; // single node
        if (getm() % 2 != 0) return false;
        int pos;
        for (int i = 0; i < getn(); i++) {
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                if (endv[j] < i) continue;
                pos = getEdgeIndex(endv[j], i);
                if (-1 == pos) return false;
//				if(we[pos] != we[j]) 
//					return false;
            }
        }
        return true;
    }

    public boolean isConnected() {
        if (null == ind) return true;
        Queue<Integer> q = new PriorityQueue<Integer>();
        boolean[] visited = new boolean[getn()];
        Arrays.fill(visited, false);
        q.add(0);
        int n;
        while (!q.isEmpty()) {
            n = q.poll();
            visited[n] = true;
            for (int i = ind[n]; i < ind[n + 1]; i++) {
                if (!q.contains(endv[i]) && !visited[endv[i]])
                    q.add(endv[i]);
            }
        }
        for (int i = 0; i < getn(); i++)
            if (!visited[i])
                return false;
        return true;
    }

    public int[] getConnectedComponents() {
        Queue<Integer> q = new PriorityQueue<Integer>();
        int[] visited = new int[getn()];
        Arrays.fill(visited, -1);
        int n;
        int cid = 0;
        while (true) {
            for (int i = 0; i < getn(); i++) {
                if (visited[i] == -1) {
                    q.add(i);
                    break;
                }
            }
            if (q.isEmpty()) break;
            while (!q.isEmpty()) {
                n = q.poll();
                visited[n] = cid;
                for (int i = ind[n]; i < ind[n + 1]; i++) {
                    if (!q.contains(endv[i]) && visited[endv[i]] == -1)
                        q.add(endv[i]);
                }
            }
            cid++;
        }
        return visited;
    }

    public int getn() {
        return names.length;
    }

    public int getm() {
        if (null == endv) return 0;
        return endv.length;
    }

    public int getEdgeIndex(int v1, int v2) {
        int pos = -1;
        for (int i = ind[v1]; i < ind[v1 + 1]; i++) {
            if (endv[i] == v2) pos = i;
        }
        return pos;
    }

    public int getReverseEdgeIndex(int i) {
        //assert(isUndirected());
        return getEdgeIndex(endv[i], getOrigin(i));
    }

    public int getNodeIndex(String name) {
        if (n2i.size() == 0) buildN2I();
        return n2i.get(name);
//		for (int i = 0; i < names.length; i++) {
//			if (names[i].equals(name)) return i;
//		}
//		return -1;
    }

    public int getOrigin(int mind) {
        assert (mind < endv.length);
        int res = 0;
        while (mind >= ind[res + 1]) res++;
        return res;
    }

    public double getWe(int v1, int v2) {
        return we[getEdgeIndex(v1, v2)];
    }

    public void setWe(int v1, int v2, double w) {
        int pos = getEdgeIndex(v1, v2);
        we[pos] = w;
    }

    public double getWnByName(String name) {
        return wn[getNodeIndex(name)];
    }

    public int getMaxOutDegree() {
        int idx = 0, max = -1;
        if (null == ind) return 0;
        for (int i = 0; i < ind.length - 1; i++) {
            if (ind[i + 1] - ind[i] > max) {
                max = ind[i + 1] - ind[i];
                idx = i;
            }
        }
        return idx;
    }

    // ------------------------ TO STRING --------------------------

    public String edgesToString(int[] edge_indices) {
        String res = "";
        for (int i : edge_indices) {
            res += getOrigin(i) + "#" + endv[i] + "\n";
        }
        return res;
    }

    public String toStringStruct() {
        return "GraphFast \n [names= " + Arrays.toString(names) + "\n ind=  "
                + Arrays.toString(ind) + "\n endv=" + Arrays.toString(endv) + "\n we=" + Arrays.toString(we)
                + "\n wn=" + Arrays.toString(wn) + "]\n";
    }

    public String toStringNodes() {
        String res = "Nodes:\n";
        for (int i = 0; i < names.length; i++) {
            res += names[i] + "(" + wn[i] + "),";
        }
        return res;
    }

    @Override
    public String toString() {
        String res = "Nodes:\n";
        for (int i = 0; i < names.length; i++) {
            res += names[i] + "(" + wn[i] + "),";
        }
        res += "\nEdges:\n";
        if (getn() == 1) return res;
        for (int i = 0; i < ind.length - 1; i++) {
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                res += j + ":" + names[i] + "\t" + names[endv[j]] + "\t" + we[j] + "\n";
            }
        }
        return res;
    }

    // just for debugging and comparison
    public String getParsedOrderedNames() {
        TreeSet<String> res = new TreeSet<String>();
        for (int i = 0; i < getn(); i++) {
            String[] na = names[i].split(",");
            for (String n : na) res.add(n);
        }

        String result = "";
        for (String n : res) result += n + " ";
        return result;
    }

    public String getNames() {
        String name = names[0];
        for (int i = 1; i < getn(); i++) {
            name += "," + names[i];
        }
        return name;
    }


    public ArrayList<Integer> getIndices(ArrayList<String> nnames) {
        ArrayList<Integer> result = new ArrayList<Integer>();
        for (String nname : nnames) result.add(getNodeIndex(nname));
        return result;
    }

    public ArrayList<String> getPositives() {
        ArrayList<String> res = new ArrayList<String>();
        for (int i = 0; i < getn(); i++) {
            if (wn[i] > 0) res.add(names[i]);
        }
        return res;
    }

    public double getMaxScore() {
        double max_score = -1;
        for (int i = 0; i < getn(); i++) {
            if (wn[i] > max_score) max_score = wn[i];
        }
        return max_score;
    }

    public String getMaxScoreName() {
        double max_score = -1;
        String max_score_name = null;
        for (int i = 0; i < getn(); i++) {
            if (wn[i] > max_score) {
                max_score = wn[i];
                max_score_name = names[i];
            }
        }
        return max_score_name;
    }

    public ArrayList<String> getAllNames() {
        ArrayList<String> result = new ArrayList<String>();
        for (int i = 0; i < getn(); i++) {
            result.add(names[i]);
        }
        return result;
    }


    public HashSet<Integer> getNeighEdges(int j) {
        HashSet<Integer> res = new HashSet<Integer>();
        int from = getOrigin(j);
        int to = endv[j];
        for (int i = ind[from]; i < ind[from + 1]; i++) {
            res.add(i);
        }
        for (int i = ind[to]; i < ind[to + 1]; i++) {
            res.add(i);
        }
        res.remove(j);
        res.remove(getReverseEdgeIndex(j));
        return res;
    }

    public String neighToString(String name) {
        String res = name;
        int i = getNodeIndex(name);
        for (int j = ind[i]; j < ind[i + 1]; j++) {
            if (i < endv[j]) {
                res += " " + names[endv[j]] + ":" + we[j] + " ";
            }
        }
        return res;
    }

    public ArrayList<Integer> obtainUndirectedEdges() {
        ArrayList<Integer> res = new ArrayList<Integer>();
        for (int i = 0; i < getn(); i++) {
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                if (i < endv[j]) res.add(j);
            }
        }
        return res;
    }

    public int translateNegEdge(int edgeIndex) {
        return aggedge_emap.get(edgeIndex);
    }

    //---------------------------- STATS ----------------------------
    public TreeMap<Double, Integer> getNodeDistribution() {
        TreeMap<Double, Integer> tm = new TreeMap<Double, Integer>();
        for (int i = 0; i < getn(); i++) {
            if (!tm.containsKey(wn[i])) tm.put(wn[i], 1);
            else tm.put(wn[i], tm.get(wn[i]) + 1);
        }
        return tm;
    }

    public TreeMap<Double, Integer> getEdgeScoreDistribution() {
        TreeMap<Double, Integer> tm = new TreeMap<Double, Integer>();
        for (int i = 0; i < getm(); i++) {
            if (!tm.containsKey(we[i])) tm.put(we[i], 1);
            else tm.put(we[i], tm.get(we[i]) + 1);
        }
        return tm;
    }

    public TreeMap<Double, Integer> getNegEdgeScoreDistribution() {
        TreeMap<Double, Integer> tm = new TreeMap<Double, Integer>();
        for (int i = 0; i < getm(); i++) {
            if (we[i] >= 0) continue;
            if (!tm.containsKey(we[i])) tm.put(we[i], 1);
            else tm.put(we[i], tm.get(we[i]) + 1);
        }
        return tm;
    }

    public TreeMap<Double, Integer> getPosEdgeScoreDistribution() {
        TreeMap<Double, Integer> tm = new TreeMap<Double, Integer>();
        for (int i = 0; i < getm(); i++) {
            if (we[i] < 0) continue;
            if (!tm.containsKey(we[i])) tm.put(we[i], 1);
            else tm.put(we[i], tm.get(we[i]) + 1);
        }
        return tm;
    }

    public TreeMap<Double, Integer> getEdgePairsScoreDistribution() {
        TreeMap<Double, Integer> tm = new TreeMap<Double, Integer>();
        double key;
        for (int i = 0; i < getn(); i++) {
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                for (int k = j + 1; k < ind[i + 1]; k++) {
                    if (we[j] > we[k]) continue;
                    key = we[j] + we[k];
                    if (!tm.containsKey(key)) tm.put(key, 1);
                    else tm.put(key, tm.get(key) + 1);
                }
            }

        }
        return tm;
    }

    // For every node reports the sum of its adjacent positive edges
    public TreeMap<Double, Integer> getNodeScoreDistribution() {
        TreeMap<Double, Integer> tm = new TreeMap<Double, Integer>();
        double key;
        for (int i = 0; i < getn(); i++) {
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                for (int k = j + 1; k < ind[i + 1]; k++) {
                    if (we[j] > we[k]) continue;
                    key = we[j] + we[k];
                    if (!tm.containsKey(key)) tm.put(key, 1);
                    else tm.put(key, tm.get(key) + 1);
                }
            }
        }
        return tm;
    }

}
