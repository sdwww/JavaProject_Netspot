import ilog.concert.IloException;
import ilog.concert.IloIntVarMap;
import ilog.concert.IloTupleSet;
import ilog.cplex.IloCplex;
import ilog.opl.*;
import ilog.opl_core.cppimpl.IloTuple;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;


public class Graph extends HashMap<Edge, Float> {

    HashMap<Edge, Float> neigh = null;
    HashMap<Edge, ArrayList<Edge>> localEdgeLists = null;
    float totalWeight = 0;

    public void add(Graph other) {
        for (Edge e : other.keySet()) {
            float w = 0;
            if (this.containsKey(e))
                w = this.get(e);
            w += other.get(e);
            this.put(e, w);
            totalWeight += other.get(e);
        }
    }

    public HashMap<Integer, HashSet<Edge>> getNodeEdgeAssociation() {
        HashSet<Edge> cur = null;
        HashMap<Integer, HashSet<Edge>> rez = new HashMap<Integer, HashSet<Edge>>();
        for (Edge e : this.keySet()) {
            int n1 = e.source();
            int n2 = e.dest();

            if (rez.containsKey(n1))
                cur = rez.get(n1);
            else
                cur = new HashSet<Edge>();
            cur.add(e);
            rez.put(n1, cur);

            if (rez.containsKey(n2))
                cur = rez.get(n2);
            else
                cur = new HashSet<Edge>();
            cur.add(e);
            rez.put(n2, cur);

        }

        return rez;
    }

    public HashSet<Edge> getNHopsAway(int n, HashSet<Edge> cur, HashMap<Edge, ArrayList<Edge>> edgeLists) {
        HashSet<Edge> rez = new HashSet<Edge>();
        rez.addAll(cur);
        int hops = 0;
        while (hops < n) {
            HashSet<Edge> next = new HashSet<Edge>();
            for (Edge e : cur) {
                for (Edge ne : edgeLists.get(e)) {
                    if (!rez.contains(ne)) {
                        rez.add(ne);
                        next.add(ne);
                    }
                }
            }
            cur = next;
            hops++;
        }
        return rez;
    }

    public void getLocalEdgeLists() {
        localEdgeLists = new HashMap<Edge, ArrayList<Edge>>();
        Graph g = this;
        for (Edge e : g.keySet()) {
            ArrayList<Edge> cur = new ArrayList<Edge>();
            for (Edge n : g.keySet()) {
                if (n.dest() == e.source() || n.source() == e.dest() || n.source() == e.source() || n.dest() == e.dest())
                    cur.add(n);
            }
            localEdgeLists.put(e, cur);
        }

    }

    @Override
    public Graph clone() {
        Graph g = new Graph();
        Iterator<Edge> it = this.keySet().iterator();
        while (it.hasNext()) {
            Edge edge = it.next();
            g.put(edge, new Float(get(edge)));
        }
        return g;
    }

    public float GetTotalWeight() {
        float w = 0;
        for (Edge e : this.keySet()) {
            w += this.get(e);
        }
        return w;
    }

    public float GetFastTotalWeight() {

        return totalWeight;
    }

    public Graph GetFastSubgraph(Subgraph subgraph) {
        Graph g = new Graph();
        float w = 0;
        for (Edge e : subgraph) {
            float val = this.get(e);
            g.put(e, val);
            w += val;
        }
        g.totalWeight = w;
        return g;
    }

    public float GetSubgraphWeight(Subgraph subgraph) {
        float w = 0;
        for (Edge e : subgraph) {
            w += this.get(e);
        }
        return w;
    }

    public Graph GetSubgraph(Subgraph subgraph) {
        Graph g = new Graph();
        for (Edge e : subgraph) {
            g.put(e, this.get(e));
        }
        return g;
    }


    public HashMap<Edge, Float> getEdgePositiveNeighbourScores(ArrayList<Edge> changes, HashMap<Edge, ArrayList<Edge>> edgeLists) {


        if (neigh == null) {
            changes = new ArrayList<Edge>(keySet());
            neigh = new HashMap<Edge, Float>();
        }

        if (changes == null || changes.size() == 0)
            return neigh;

        for (Edge e : changes) {
            float score = 0;
            for (Edge n : edgeLists.get(e))
                score += this.get(n);
            neigh.put(e, score);
        }

        return neigh;

    }

    // public HashMap<Edge,Float> getEdgePositiveNeighbourScores(HashMap<Edge, ArrayList<Edge>> edgeLists) {
    //    return getEdgePositiveNeighbourScores(null, edgeLists);
    // }

    public HashMap<Edge, Float> getEdgePositiveNeighbourScores() {


        if (neigh == null) {
            neigh = new HashMap<Edge, Float>();
        } else
            return neigh;

        if (localEdgeLists == null)
            getLocalEdgeLists();

        for (Edge e : keySet()) {
            float score = 0;
            for (Edge n : localEdgeLists.get(e))
                score += this.get(n);
            //if( neigh.containsKey(e) == false)
            //System.out.println("%" + neigh.containsKey(e));
            neigh.put(e, score);
        }

  /*  	for (Edge e: keySet()) {
            float score  = 0;
    		for (Edge n:keySet())
          if(e.source() == n.source() || e.source() == n.dest() || e.dest() == n.source() || e.dest() == n.dest())			
    			  score += this.get(n);
		//if( neigh.containsKey(e) == false)
			//System.out.println("%" + neigh.containsKey(e));
    		neigh.put(e, score);
    	}*/

        return neigh;

    }

    // Returns a HashMap Edge->score, where the score is the
    // score of the positive component involving the edge
    // First very inefficient version
    public void getEdgePositiveComponentScores(HashMap<Edge, Float> res) {
        // this stores a mapping from a representative node to the set of nodes in the group
        HashMap<Integer, HashSet<Integer>> rep2group = new HashMap<Integer, HashSet<Integer>>();
        // this stores reverse mapping from nodes to their representative
        HashMap<Integer, Integer> node2rep = new HashMap<Integer, Integer>();
        // this stores the mapping from representatives to the scores of their groups
        HashMap<Integer, Float> rep2score = new HashMap<Integer, Float>();
        int s, d;
        boolean changed = true;
        while (changed) {
            changed = false;
            // init hash maps
            for (Edge e : this.keySet()) {
                s = e.source();
                d = e.dest();
                // TODO(misael,razvan): we are working with +/- weights right?
                // in case we work with 0/1 dataset, we should change this condition
                // and later how w score
                if (this.get(e) >= 0) {
                    if (node2rep.containsKey(s) &&
                            node2rep.containsKey(d)) {
                        if (node2rep.get(s) != node2rep.get(d)) {  // only if they are in diff grp
                            // merge d's group in s's group
                            rep2group.get(node2rep.get(s)).addAll(rep2group.get(node2rep.get(d)));
                            // set the new representative to d's group
                            int cd = node2rep.get(d);
                            for (Integer n : rep2group.get(node2rep.get(d))) {
                                node2rep.put(n, node2rep.get(s));
                            }
                            // remove d's group
                            //rep2group.remove(cd);  it gets a null pointer on tgFull264 if this isn't commented out
                            changed = true;
                        }
                    } else if (node2rep.containsKey(s)) { // only s is in a group
                        //System.out.println(rep2group.containsKey(node2rep.get(s)));
                        rep2group.get(node2rep.get(s)).add(d);
                        node2rep.put(d, node2rep.get(s));
                        changed = true;
                    } else if (node2rep.containsKey(d)) { // only d is in a group
                        rep2group.get(node2rep.get(d)).add(s);
                        node2rep.put(s, node2rep.get(d));
                        changed = true;
                    } else { // none of the two are in a group
                        rep2group.put(s, new HashSet<Integer>());
                        rep2group.get(s).add(s);
                        rep2group.get(s).add(d);
                        node2rep.put(s, s);
                        node2rep.put(d, s);
                        changed = true;
                    }
                }
                // we do not consider mapping nodes with negative edges
                // later we will just set edges with non-mapped nodes
                // to their actual value
            }
        }
        // now find the scores for each cluster
        for (Integer n : rep2group.keySet()) {
            rep2score.put(n, 0.0f);
        }
        for (Edge e : keySet()) {
            s = e.source();
            d = e.dest();
            // TODO(misael,razvan): same check as above
            if (this.get(e) >= 0) {
                rep2score.put(node2rep.get(s), rep2score.get(node2rep.get(s)) + get(e));
            }
        }
        // fill the mapping from edges to scores of components in res
        res.clear();
        for (Edge e : keySet()) {
            s = e.source();
            // TODO(misael,razvan): same check as above
            if (this.get(e) >= 0) {
                res.put(e, rep2score.get(node2rep.get(s)));
            } else {
                res.put(e, get(e));
            }
        }
    }

    public Subgraph MaximumScoreSubgraph() throws IOException {
        Subgraph subgraph = new Subgraph();
        String datFileName = "graph.dat";
        Writer datFile = new OutputStreamWriter(new FileOutputStream(datFileName), Config.fEncoding);
        try {
            datFile.write("tuples = {\n");
            int n = -1;
            for (Edge e : this.keySet()) {
                datFile.write("  <" + e.source() + " " + e.dest() + " " + this.get(e) + ">,\n");
                if (e.source() > n)
                    n = e.source();
                if (e.dest() > n)
                    n = e.dest();
            }
            n++;
            datFile.write("};\n");
            datFile.write("n=" + n + ";\n");
        } finally {
            datFile.close();
        }

        int status = 127;
        try {
            IloOplFactory.setDebugMode(false);
            IloOplFactory oplF = new IloOplFactory();
            IloOplErrorHandler errHandler = oplF.createOplErrorHandler();
            IloOplModelSource modelSource = oplF.createOplModelSource("MS.mod");
            IloOplSettings settings = oplF.createOplSettings(errHandler);
            IloOplModelDefinition def = oplF.createOplModelDefinition(modelSource, settings);
            IloCplex cplex = oplF.createCplex();
            cplex.setOut(null);
            IloOplModel opl = oplF.createOplModel(def, cplex);
            IloOplDataSource dataSource = oplF.createOplDataSource(datFileName);
            opl.addDataSource(dataSource);
            opl.generate();
            if (cplex.solve()) {
                System.out.println("OBJECTIVE: " + opl.getCplex().getObjValue());
                //opl.postProcess();
                IloTupleSet edges = opl.getElement("edges").asTupleSet();
                IloIntVarMap y = opl.getElement("y").asIntVarMap();
                for (Iterator it = edges.iterator(); it.hasNext(); ) {
                    IloTuple tuple = (IloTuple) it.next();
                    int u = tuple.getIntValue("u");
                    int v = tuple.getIntValue("v");
                    if ((opl.getCplex().getValue(y.get(tuple)) == 1) && (u <= v)) {
                        subgraph.add(new Edge(u - 1, v - 1));
                        //System.out.println("" + (u-1) + "," + (v-1));
                    }
                }
                //opl.printSolution(System.out);
            } else {
                System.out.println("No solution!");
            }
            oplF.end();
            status = 0;
        } catch (IloOplException e) {
            e.printStackTrace();
            status = 2;
        } catch (IloException e) {
            e.printStackTrace();
            status = 3;
        } catch (Exception e) {
            e.printStackTrace();
            status = 4;
        }

        return subgraph;
    }

    public Subgraph MaximumScoreSubgraph(Edge pivot) throws IOException {
        Subgraph subgraph = new Subgraph();
        String datFileName = "graph.dat";
        Writer datFile = new OutputStreamWriter(new FileOutputStream(datFileName), Config.fEncoding);
        try {
            datFile.write("tuples = {\n");
            int n = -1;
            for (Edge e : this.keySet()) {
                datFile.write("  <" + e.source() + " " + e.dest() + " " + this.get(e) + ">,\n");
                if (e.source() > n)
                    n = e.source();
                if (e.dest() > n)
                    n = e.dest();
            }
            n++;
            datFile.write("};\n");
            datFile.write("n=" + n + ";\n");
            datFile.write("pivot=<" + Math.min(pivot.source(), pivot.dest()) + " " + Math.max(pivot.source(), pivot.dest()) + ">;\n");
        } finally {
            datFile.close();
        }

        int status = 127;
        try {
            IloOplFactory.setDebugMode(false);
            IloOplFactory oplF = new IloOplFactory();
            IloOplErrorHandler errHandler = oplF.createOplErrorHandler();
            IloOplModelSource modelSource = oplF.createOplModelSource("MS_pivot.mod");
            IloOplSettings settings = oplF.createOplSettings(errHandler);
            IloOplModelDefinition def = oplF.createOplModelDefinition(modelSource, settings);
            IloCplex cplex = oplF.createCplex();
            cplex.setOut(null);
            IloOplModel opl = oplF.createOplModel(def, cplex);
            IloOplDataSource dataSource = oplF.createOplDataSource(datFileName);
            opl.addDataSource(dataSource);
            opl.generate();
            if (cplex.solve()) {
                System.out.println("OBJECTIVE: " + opl.getCplex().getObjValue());
                //opl.postProcess();
                IloTupleSet edges = opl.getElement("edges").asTupleSet();
                IloIntVarMap y = opl.getElement("y").asIntVarMap();
                for (Iterator it = edges.iterator(); it.hasNext(); ) {
                    IloTuple tuple = (IloTuple) it.next();
                    int u = tuple.getIntValue("u");
                    int v = tuple.getIntValue("v");
                    if ((opl.getCplex().getValue(y.get(tuple)) == 1) && (u <= v)) {
                        subgraph.add(new Edge(u - 1, v - 1));
                        //System.out.println("" + (u-1) + "," + (v-1));
                    }
                }
                //opl.printSolution(System.out);
            } else {
                System.out.println("No solution!");
            }
            oplF.end();
            status = 0;
        } catch (IloOplException e) {
            e.printStackTrace();
            status = 2;
        } catch (IloException e) {
            e.printStackTrace();
            status = 3;
        } catch (Exception e) {
            e.printStackTrace();
            status = 4;
        }

        return subgraph;
    }

    public Subgraph MaximumScoreSubgraph_() throws IOException {
        Subgraph subgraph = new Subgraph();
        String datFileName = "graph.dat";
        Writer datFile = new OutputStreamWriter(new FileOutputStream(datFileName), Config.fEncoding);
        try {
            datFile.write("tuples = {\n");
            int n = -1;
            for (Edge e : this.keySet()) {
                datFile.write("  <" + e.source() + " " + e.dest() + " " + this.get(e) + ">,\n");
                if (e.source() > n)
                    n = e.source();
                if (e.dest() > n)
                    n = e.dest();
            }
            n++;
            datFile.write("};\n");
            datFile.write("n=" + n + ";\n");
        } finally {
            datFile.close();
        }

        try {
            String cmd = Config.oplrun + " MS.mod " + datFileName;
            Runtime run = Runtime.getRuntime();
            Process pr = run.exec(cmd);
            pr.waitFor();
            BufferedReader buf = new BufferedReader(new InputStreamReader(pr.getInputStream()));
            String line;
            int status = 0;
            while ((line = buf.readLine()) != null && status != 3) {
                switch (status) {
                    case 0:
                        if (line.compareTo("<<< solve") == 0)
                            status = 1;
                        break;
                    case 1:
                        if (line.length() > 9 && line.substring(0, 9).compareTo("OBJECTIVE") == 0)
                            status = 2;
                        break;
                    case 2:
                        if (line.compareTo("<<< post process") == 0)
                            status = 3;
                        else {
                            Edge e = Subgraph.processLine(line);
                            if (e != null)
                                subgraph.add(e);
                        }
                        break;
                }
            }
        } catch (Exception e) {
            System.out.println("Unable to run command oplrun. Check Config.java property oplrun");
        } finally {
            datFile.close();
        }
        return subgraph;
    }

    public Subgraph MaximumScoreSubgraph_(Edge pivot) throws IOException {
        Subgraph subgraph = new Subgraph();
        String datFileName = "graph.dat";
        Writer datFile = new OutputStreamWriter(new FileOutputStream(datFileName), Config.fEncoding);
        try {
            datFile.write("tuples = {\n");
            int n = -1;
            for (Edge e : this.keySet()) {
                datFile.write("  <" + e.source() + " " + e.dest() + " " + this.get(e) + ">,\n");
                if (e.source() > n)
                    n = e.source();
                if (e.dest() > n)
                    n = e.dest();
            }
            n++;
            datFile.write("};\n");
            datFile.write("n=" + n + ";\n");
            datFile.write("pivot=<" + Math.min(pivot.source(), pivot.dest()) + " " + Math.max(pivot.source(), pivot.dest()) + ">;\n");
        } finally {
            datFile.close();
        }

        try {
            String cmd = Config.oplrun + " MS_pivot.mod " + datFileName;
            Runtime run = Runtime.getRuntime();
            Process pr = run.exec(cmd);
            pr.waitFor();
            BufferedReader buf = new BufferedReader(new InputStreamReader(pr.getInputStream()));
            String line;
            int status = 0;
            while ((line = buf.readLine()) != null && status != 3) {
                switch (status) {
                    case 0:
                        if (line.compareTo("<<< solve") == 0)
                            status = 1;
                        break;
                    case 1:
                        if (line.length() > 9 && line.substring(0, 9).compareTo("OBJECTIVE") == 0)
                            status = 2;
                        break;
                    case 2:
                        if (line.compareTo("<<< post process") == 0)
                            status = 3;
                        else {
                            Edge e = Subgraph.processLine(line);
                            if (e != null)
                                subgraph.add(e);
                        }
                        break;
                }
            }
        } catch (Exception e) {
            System.out.println("Unable to run command oplrun. Check Config.java property oplrun");
        } finally {
            datFile.close();
        }
        return subgraph;
    }

    // this is a unit test for positive component scores
    private static boolean testGetPositiveComponentScores() {
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

        HashMap<Edge, Float> poscompscores = new HashMap<Edge, Float>();
        g.getEdgePositiveComponentScores(poscompscores);

        // test if actpcs is the same as poscompres
        assert (poscompscores.size() == actpcs.size()) : "Wrong number of positive components";
        for (Edge ed : g.keySet()) {
            assert (poscompscores.containsKey(ed)) :
                    "Edge " + ed.toString() + " missing from positive components\n";
            assert (poscompscores.get(ed).equals(actpcs.get(ed))) :
                    "Wrong pos comp score for edge " + ed.toString() + "= " + poscompscores.get(ed) +
                            ". Should be " + actpcs.get(ed) + " instead.\n";
        }

        return true;
    }

    public static void main(String[] args) throws IOException {
        if (testGetPositiveComponentScores()) {
            System.out.print("Tested GraphLib.Graph successfully.\n");
        }
        ;
    }
}
