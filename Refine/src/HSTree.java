import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class HSTree {
    private HashMap<Integer, HashSet<Edge>> aggEdges;
    private int root2; // the other endpoint of the root edge

    private Node root;
    private HashMap<Integer, Node> mapNodes = new HashMap<Integer, Node>();
    private int components_counter = 0;

    class Node {
        int id;
        double weight;
        int component;
        public double score;

        Node parent = null;

        public class Child {
            Node node;
            double weight;
            boolean isTakenUp = false;
            double scoreUp;
            boolean isTakenDown = false;
            double scoreDown;

            Child(Node _n, double _w) {
                weight = _w;
                node = _n;
            }
        }

        public HashSet<Child> children = new HashSet<Child>();

        Node(int _id, double w) {
            id = _id;
            weight = w;
            mapNodes.put(id, this);
        }

        public void setChildren(TreeFast tree) {
            HashSet<Integer> visited = new HashSet<Integer>();
            setChildren(tree, visited);
        }

        public void setChildren(TreeFast tree, Set<Integer> visited) {
            visited.add(id);
            int u = tree.n2i.get("" + id);
            for (int j = tree.ind[u]; j < tree.ind[u + 1]; j++) {
                assert (tree.levels[tree.endv[j]] == tree.levels[u] + 1) || (tree.levels[tree.endv[j]] == tree.levels[u] - 1);
                int v = tree.endv[j];
                double wn = tree.wn[v];
                double we = tree.we[j];

                int id = Integer.parseInt(tree.names[v]);
                if (visited.contains(id))
                    continue;
                Node child = new Node(id, wn);
                child.parent = this;
                children.add(new Child(child, we));
                child.setChildren(tree, visited);
            }
        }

        public double computeHSup() {
            double totalScore = 0.0;
            for (Child c : children) {
                double score = c.node.computeHSup() - c.weight;
                c.scoreUp = score;
                if (score > 0.0) {
                    totalScore += score;
                    c.isTakenUp = true;
                }
            }
            totalScore += weight;
            return totalScore;
        }

        public void computeHSdown(double score) {
            computeHSdown(score, components_counter++);
        }

        public void computeHSdown(double _score, int _component) {
            score = _score;
            component = _component;
            for (Child c : children) {
                double scoreC = score;
                if (c.isTakenUp)
                    scoreC -= c.scoreUp;
                scoreC -= c.weight;
                c.scoreDown = scoreC;
                if (scoreC > 0)
                    c.isTakenDown = true;

                scoreC = Math.max(0.0, scoreC);

                double totalChildrenScore = 0.0;
                for (Child c2 : c.node.children)
                    totalChildrenScore += Math.max(0.0, c2.scoreUp);

                int newcomp = component;
                if (!c.isTakenDown)
                    newcomp = components_counter++;

                c.node.computeHSdown(scoreC + totalChildrenScore + c.node.weight, newcomp);
            }
        }

        public void collect(HashMap<Integer, HashSet<Edge>> _aggEdges, Subgraph subgraph) {
            collectUp(_aggEdges, subgraph);
            collectDown(_aggEdges, subgraph);
        }

        public void collectUp(HashMap<Integer, HashSet<Edge>> _aggEdges, Subgraph subgraph) {
            if (parent == null)
                return;

            for (Child c : parent.children) {
                if (c.node == this) {
                    if (!c.isTakenDown)
                        return;
                    break;
                }
            }

            subgraph.add(new Edge(parent.id, id));
            subgraph.addAll(_aggEdges.get(parent.id));

            for (Child c : parent.children) {
                if (c.node != this && c.isTakenUp)
                    c.node.collectDown(_aggEdges, subgraph);
            }

            parent.collectUp(_aggEdges, subgraph);
        }

        public void collectDown(HashMap<Integer, HashSet<Edge>> _aggEdges, Subgraph subgraph) {
            subgraph.addAll(_aggEdges.get(id));
            for (Node.Child c : children) {
                if (c.isTakenUp) {
                    subgraph.add(new Edge(id, c.node.id));
                    c.node.collectDown(_aggEdges, subgraph);
                }
            }
        }

    }

    public void build(TreeFast tree, HashMap<Integer, HashSet<Edge>> _aggEdges, Edge e) {
        aggEdges = _aggEdges;
        root2 = e.dest();

        root = new Node(e.source(), tree.wn[tree.n2i.get("" + e.source())]);
        root.setChildren(tree);
    }

    public void build(TreeFast tree, HashMap<Integer, HashSet<Edge>> _aggEdges) {
        aggEdges = _aggEdges;
        root2 = -1;

        root = new Node(0, tree.wn[tree.n2i.get("" + 0)]);
        root.setChildren(tree);
    }

    private Node getNode(int i) {
        return mapNodes.get(i);
    }

    private Node getSourceNode(Edge e) {
        Node u = mapNodes.get(e.source());
        Node v = mapNodes.get(e.dest());

        Node sourcenode = null;

        for (Node.Child c : u.children)
            if (c.node == v)
                sourcenode = u;

        if (sourcenode == null)
            for (Node.Child c : v.children)
                if (c.node == u)
                    sourcenode = v;

        assert sourcenode != null : "getSourceNode it should be called only after getEdge return a valid object";

        return sourcenode;
    }

    private Node.Child getEdge(Edge e) {
        Node u = mapNodes.get(e.source());
        Node v = mapNodes.get(e.dest());

        Node.Child nodechild = null;

        for (Node.Child c : u.children)
            if (c.node == v)
                nodechild = c;

        if (nodechild == null)
            for (Node.Child c : v.children)
                if (c.node == u)
                    nodechild = c;

        return nodechild;
    }

    public Subgraph computeHS(Edge e) {
        //compute scores
        double score = root.computeHSup();
        root.computeHSdown(score);
        // collect edges
        Subgraph subgraph = new Subgraph();
        Node.Child nodechild = getEdge(e);
        double scoreE = 0.0;
        if (nodechild == null) {
            Node source = getNode(e.source());
            Node dest = getNode(e.dest());
            source.collect(aggEdges, subgraph);
            if (source.component != dest.component)
                dest.collect(aggEdges, subgraph);
        } else {
            int ndir = 0;
            if (nodechild.isTakenDown) {
                nodechild.node.collectUp(aggEdges, subgraph);
            }
            if (nodechild.isTakenUp) {
                nodechild.node.collectDown(aggEdges, subgraph);
            }
        }
        subgraph.add(e);
        subgraph.addAll(aggEdges.get(e.source()));
        subgraph.addAll(aggEdges.get(e.dest()));

        return subgraph;
    }

    public Graph computeHSall(Graph graph) {
        //compute scores
        double score = root.computeHSup();
        root.computeHSdown(score);

        // build graph scores
        Graph graphScore = new Graph();
        for (Edge e : graph.keySet()) {
            Node.Child nodechild = getEdge(e);
            double scoreE = 0.0;
            if (nodechild == null) {
                Node source = getNode(e.source());
                Node dest = getNode(e.dest());
                if (source.component == dest.component)
                    scoreE = Math.max(source.score, dest.score) + Math.min(0, graph.get(e));
                else
                    scoreE = source.score + dest.score + Math.min(0, graph.get(e));
            } else {
                int ndir = 0;
                if (nodechild.isTakenDown) {
                    ndir++;
                    scoreE += nodechild.scoreDown;
                }
                if (nodechild.isTakenUp) {
                    ndir++;
                    scoreE += nodechild.scoreUp;
                }
                if (ndir == 2)
                    scoreE += nodechild.weight; // this cost was considered twice
                if (ndir == 0) {
                    Node source = getSourceNode(e);
                    Node dest = nodechild.node;
                    scoreE = source.weight + dest.weight - nodechild.weight;
                }
            }
            graphScore.put(e, (float) scoreE);
        }

        return graphScore;
    }
}
