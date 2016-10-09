
public class Pattern {
    public Subgraph subgraph;
    public Subsequence subsequence;
    public double score;

    public Pattern(Subgraph g, Subsequence s) {
        subgraph = g;
        subsequence = s;
    }

    public Pattern() {

    }

    public boolean intersects(Pattern p) {
        for (int i : subsequence) {
            if (p.subsequence.contains(i)) {
                for (Edge e : subgraph)
                    if (p.subgraph.contains(e))
                        return true;
            }
        }
        return false;
    }

    public boolean equals(Object ot) {
        Pattern other = (Pattern) ot;
        return subgraph.equals(other.subgraph) && subsequence.equals(other.subsequence);
    }

    public int hashCode() {
        return subgraph.hashCode() * subsequence.hashCode();
    }
}
