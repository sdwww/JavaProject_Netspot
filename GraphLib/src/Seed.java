import java.util.*;

public class Seed extends HashMap<Subsequence, Subgraph> {

    public boolean equals(Object that) {
        if (this == that)
            return true;

        if (!(that instanceof Seed))
            return false;

        Seed thatS = (Seed) that;

        Edge e1 = null;
        Edge e2 = null;

        for (Subsequence cs : this.keySet())  // there should only be one edge in each of these
            for (Edge e : this.get(cs))
                e1 = e;

        for (Subsequence cs : thatS.keySet())
            for (Edge e : thatS.get(cs))
                e2 = e;

        return e1.equals(e2);
    }

    public int hashCode() {

        Edge e1 = null;
        for (Subsequence cs : this.keySet())  // there should only be one edge in each of these
            for (Edge e : this.get(cs))
                e1 = e;

        return e1.hashCode();
    }

}
