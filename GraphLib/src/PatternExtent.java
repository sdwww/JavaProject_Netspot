import java.util.HashSet;
import java.util.Iterator;

public class PatternExtent extends HashSet<PatternExtent.Entry> {
    public class Entry {
        public Edge edge;
        public int time;

        public Entry(Edge ed, int t) {
            edge = ed;
            time = t;
        }

        @Override
        public boolean equals(Object obj) {
            Entry other = (Entry) obj;
            return edge.equals(other.edge) && (time == other.time);    //To change body of overridden methods use File | Settings | File Templates.
        }

        @Override
        public int hashCode() {
            return edge.hashCode() + 878765879 * time;
        }
    }

    private DynaGraph network;

    public PatternExtent(DynaGraph net, Subgraph patGraph, Subsequence patSeq) {
        network = net;
        for (Iterator<Edge> itEd = patGraph.iterator(); itEd.hasNext(); ) {
            Edge ed = itEd.next();
            for (Iterator<Integer> itT = patSeq.iterator(); itT.hasNext(); ) {
                int t = itT.next();
                if (network.get(t).get(ed) > 0.0)
                    add(new Entry(ed, t));
            }
        }
    }

    public float GetEnergy() {
        float energy = 0;
        for (Iterator<Entry> it = iterator(); it.hasNext(); ) {
            Entry entry = it.next();
            float w = network.get(entry.time).get(entry.edge);
            if (w > 0) {
                energy += w;
            }
        }
        return energy;
    }

    public float GetEnergy(Subgraph patGraph, Subsequence patSeq) {
        float energy = 0;
        for (Iterator<Edge> itEd = patGraph.iterator(); itEd.hasNext(); ) {
            Edge ed = itEd.next();
            for (Iterator<Integer> itT = patSeq.iterator(); itT.hasNext(); ) {
                int t = itT.next();
                if (this.contains(new Entry(ed, t))) {
                    float w = network.get(t).get(ed);
                    if (w > 0) {
                        energy += w;
                    }
                }
            }
        }
        return energy;
    }

    public Entry PickEntry() {
        Iterator<Entry> it = iterator();
        return it.next();
    }

    public void Subtract(Subgraph patGraph, Subsequence patSeq) {
        for (Iterator<Edge> itEd = patGraph.iterator(); itEd.hasNext(); ) {
            Edge ed = itEd.next();
            for (Iterator<Integer> itT = patSeq.iterator(); itT.hasNext(); ) {
                int t = itT.next();
                Iterator<Entry> itEntry = this.iterator();
                remove(new Entry(ed, t));
            }
        }
    }
}
