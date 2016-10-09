import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class QueueEdgesTimeticks {

    public QueueEdgesTimeticks(int _edges, int _timeticks) {
        edges = _edges;
        timeticks = _timeticks;
        map = new HashMap<DynaGraph.EdgeTimetick, Values>();
        //fastMax = Float.NEGATIVE_INFINITY;
        //fastEd = new TreeSet<DynaGraph.EdgeTimetick>(new EdgeTimeTickComparator<DynaGraph.EdgeTimetick>());
    }


    public HashMap<Subsequence, Subgraph> GenerateHSMSS(DynaGraph dynaGraph, HSMSSTimes times) throws IOException {
        return GenPattern.multiHSMSSPattern(dynaGraph, times, this);
    }

    public HashMap<Subsequence, Subgraph> UpdateHSMSS(DynaGraph dynaGraph, Subgraph patGraph, Subsequence patSeq, HSMSSTimes times) throws IOException {
        return GenPattern.updateHSMSSPattern(dynaGraph, patGraph, patSeq, times, this);
    }

    public HashMap<Subsequence, Subgraph> GenerateHSMSS2(DynaGraph dynaGraph, HSMSSTimes times) throws IOException {
        return GenPattern.multiHSMSSPattern2(dynaGraph, times, this);
    }

    public HashMap<Subsequence, Subgraph> UpdateHSMSS2(DynaGraph dynaGraph, Subgraph patGraph, Subsequence patSeq, HSMSSTimes times) throws IOException {
        return GenPattern.updateHSMSSPattern2(dynaGraph, patGraph, patSeq, times, this);
    }

    public HashMap<Subsequence, Subgraph> GenerateHSMSS3(DynaGraph dynaGraph, HSMSSTimes times) throws IOException {
        return GenPattern.multiHSMSSPattern3(dynaGraph, times, this);
    }

    public HashMap<Subsequence, Subgraph> UpdateHSMSS3(DynaGraph dynaGraph, Subgraph patGraph, Subsequence patSeq, HSMSSTimes times) throws IOException {
        return GenPattern.updateHSMSSPattern3(dynaGraph, patGraph, patSeq, times, this);
    }

    public class Values {
        Values(float _valHS, float _valMSS) {
            valHS = _valHS;
            valMSS = _valMSS;
        }

        public float valHS;
        public float valMSS;

        public float getScore() {
            return Math.max(valHS, 0) * Math.max(valMSS, 0);
        }

        public float getHSScore() {
            return valHS;
        }

        public float getMSSScore() {
            return valMSS;
        }
    }

    private static HashMap<DynaGraph.EdgeTimetick, Values> map;
    int edges;
    int timeticks;
    // float fastMax;
    // TreeSet<DynaGraph.EdgeTimetick> fastEd;

    public void set(Edge ed, int t, float valHS, float valMSS) {
        DynaGraph.EdgeTimetick ced = new DynaGraph.EdgeTimetick(ed, t);
        map.put(ced, new Values(valHS, valMSS));
        //fastEd.add(ced);
    }

    public void changeHS(Edge ed, int t, float val) {
        DynaGraph.EdgeTimetick ced = new DynaGraph.EdgeTimetick(ed, t);
        //fastEd.remove(ced);
        Values vals = map.get(ced);
        vals.valHS = val;
        //fastEd.add(ced);
    }

    public void changeMSS(Edge ed, int t, float val) {
        DynaGraph.EdgeTimetick ced = new DynaGraph.EdgeTimetick(ed, t);
        //fastEd.remove(ced);
        map.get(ced).valMSS = val;
        //fastEd.add(ced);
    }

    public float getMaxValue() {
        float max = Float.NEGATIVE_INFINITY;
        DynaGraph.EdgeTimetick ed = null;
        for (DynaGraph.EdgeTimetick e : map.keySet()) {
            float sc = map.get(e).getScore();
            if (sc > max) {
                max = sc;
                ed = e;
            }
        }
        //System.out.println("Max score:" + max);
        return max;
    }

    public float getValueHS(DynaGraph.EdgeTimetick e) {
        return map.get(e).getHSScore();
    }

    public float getValueMSS(DynaGraph.EdgeTimetick e) {
        return map.get(e).getMSSScore();
    }

    public DynaGraph.EdgeTimetick getMax() {
        float max = Float.NEGATIVE_INFINITY;
        DynaGraph.EdgeTimetick ed = null;
        for (DynaGraph.EdgeTimetick e : map.keySet()) {
            float sc = map.get(e).getScore();
            if (sc > max) {
                max = sc;
                ed = e;
            }
        }
        //System.out.println("Max score:" + max);
        return ed;
    }


    public DynaGraph.EdgeTimetick getMaxMSS() {
        float max = Float.NEGATIVE_INFINITY;
        DynaGraph.EdgeTimetick ed = null;
        for (DynaGraph.EdgeTimetick e : map.keySet()) {
            float sc = map.get(e).getMSSScore();
            if (sc > max) {
                max = sc;
                ed = e;
            }
        }
        return ed;
    }


    public HashMap<Edge, Sequence> getSequence() {
        HashMap<Edge, Sequence> subs = new HashMap<Edge, Sequence>();

        for (DynaGraph.EdgeTimetick et : map.keySet()) {
            if (!subs.containsKey(et.ed)) {
                Sequence seq = new Sequence();
                seq.setSize(timeticks);
                subs.put(et.ed, seq);
            }

            subs.get(et.ed).set(et.t, map.get(et).valHS);
        }
        return subs;
    }
}
