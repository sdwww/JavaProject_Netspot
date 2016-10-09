
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class GenPattern {

    public static TreeSet<SeedEntry> seeds = null;

    public static void maxPattern(DynaGraph g, int nrGr, int nrSeq, String outFileGr, String outFileSeq) throws IOException {
        BufferedWriter out1 = new BufferedWriter(new FileWriter(new File(outFileGr)));
        BufferedWriter out2 = new BufferedWriter(new FileWriter(new File(outFileSeq)));

        HashMap<Integer, Float> seqs = new HashMap<Integer, Float>();
        HashMap<Edge, Float> grs = new HashMap<Edge, Float>();

        seqs.put(0, g.get(0).GetTotalWeight());

        for (Edge e : g.get(0).keySet())
            grs.put(e, g.get(0).get(e));

        for (int i = 1; i < g.size(); i++) {
            seqs.put(i, g.get(i).GetTotalWeight());

            for (Edge e : g.get(i).keySet())
                grs.put(e, grs.get(e) + g.get(i).get(e));
        }
        List<Map.Entry<Integer, Float>> listSeq = new ArrayList<Map.Entry<Integer, Float>>(seqs.entrySet());
        List<Map.Entry<Edge, Float>> listGr = new ArrayList<Map.Entry<Edge, Float>>(grs.entrySet());

        Collections.sort(listSeq, new ValueComparator<Integer, Float>());
        Collections.sort(listGr, new ValueComparator<Edge, Float>());


        for (int i = 0; i < nrGr; i++) {
            Edge ce = listGr.get(i).getKey();
            out1.append(ce.source() + "," + ce.dest() + "\n");
        }
        for (int i = 0; i < nrSeq; i++) {
            out2.append(listSeq.get(i).getKey() + "\n");
        }

        out1.close();
        out2.close();
    }

    public static class ValueComparator<K, V extends Comparable<? super V>> implements Comparator<Map.Entry<K, V>> {

        public int compare(Map.Entry<K, V> a, Map.Entry<K, V> b) {
            return a.getValue().compareTo(b.getValue()) * -1;

        }
    }

    public static void randPattern(DynaGraph g, String outFileGr, String outFileSeq) throws IOException {
        Random r = new Random();
        int nrSeq = r.nextInt(g.size());
        int nrGr = r.nextInt(g.get(0).size());
        randPattern(g, nrGr, nrSeq, outFileGr, outFileSeq, r);
    }

    public static void randPattern(DynaGraph g, int nrGr, int nrSeq, String outFileGr, String outFileSeq) throws IOException {
        randPattern(g, nrGr, nrSeq, outFileGr, outFileSeq, new Random());
    }


    public static HashMap<Subsequence, Subgraph> randPattern(DynaGraph g, Random r) throws IOException {
        ArrayList<Integer> seqs = new ArrayList<Integer>();
        ArrayList<Edge> grs = new ArrayList<Edge>();

        for (Edge e : g.get(0).keySet()) {
            grs.add(e);
        }

        Subsequence ss = new Subsequence();
        int t = r.nextInt(g.size());
        ss.add(t);

        int cur = r.nextInt(grs.size());
        Edge e = grs.get(cur);

        Subgraph mgraph = new Subgraph();
        mgraph.add(e);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(ss, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> maxPattern(DynaGraph g) throws IOException {
        HashMap<Integer, Float> seqs = new HashMap<Integer, Float>();
        HashMap<Edge, Float> grs = new HashMap<Edge, Float>();

        seqs.put(0, g.get(0).GetTotalWeight());

        for (Edge e : g.get(0).keySet())
            grs.put(e, g.get(0).get(e));

        for (int i = 1; i < g.size(); i++) {
            seqs.put(i, g.get(i).GetTotalWeight());

            for (Edge e : g.get(i).keySet())
                grs.put(e, grs.get(e) + g.get(i).get(e));
        }
        List<Map.Entry<Integer, Float>> listSeq = new ArrayList<Map.Entry<Integer, Float>>(seqs.entrySet());
        List<Map.Entry<Edge, Float>> listGr = new ArrayList<Map.Entry<Edge, Float>>(grs.entrySet());

        int t = Collections.min(listSeq, new ValueComparator<Integer, Float>()).getKey();// with this comparator it corresponds to the max
        Edge e = Collections.min(listGr, new ValueComparator<Edge, Float>()).getKey();

        Subsequence ss = new Subsequence();
        ss.add(t);

        Subgraph mgraph = new Subgraph();
        mgraph.add(e);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(ss, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> max2Pattern(DynaGraph g) throws IOException {
        float max = Float.NEGATIVE_INFINITY;
        Edge maxE = null;
        int maxT = 0;
        for (int i = 0; i < g.size(); i++) {
            for (Edge e : g.get(i).keySet()) {
                float w = g.get(i).get(e);
                if (w > max) {
                    max = w;
                    maxE = e;
                    maxT = i;
                }
            }
        }

        Subsequence ss = new Subsequence();
        ss.add(maxT);
        Subgraph mgraph = new Subgraph();
        mgraph.add(maxE);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(ss, mgraph);

        return hm;
    }

    public static void randPattern(DynaGraph g, int nrGr, int nrSeq, String outFileGr, String outFileSeq, Random r) throws IOException {
        BufferedWriter out1 = new BufferedWriter(new FileWriter(new File(outFileGr)));
        BufferedWriter out2 = new BufferedWriter(new FileWriter(new File(outFileSeq)));

        ArrayList<Integer> seqs = new ArrayList<Integer>();
        ArrayList<Edge> grs = new ArrayList<Edge>();

        for (int i = 0; i < g.size(); i++) {
            seqs.add(i);
        }

        for (Edge e : g.get(0).keySet()) {
            grs.add(e);
        }

        for (int i = 0; i < nrSeq; i++) {
            int cur = r.nextInt(seqs.size());
            out2.append(seqs.get(cur) + "\n");
            seqs.remove(cur);
        }

        for (int i = 0; i < nrGr; i++) {
            int cur = r.nextInt(grs.size());
            Edge ce = grs.get(cur);
            out1.append(ce.source() + "," + ce.dest() + "\n");
            grs.remove(cur);
        }

        out1.close();
        out2.close();
    }


    public static HashMap<Subsequence, Float> allMSS(float[] a) {
        //System.out.println("MSS Size: " + a.length);
        float[] rez1 = new float[a.length];
        int[] start1 = new int[a.length];
        rez1[0] = a[0];
        start1[0] = 0;
        for (int i = 1; i < a.length; i++) {
            if (rez1[i - 1] + a[i] > a[i]) {
                rez1[i] = rez1[i - 1] + a[i];
                start1[i] = start1[i - 1];
            } else {
                rez1[i] = a[i];
                start1[i] = i;
            }
        }

        float[] rez2 = new float[a.length];
        int[] start2 = new int[a.length];
        rez2[a.length - 1] = a[a.length - 1];
        start2[a.length - 1] = a.length - 1;
        for (int i = a.length - 2; i >= 0; i--) {
            if (rez2[i + 1] + a[i] > a[i]) {
                rez2[i] = rez2[i + 1] + a[i];
                start2[i] = start2[i + 1];
            } else {
                rez2[i] = a[i];
                start2[i] = i;
            }

        }
        float mq = Float.NEGATIVE_INFINITY;
        int mi = 0;
        for (int i = 0; i < a.length; i++) {
            if (rez1[i] + rez2[i] - a[i] > mq) {
                mq = rez1[i] + rez2[i] - a[i];
                mi = i;
                //System.out.println("mq is: " + mq);
            }

        }
        Subsequence seq = new Subsequence();


        for (int i = start1[mi]; i <= start2[mi]; i++) {
            seq.add(i);
        }

        HashMap<Subsequence, Float> hm = new HashMap<Subsequence, Float>();
        hm.put(seq, mq);
        return hm;
    }

    public static float[] allPivotMSS(Sequence a) {
        float[] rez1 = new float[a.size()];

        rez1[0] = a.get(0);
        for (int i = 1; i < a.size(); i++)
            rez1[i] = Math.max(rez1[i - 1] + a.get(i), a.get(i));

        float[] rez2 = new float[a.size()];
        rez2[a.size() - 1] = a.get(a.size() - 1);
        for (int i = a.size() - 2; i >= 0; i--) {
            rez2[i] = Math.max(rez2[i + 1] + a.get(i), a.get(i));

        }
        float mq = Float.NEGATIVE_INFINITY;
        int mi = 0;
        for (int i = 0; i < a.size(); i++)
            rez1[i] = rez1[i] + rez2[i] - a.get(i);

        return rez1;

    }

    public static float[] allPivotMSS(float[] a) {
        float[] rez1 = new float[a.length];

        rez1[0] = a[0];
        for (int i = 1; i < a.length; i++)
            rez1[i] = Math.max(rez1[i - 1] + a[i], a[i]);

        float[] rez2 = new float[a.length];
        rez2[a.length - 1] = a[a.length - 1];
        for (int i = a.length - 2; i >= 0; i--) {
            rez2[i] = Math.max(rez2[i + 1] + a[i], a[i]);

        }
        float mq = Float.NEGATIVE_INFINITY;
        int mi = 0;
        for (int i = 0; i < a.length; i++)
            rez1[i] = rez1[i] + rez2[i] - a[i];

        return rez1;

    }

    public static HashMap<Subsequence, Float> bruteMSS(float[] a) {
        HashMap<Subsequence, Float> rez = new HashMap<Subsequence, Float>();
        Sequence b = new Sequence();
        for (int i = 0; i < a.length; i++)
            b.add(a[i]);
        //System.out.println(b);
        float mq = Float.NEGATIVE_INFINITY;
        Subsequence ms = null;
        for (int i = 0; i < a.length; i++) {
            float q = 0;
            Subsequence seq = b.MaximumScoreSubsequence(i);
            //System.out.println(seq);
            for (int s : seq) {
                //	System.out.print(a[s] + " ");
                q += a[s];
            }
            //System.out.println(q + " " + mq);
            if (q > mq) {
                mq = q;
                ms = seq;
            }
        }
        rez.put(ms, mq);
        return rez;
    }

    public static void conPattern(DynaGraph g, String outFileGr, String outFileSeq) throws IOException {
        BufferedWriter out1 = new BufferedWriter(new FileWriter(new File(outFileGr)));
        BufferedWriter out2 = new BufferedWriter(new FileWriter(new File(outFileSeq)));
        HashMap<Edge, float[]> comps = g.GetMatrixPositiveConnectedComponent();

        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        HashMap<Subsequence, Float> curm = null;
        for (Edge e : comps.keySet()) {
            curm = allMSS(comps.get(e));
            for (Subsequence cs : curm.keySet())
                if (curm.get(cs) > ms) {
                    ms = curm.get(cs);
                    mseq = cs;
                    me = e;
                }
        }

        out1.append(me.source() + "," + me.dest() + "\n");

        for (int i : mseq) {
            out2.append(i + "\n");
        }

        out1.close();
        out2.close();
    }

    public static HashMap<Subsequence, Subgraph> conPattern(DynaGraph g) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveConnectedComponent();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        HashMap<Subsequence, Float> curm = null;
        for (Edge e : comps.keySet()) {
            curm = allMSS(comps.get(e));
            for (Subsequence cs : curm.keySet())
                if (curm.get(cs) > ms) {
                    ms = curm.get(cs);
                    mseq = cs;
                    me = e;
                }
        }
        System.out.println(ms);
        for (int i : mseq)
            System.out.print(comps.get(me)[i] + "f, ");
        System.out.println();

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;

    }

    public static HashMap<Subsequence, Subgraph> multConPattern(DynaGraph g) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveConnectedComponent();
        HashMap<Edge, Sequence> subs = g.getAllSubsequences();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        float[] curm = null;

        for (Edge e : subs.keySet()) {
            curm = allPivotMSS(subs.get(e));
            for (int i = 0; i < curm.length; i++)
                if (curm[i] * comps.get(e)[i] > ms) {
                    ms = curm[i] * comps.get(e)[i];
                    mseq = subs.get(e).MaximumScoreSubsequence(i);
                    me = e;
                }
        }

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> multCon2Pattern(DynaGraph g) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveConnectedComponent();
        HashMap<Edge, Sequence> subs = g.getAllSubsequences();

        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        int t = 0;

        for (Edge e : subs.keySet()) {
            for (int i = 0; i < g.size(); i++)
                if (comps.get(e)[i] > ms) {
                    ms = comps.get(e)[i];
                    //mseq = subs.get(e).MaximumScoreSubsequence(i);
                    me = e;
                    t = i;
                }
        }

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        mseq = new Subsequence();
        mseq.add(t);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> multCon3Pattern(DynaGraph g) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveConnectedComponent();
        HashMap<Edge, Sequence> subs = g.getAllSubsequences();

        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;

        for (Edge e : subs.keySet()) {
            for (int i = 0; i < g.size(); i++)
                if (comps.get(e)[i] > ms) {
                    ms = comps.get(e)[i];
                    mseq = subs.get(e).MaximumScoreSubsequence(i);
                    me = e;
                }
        }

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> multiHSMSSPattern(DynaGraph g, HSMSSTimes times, QueueEdgesTimeticks queue) throws IOException {

        ElapsedTime et;
        float elapsedTime;
        long before, after, hs, sq, proc, begin, hhs, que;
        begin = System.currentTimeMillis();

        Vector<Graph> scoresHS = new Vector<Graph>();
        scoresHS.setSize(g.size());

        before = System.currentTimeMillis();
        for (int i = 0; i < g.size(); i++) {
            scoresHS.set(i, NetAlgorithms.AllHookedHS(g.get(i), times));
        }
        after = System.currentTimeMillis();
        hs = (after - before);

        et = ElapsedTime.Start();

        before = System.currentTimeMillis();
        HashMap<Edge, Sequence> subs = g.getAllFastSubsequences();
        after = System.currentTimeMillis();
        sq = (after - before);

        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSScopy += elapsedTime;

        et = ElapsedTime.Start();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        float[] curm = null;

        before = System.currentTimeMillis();
        for (Edge e : subs.keySet()) {
            curm = allPivotMSS(subs.get(e));
            for (int i = 0; i < curm.length; i++) {
                queue.set(e, i, scoresHS.get(i).get(e), curm[i]);
            }
        }
        after = System.currentTimeMillis();
        proc = (after - before);

        before = System.currentTimeMillis();
        DynaGraph.EdgeTimetick ed = queue.getMax();
        after = System.currentTimeMillis();
        que = (after - before);


        Subgraph mgraph = NetAlgorithms.HookedHS(g.get(ed.t), ed.ed);

        before = System.currentTimeMillis();
        Sequence seq = g.GetFastInducedSubsequence(mgraph);
        after = System.currentTimeMillis();
        hhs = (after - before);

        mseq = seq.MaximumScoreSubsequence();


        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSS += elapsedTime;

        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        after = System.currentTimeMillis();
        System.out.print(hs + "\t\t" + sq + "\t\t" + proc + "\t\t" + que + "\t\t" + hhs + "\t\t\t" + (after - begin) + "\t");
        return hm;
    }

    public static HashMap<Subsequence, Subgraph> updateHSMSSPattern(DynaGraph g, Subgraph patGraph, Subsequence patSeq, HSMSSTimes times, QueueEdgesTimeticks queue) throws IOException {
        ElapsedTime et;
        float elapsedTime;
        long before, after, hs, sq, proc, begin, hhs, que;
        begin = System.currentTimeMillis();
        Graph scoresHS;

        before = System.currentTimeMillis();
        for (int i : patSeq) {
            scoresHS = NetAlgorithms.AllHookedHS(g.get(i), times);
            for (Edge e : g.get(i).keySet()) {
                queue.changeHS(e, i, scoresHS.get(e));
            }
        }
        after = System.currentTimeMillis();
        hs = (after - before);

        et = ElapsedTime.Start();

        before = System.currentTimeMillis();
        HashMap<Edge, Sequence> subs = g.getAllFastSubsequences();
        after = System.currentTimeMillis();
        sq = (after - before);

        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSScopy += elapsedTime;

        et = ElapsedTime.Start();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        float[] curm = null;

        before = System.currentTimeMillis();
        for (Edge e : patGraph) {
            curm = allPivotMSS(subs.get(e));
            for (int i = 0; i < curm.length; i++) {
                float sc = curm[i];
                queue.changeMSS(e, i, curm[i]);
            }
        }
        after = System.currentTimeMillis();
        proc = (after - before);

        before = System.currentTimeMillis();
        DynaGraph.EdgeTimetick ed = queue.getMax();
        after = System.currentTimeMillis();
        que = (after - before);
        //double score = queue.getValueHS(ed);

        Subgraph mgraph = NetAlgorithms.HookedHS(g.get(ed.t), ed.ed);

        //double score2 = g.get(ed.t).GetSubgraph(mgraph).GetTotalWeight() ;
        //System.out.println(score + "," + score2);

        before = System.currentTimeMillis();
        Sequence seq = g.GetFastInducedSubsequence(mgraph);
        after = System.currentTimeMillis();
        hhs = (after - before);

        mseq = seq.MaximumScoreSubsequence();


        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSS += elapsedTime;

        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        after = System.currentTimeMillis();
        System.out.print(hs + "\t\t" + sq + "\t\t" + proc + "\t\t" + que + "\t\t" + hhs + "\t\t\t" + (after - begin) + "\t");
        return hm;
    }

    public static HashMap<Subsequence, Subgraph> multiHSMSSPattern2(DynaGraph g, HSMSSTimes times, QueueEdgesTimeticks queue) throws IOException {

        ElapsedTime et;
        float elapsedTime;

        Vector<Graph> scoresHS = new Vector<Graph>();
        scoresHS.setSize(g.size());
        for (int i = 0; i < g.size(); i++) {
            scoresHS.set(i, NetAlgorithms.AllHookedHS(g.get(i), times));
        }

        et = ElapsedTime.Start();
        //HashMap<Edge,Sequence> subs = g.getAllSubsequences();
        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSScopy += elapsedTime;

        et = ElapsedTime.Start();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        //float[] curm = null;

        for (Edge e : g.get(0).keySet()) {
            //curm = allPivotMSS(subs.get(e));
            for (int i = 0; i < g.size(); i++) {
                queue.set(e, i, scoresHS.get(i).get(e), 1);
            }
        }

        DynaGraph.EdgeTimetick ed = queue.getMax();

        TopDown mss = new TopDown();
        Subgraph subgraph = mss.Run(g.get(ed.t));
        Sequence sequence = g.GetInducedSubsequence(subgraph);
        mseq = sequence.MaximumScoreSubsequence();
        //mseq = subs.get(ed.ed).MaximumScoreSubsequence(ed.t);
        me = ed.ed;

        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSS += elapsedTime;

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> updateHSMSSPattern2(DynaGraph g, Subgraph patGraph, Subsequence patSeq, HSMSSTimes times, QueueEdgesTimeticks queue) throws IOException {

        ElapsedTime et;
        float elapsedTime;

        Graph scoresHS;
        for (int i : patSeq) {
            scoresHS = NetAlgorithms.AllHookedHS(g.get(i), times);
            for (Edge e : g.get(i).keySet()) {
                queue.changeHS(e, i, scoresHS.get(e));
            }
        }

        et = ElapsedTime.Start();
        HashMap<Edge, Sequence> subs = g.getAllSubsequences();
        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSScopy += elapsedTime;

        et = ElapsedTime.Start();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;


        DynaGraph.EdgeTimetick ed = queue.getMax();
        TopDown mss = new TopDown();
        Subgraph subgraph = mss.Run(g.get(ed.t));
        Sequence sequence = g.GetInducedSubsequence(subgraph);
        mseq = sequence.MaximumScoreSubsequence();
        //mseq = subs.get(ed.ed).MaximumScoreSubsequence(ed.t);
        me = ed.ed;

        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSS += elapsedTime;

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }


    public static HashMap<Subsequence, Subgraph> multiHSMSSPattern3(DynaGraph g, HSMSSTimes times, QueueEdgesTimeticks queue) throws IOException {

        ElapsedTime et;
        float elapsedTime;

        Vector<Graph> scoresHS = new Vector<Graph>();
        scoresHS.setSize(g.size());
        for (int i = 0; i < g.size(); i++) {
            scoresHS.set(i, NetAlgorithms.AllHookedHS(g.get(i), times));
        }

        et = ElapsedTime.Start();
        HashMap<Edge, Sequence> subs = new HashMap<Edge, Sequence>();
        for (Edge e : g.get(0).keySet()) {
            Sequence seq = new Sequence();
            seq.setSize(g.size());
            for (int t = 0; t < g.size(); t++)
                seq.set(t, scoresHS.get(t).get(e));
            subs.put(e, seq);
        }
        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSScopy += elapsedTime;

        et = ElapsedTime.Start();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        float[] curm = null;

        for (Edge e : subs.keySet()) {
            curm = allPivotMSS(subs.get(e));
            for (int i = 0; i < curm.length; i++) {
                queue.set(e, i, scoresHS.get(i).get(e), curm[i]);
            }
        }

        DynaGraph.EdgeTimetick ed = queue.getMaxMSS();
        mseq = subs.get(ed.ed).MaximumScoreSubsequence(ed.t);
        me = ed.ed;

        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSS += elapsedTime;

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> updateHSMSSPattern3(DynaGraph g, Subgraph patGraph, Subsequence patSeq, HSMSSTimes times, QueueEdgesTimeticks queue) throws IOException {

        ElapsedTime et;
        float elapsedTime;

        Graph scoresHS;
        for (int i : patSeq) {
            scoresHS = NetAlgorithms.AllHookedHS(g.get(i), times);
            for (Edge e : g.get(i).keySet()) {
                queue.changeHS(e, i, scoresHS.get(e));
            }
        }

        et = ElapsedTime.Start();
        HashMap<Edge, Sequence> subs = queue.getSequence();
        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSScopy += elapsedTime;

        et = ElapsedTime.Start();
        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        float[] curm = null;

        for (Edge e : g.get(0).keySet()) {
            curm = allPivotMSS(subs.get(e));
            for (int i = 0; i < curm.length; i++) {
                float sc = curm[i];
                queue.changeMSS(e, i, curm[i]);
            }
        }

        DynaGraph.EdgeTimetick ed = queue.getMaxMSS();
        mseq = subs.get(ed.ed).MaximumScoreSubsequence(ed.t);
        me = ed.ed;

        elapsedTime = ElapsedTime.Stop(et);
        times.timeMSS += elapsedTime;

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> maxEdge(DynaGraph g) throws IOException {
        float max = Float.NEGATIVE_INFINITY;
        Edge me = null;
        int t = -1;
        for (int q = 0; q < g.size(); q++) {
            Graph cg = g.get(q);
            for (Edge e : cg.keySet()) {
                if (cg.get(e) > max) {
                    max = cg.get(e);
                    me = e;
                    t = q;
                }
            }
        }

        Subsequence mseq = new Subsequence();
        mseq.add(t);
        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static HashMap<Subsequence, Subgraph> multNeighPattern(DynaGraph g) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveNeighbourScores();
        HashMap<Edge, Sequence> subs = g.getAllSubsequences();

        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        float[] curm = null;

        for (Edge e : subs.keySet()) {
            curm = allPivotMSS(subs.get(e));
            for (int i = 0; i < curm.length; i++)
                if (curm[i] * comps.get(e)[i] > ms) {
                    ms = curm[i] * comps.get(e)[i];
                    mseq = subs.get(e).MaximumScoreSubsequence(i);
                    me = e;
                }
        }

        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);

        return hm;
    }

    public static void neighPattern(DynaGraph g, String outFileGr, String outFileSeq) throws IOException {
        BufferedWriter out1 = new BufferedWriter(new FileWriter(new File(outFileGr)));
        BufferedWriter out2 = new BufferedWriter(new FileWriter(new File(outFileSeq)));
        HashMap<Edge, float[]> comps = g.GetMatrixPositiveNeighbourScores();

        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        HashMap<Subsequence, Float> curm = null;
        int count = 0;
        for (Edge e : comps.keySet()) {
            //System.out.println(count++);
            curm = allMSS(comps.get(e));
            for (Subsequence cs : curm.keySet())
                if (curm.get(cs) > ms) {
                    ms = curm.get(cs);
                    mseq = cs;
                    me = e;
                }
        }

        out1.append(me.source() + "," + me.dest() + "\n");
        // System.out.print("printing ");
        for (int i : mseq) {
            //	System.out.print(i + " ");
            out2.append(i + "\n");
        }
        //System.out.println();
        out1.close();
        out2.close();
    }

    public static HashMap<Subsequence, Subgraph> neighPattern(DynaGraph g) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveNeighbourScores();
        HashMap<Subsequence, Float> curm = null;

        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        //if(seeds == null)

        for (Edge e : comps.keySet()) {
            //System.out.println(count++);
            curm = allMSS(comps.get(e));
            for (Subsequence cs : curm.keySet())
                if (curm.get(cs) > ms) {
                    ms = curm.get(cs);
                    mseq = cs;
                    me = e;
                }
        }
        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);
        return hm;
    }

    public static HashMap<Subsequence, Subgraph> neighPattern(DynaGraph g, HashMap<Integer, ArrayList<Edge>> changes) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveNeighbourScores(changes);
        HashMap<Subsequence, Float> curm = null;

        Edge me = null;
        float ms = Float.NEGATIVE_INFINITY;
        Subsequence mseq = null;
        //if(seeds == null)

        for (Edge e : comps.keySet()) {
            //System.out.println(count++);
            curm = allMSS(comps.get(e));
            for (Subsequence cs : curm.keySet())
                if (curm.get(cs) > ms) {
                    ms = curm.get(cs);
                    mseq = cs;
                    me = e;
                }
        }
        Subgraph mgraph = new Subgraph();
        mgraph.add(me);
        HashMap<Subsequence, Subgraph> hm = new HashMap<Subsequence, Subgraph>();
        hm.put(mseq, mgraph);
        return hm;

	      /*  if (false){//seeds == null) {
              seeds = new TreeMap<Seed, Float> ();
	          for (Edge e:comps.keySet()) {
	          	curm = allMSS(comps.get(e));
	          	for (Subsequence cs:curm.keySet()) {
	          	  Subgraph mgraph = new Subgraph();
	              mgraph.add(e);
	              Seed hm = new Seed();
	              hm.put(cs, mgraph);
	          	  seeds.put(curm.get(cs),hm); 
	          	}
	          }
	        }
	        else {
	          HashSet<Edge> changedEdges = new HashSet<Edge>();
	          for(int t: changes.keySet())
	            for(Edge e: changes.get(t)) {
	              changedEdges.add(e);
	              for(Edge ne : g.edgeLists.get(e))
	                changedEdges.add(ne);
	            }
	             seeds = new TreeSet<SeedEntry> (); 
	          for (Edge e:comps.keySet()) {
	          	curm = allMSS(comps.get(e));
	          	for (Subsequence cs:curm.keySet()) {
	            	  Subgraph mgraph = new Subgraph();
	                mgraph.add(e);
	                Seed hm = new Seed();
	                hm.put(cs, mgraph);
	                seeds.add(new SeedEntry(hm, curm.get(cs)));
	          		}
	          	 Subgraph mgraph = new Subgraph();
	              mgraph.add(e);
	              HashMap<Subsequence,Subgraph> hm = new HashMap<Subsequence,Subgraph>();
	              hm.put(cs, mgraph);
	          	  seeds.put(curm.get(cs),hm); 
	          	
	          }  
	          
	          
	        }
	            
	        return seeds.first().getKey();*/

    }

    public static Subsequence MaximumScoreSubsequence(float[] seq, int pivot) {
        // the pivot is constrained to belong to the Maximum Score Subsequence
        int i = 0;
        float s = 0;
        float sMax = Float.NEGATIVE_INFINITY;
        int iMax = 0, jMax = -1;
        for (int j = 0; (j < seq.length) && (i <= pivot); j++) {
            float w = seq[j];
            s += w;
            if ((w > s + w) && (j <= pivot)) {
                i = j;
                s = w;
            }
            if ((s > sMax) && (i <= pivot) && (pivot <= j)) {
                sMax = s;
                iMax = i;
                jMax = j;
            }
        }
        Subsequence subseq = new Subsequence();
        for (int t = iMax; t <= jMax; t++) {
            subseq.add(t);
        }
        return subseq;
    }

    public static ArrayList<Map.Entry<HashMap<Subsequence, Subgraph>, Float>> allNeighPatterns(DynaGraph g, int limit) throws IOException {

        HashMap<Edge, float[]> comps = g.GetMatrixPositiveNeighbourScores();
        HashMap<HashMap<Subsequence, Subgraph>, Float> rez = new HashMap<HashMap<Subsequence, Subgraph>, Float>();
        HashMap<Subsequence, Subgraph> temp = null;
        Subgraph tempGraph = null;
        float[] curm;

        for (Edge e : comps.keySet()) {
            curm = allPivotMSS(comps.get(e));
            for (int i = 0; i < curm.length; i++) {
                if (curm[i] > limit) {
                    temp = new HashMap<Subsequence, Subgraph>();
                    tempGraph = new Subgraph();
                    tempGraph.add(e);
                    temp.put(MaximumScoreSubsequence(comps.get(e), i), tempGraph);
                    rez.put(temp, curm[i]);
                }
            }
        }

        ArrayList<Map.Entry<HashMap<Subsequence, Subgraph>, Float>> listRez = new ArrayList<Map.Entry<HashMap<Subsequence, Subgraph>, Float>>(rez.entrySet());

        Collections.sort(listRez, new ValueComparator<HashMap<Subsequence, Subgraph>, Float>());
        //System.out.println("--------------" + rez.size() + " " + listRez.get(0).getValue() + " " + listRez.get(100).getValue() + " " + listRez.get(1000).getValue() + " " + listRez.get(10000).getValue() +"------------");
        System.out.println("%-------------" + rez.size());
        return listRez;

    }

    public static boolean checkMaxSubsequence(float[] a, boolean prt) {
        HashMap<Subsequence, Float> mOne = allMSS(a);
        HashMap<Subsequence, Float> mTwo = bruteMSS(a);
        Subsequence one = (Subsequence) mOne.keySet().toArray()[0];
        Subsequence two = (Subsequence) mTwo.keySet().toArray()[0];

        if (prt) {
            System.out.println("allMss  : " + one + " " + mOne.get(one));
            System.out.println("bruteMss: " + two + " " + mTwo.get(two));
        }

        return one.equals(two);
        /*Sequence b = new Sequence();
        for (int i=0; i<a.length; i++)
			b.add(a[i]);
		float[] r = allMSS(a);
				
		for (int i=0; i<a.length; i++){
			
			float q = 0;
			for (int s : b.MaximumScoreSubsequence(i))	{
				//System.out.print(a[s] + " ");
				q += a[s];
			}
			System.out.println(a[i] + " " + r[i] + " " + q);
			if( r[i] != q)
				return false;
		}
		return true;*/
    }

    public static void main(String[] args) {
        float[] a = {257852.0f, 235825.0f, -386057.0f, 34725.0f, 779608.0f, 255519.0f, -843690.0f, -334648.0f, -273929.0f, -265397.0f, -819798.0f, -687592.0f, -604862.0f, -374581.0f, 202123.0f, 150002.0f, 62742.0f, -451688.0f, -277530.0f, 778295.0f, 616357.0f, -318422.0f, 130664.0f, 856311.0f, -340717.0f, 821832.0f, 464999.0f, 912150.0f, 962725.0f, -683084.0f, 399457.0f, 610842.0f, -362677.0f, 516657.0f, 545005.0f, -615018.0f, -164983.0f, 715846.0f, 656695.0f, 759760.0f, -651063.0f, 316447.0f, 870844.0f, -787640.0f, 558132.0f, -833360.0f, -196747.0f, -136238.0f, 317852.0f, 244207.0f, -798916.0f, 548067.0f, 641107.0f, -613205.0f, -163100.0f, 994985.0f, 570028.0f, -452189.0f, 62039.0f, 555186.0f, 714221.0f, -23972.0f, 626380.0f, 941204.0f, 598735.0f, 454535.0f, -128228.0f, 809927.0f, 701308.0f, 850023.0f, -554267.0f, 136356.0f, 864655.0f, 973735.0f, -382781.0f, -546983.0f, 408343.0f, 105324.0f, -12009.0f, -11374.0f, 852371.0f, -147100.0f, -38609.0f, -791928.0f, 792182.0f, -10399.0f, -613816.0f, 668135.0f, 368228.0f, 906529.0f, -99762.0f, -487421.0f, -716659.0f};
        float[] b = {2.4327848f, 11.862417f, 2.8112965f, -0.0588938f, -3.6530747f, 6.610767f, 21.799898f, -0.40178147f, 13.525819f, -0.40178147f, 7.655691f, 18.413778f, 8.58525f, 8.498693f, 14.297746f, -0.2863042f, 4.3912144f, -2.1607733f, 9.0289f, -5.029616f, -3.4512112f, 24.044737f, -0.20573509f, 4.1507387f, 0.4265331f, -1.2515389f, -3.8365014f, -1.3600633f, 3.7813199f, -0.4214638f, 1.3727086f, 11.380615f, 1.9006398f};
        System.out.println(checkMaxSubsequence(b, true));
        printDebug(a);

        Random generator = new Random();
        int count = 0;
        boolean run = true;
        while (run) {
            int n = generator.nextInt(100) + 1;
            float[] rand = new float[n];
            for (int i = 0; i < n; i++) {
                rand[i] = generator.nextInt(1000000000);
                if (generator.nextBoolean())
                    rand[i] *= -1;
            }

            if (!checkMaxSubsequence(rand, false)) {
                System.out.println("FAIL on:");
                printDebug(rand);
                break;
            } else {
                count++;
                if (count % 10000 == 0)
                    System.out.println(count);
            }
        }
    }

    public static void printDebug(float[] rand) {
        int n = rand.length;
        for (int i = 0; i < n; i++)
            System.out.print(i + "\t\t");
        System.out.println();
        for (int i = 0; i < n; i++)
            System.out.print(rand[i] + "f,\t");
    }

    public static boolean checkAllSubsequences(DynaGraph g) {
        HashMap<Edge, Sequence> slow = g.getAllSubsequences();
        HashMap<Edge, Sequence> fast = g.getAllSubsequences();

        return slow.equals(fast);
    }

/*	public static boolean checkQueueMax(QueueEdgesTimeticks queue) {
        DynaGraph.EdgeTimetick slow = queue.getMax();
		DynaGraph.EdgeTimetick fast = queue.getFastMax();
		
		System.out.println(slow + " " + fast);
		return slow.equals(fast);
	} */

    public static boolean checkTotalWeight(DynaGraph dg) {
        for (Graph g : dg) {
            float slow = g.GetTotalWeight();
            float fast = g.GetTotalWeight();
            if (slow != fast)
                return false;
        }

        return true;
    }
}
