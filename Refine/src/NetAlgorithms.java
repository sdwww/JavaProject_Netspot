import exp.RunMesasurements;
import graph.TimeGraph;
import graph.WeightedTimeGraph;

import java.io.IOException;
import java.util.*;


abstract class MaximumScoreSubgraphMethod {
    public abstract Subgraph Run(Graph g) throws Exception;
}

class MaximumScoreSubgraphILPMethod extends MaximumScoreSubgraphMethod {
    public Subgraph Run(Graph g) throws IOException {
        return g.MaximumScoreSubgraph();
    }
}

class TopDown extends MaximumScoreSubgraphMethod {
    public Subgraph Run(Graph g) {
        return NetAlgorithms.MaximumScoreSubgraphTopDown(g);
    }
}

class NetAlgorithms {

    enum Type {
        neigh,
        con,
        multcon,
        multneigh,
        hsmss,
        baseline,
        maxagg,
        max,
        mf,
        multcon2,
        multcon3,
        max2,
        hsmss2, hsmss3, rand,
        baseline_opt, randed, maxed
    }

    ;

    private static Random rand = new Random();

    public static DynaGraph NetworkShuffling(DynaGraph network, int nSlices) {
        DynaGraph newNet = new DynaGraph();
        newNet.setSize(nSlices);

        for (int i = 0; i < nSlices; i++) {
            Graph slice = network.get(0).clone();
            for (Edge e : slice.keySet()) {
                slice.put(e, network.get(rand.nextInt(network.size())).get(e));
            }
            newNet.set(i, slice);
        }
        return newNet;
    }

    public static void Refine(DynaGraph network, Ref<Subgraph> patGraphRef, Ref<Subsequence> patSeqRef, MaximumScoreSubgraphMethod MSS) throws Exception {
        Subgraph patGraph = patGraphRef.get();
        Subsequence patSeq = patSeqRef.get();

        Sequence seq; //= network.GetInducedSubsequence(patGraph);
        //patSeq = seq.MaximumScoreSubsequence();

        float oldscore = 0F;
        float score = network.GetScore(patGraph, patSeq);
        Subsequence patSeqOld;
        do {
            patSeqOld = patSeq;
            oldscore = score;
            Graph g = network.GetInducedSubgraph(patSeq);
            patGraph = MSS.Run(g); // computing the MaximumScoreSubgraph
            seq = network.GetInducedSubsequence(patGraph);
            patSeq = seq.MaximumScoreSubsequence();
            score = network.GetScore(patGraph, patSeq);
        } while (score > oldscore);

        patGraphRef.set(patGraph);
        patSeqRef.set(patSeq);
    }

    public static void RefineAll(DynaGraph network, List<Subgraph> patsGraphOut, List<Subsequence> patsSeqOut, DynaGraph networkClone, MaximumScoreSubgraphMethod MSS, float coverage, int stopNo, float stopScore, int hardLimit, Type type, ArrayList<Double> scores) throws IOException, Exception {

        DynaGraph oldnetwork = network;
        network = networkClone; //network.clone();

        Subgraph patGraph = new Subgraph();
        Subsequence patSeq = new Subsequence();

        patsGraphOut.clear();
        patsSeqOut.clear();

        PatternExtent ext = new PatternExtent(network, patGraph, patSeq);

        Random rnd = new Random();

        HSMSSTimes times = new HSMSSTimes();
        float seedTime = 0f;
        float refineTime = 0f;
        float subtractTime = 0f;
        float cSum = 0;
        int i = 0;
        float area = 0f;
        ElapsedTime et;
        ArrayList<Float> sums = new ArrayList<Float>();
        ArrayList<Float> areas = new ArrayList<Float>();
        float elapsedTime;

        int patternsOver = 0;
        ArrayList<Map.Entry<HashMap<Subsequence, Subgraph>, Float>> patterns = null;
        HashMap<Subsequence, Subgraph> nextPat = null;
        int neighLimit = 80;
        HashMap<Integer, ArrayList<Edge>> changes = new HashMap<Integer, ArrayList<Edge>>();
        if (type == Type.baseline) {
            et = ElapsedTime.Start();
            patterns = GenPattern.allNeighPatterns(networkClone, neighLimit);
            elapsedTime = ElapsedTime.Stop(et);
            seedTime += elapsedTime;
        }
        int nextIndex = 0;
        int count = 0;
        int stopCount = 0;
        QueueEdgesTimeticks queueSeeds = new QueueEdgesTimeticks(network.get(0).size(), network.size());
        Subgraph patGraphLast = null;
        Subsequence patSeqLast = null;
        while (patternsOver < stopNo && stopCount <= hardLimit) {
            stopCount++;

            et = ElapsedTime.Start();
            switch (type) {
                case neigh:
                    nextPat = GenPattern.neighPattern(networkClone, changes);
                    break;
                case con:
                    nextPat = GenPattern.conPattern(networkClone);
                    break;
                case multcon:
                    nextPat = GenPattern.multConPattern(networkClone);
                    break;
                case multneigh:
                    nextPat = GenPattern.multNeighPattern(networkClone);
                    break;
                case hsmss:
                    if (stopCount == 1)
                        nextPat = queueSeeds.GenerateHSMSS(networkClone, times);
                    else
                        nextPat = queueSeeds.UpdateHSMSS(networkClone, patGraphLast, patSeqLast, times);
                    break;
                case baseline:
                    nextPat = patterns.get(nextIndex).getKey();
                    nextIndex++;
                    break;
                case rand:
                    nextPat = GenPattern.randPattern(networkClone, rnd);
                    break;
                case max:
                    nextPat = GenPattern.maxPattern(networkClone);
                    break;
                case max2:
                    nextPat = GenPattern.max2Pattern(networkClone);
                    break;
                case multcon2:
                    nextPat = GenPattern.multCon2Pattern(networkClone);
                    break;
                case multcon3:
                    nextPat = GenPattern.multCon3Pattern(networkClone);
                    break;
                case hsmss2:
                    if (stopCount == 1)
                        nextPat = queueSeeds.GenerateHSMSS2(networkClone, times);
                    else
                        nextPat = queueSeeds.UpdateHSMSS2(networkClone, patGraphLast, patSeqLast, times);
                    break;
                case hsmss3:
                    if (stopCount == 1)
                        nextPat = queueSeeds.GenerateHSMSS3(networkClone, times);
                    else
                        nextPat = queueSeeds.UpdateHSMSS3(networkClone, patGraphLast, patSeqLast, times);
                    break;
                case maxed:
                    nextPat = GenPattern.maxPattern(networkClone);
                    break;
                default:
                    throw new Exception("invalid pattern type, see usage");

            }
            //System.out.println("%" + i + " " + elapsedTime);
            for (Subsequence s : nextPat.keySet())
                patSeq = s;
            patGraph = nextPat.get(patSeq);


            Ref<Subgraph> patGraphRef = new Ref<Subgraph>((Subgraph) patGraph.clone());
            Ref<Subsequence> patSeqRef = new Ref<Subsequence>((Subsequence) patSeq.clone());

            elapsedTime = ElapsedTime.Stop(et);
            seedTime += elapsedTime;

            et = ElapsedTime.Start();
            Refine(network, patGraphRef, patSeqRef, MSS);
            elapsedTime = ElapsedTime.Stop(et);
            refineTime += elapsedTime;

            Subgraph patGraphNew = patGraphRef.get();
            Subsequence patSeqNew = patSeqRef.get();

            //patGraphNew.print();
            //patSeqNew.print();

            //System.out.print("Seed found: " + patGraph.size() + " " + patSeq.size());
            //debugRefine(patGraph,patSeq);
            //System.out.print("Pattern found: " + patGraphNew.size() + " " + patSeqNew.size());
            //debugRefine(patGraphNew,patSeqNew);

            float score = network.GetScore(patGraphNew, patSeqNew);
            System.out.println("%" + i + " " + score);

            if (score < stopScore)
                patternsOver++;
            else {
                patternsOver = 0;
                patsGraphOut.add(patGraphNew);
                patsSeqOut.add(patSeqNew);
                cSum += score;
                area += cSum;
                scores.add((double) score);
                sums.add(cSum);
                areas.add(area);
            }

            et = ElapsedTime.Start();
            ext.Subtract(patGraphNew, patSeqNew);
            changes = network.Erase(patGraphNew, patSeqNew);
            elapsedTime = ElapsedTime.Stop(et);
            subtractTime += elapsedTime;

            if (type == Type.baseline) {
                while (nextIndex < patterns.size()) {
                    for (Subsequence s : patterns.get(nextIndex).getKey().keySet())
                        patSeq = s;
                    patGraph = patterns.get(nextIndex).getKey().get(patSeq);
                    if (Intersects(patGraphNew, patSeqNew, patGraph, patSeq))
                        nextIndex++;
                    else
                        break;
                }
                if (nextIndex > 10) {
                    count++;
                    et = ElapsedTime.Start();
                    patterns = GenPattern.allNeighPatterns(networkClone, Math.max(neighLimit - 10 * count, 10));
                    elapsedTime = ElapsedTime.Stop(et);
                    seedTime += elapsedTime;
                    nextIndex = 0;
                }
            }
            i++;
            patGraphLast = patGraphNew;
            patSeqLast = patSeqNew;
        }

        System.out.println("# Seed_used: " + type);
        System.out.println("# Seed_time: " + seedTime);
        System.out.println("#   MST: " + times.timeMST);
        System.out.println("#   HS: " + times.timeHS);
        System.out.println("#   MSScopy: " + times.timeMSScopy);
        System.out.println("#   MSS: " + times.timeMSS);
        System.out.println("# Refine_time: " + refineTime);
        System.out.println("# Substract_time: " + subtractTime);

        System.out.println("# Subpattern\tScore\tSum_of_Scores\tArea");
        //Collections.sort(scores);
        for (i = 0; i < sums.size(); i++) {
            //System.out.println("" + i + "\t" + scores.get(scores.size() - i - 1) + "\t" + sums.get(i) + "\t" + areas.get(i));
            System.out.println("& " + i + "\t" + scores.get(i) + "\t" + sums.get(i) + "\t" + areas.get(i));
        }
    }

    public static void RefineExperiment(DynaGraph network,
                                        List<Subgraph> patsGraphOut,
                                        List<Subsequence> patsSeqOut,
                                        DynaGraph networkClone,
                                        MaximumScoreSubgraphMethod MSS,
                                        int stopNo,
                                        float stopScore,
                                        int hardLimit, Type type,
                                        ComputeJaccard compJaccard,
                                        ArrayList<Double> time,
                                        ArrayList<Double> scores,
                                        ArrayList<Double> jaccards) throws IOException, Exception {

        DynaGraph oldnetwork = network;
        network = networkClone; //network.clone();

        Subgraph patGraph = new Subgraph();
        Subsequence patSeq = new Subsequence();

        patsGraphOut.clear();
        patsSeqOut.clear();

        //PatternExtent ext = new PatternExtent(network,patGraph,patSeq);

        Random rnd = new Random();

        HSMSSTimes times = new HSMSSTimes();
        float cSum = 0;
        int patternsOver = 0;
        int stopCount = 0;

        HashMap<Subsequence, Subgraph> nextPat = null;
        HashMap<Integer, ArrayList<Edge>> changes = new HashMap<Integer, ArrayList<Edge>>();

        QueueEdgesTimeticks queueSeeds = new QueueEdgesTimeticks(network.size(), network.get(0).size());
        Subgraph patGraphLast = null;
        Subsequence patSeqLast = null;

        long start = System.currentTimeMillis();
        long before, after, hsmss = 0, refine, begin;

        Random rand = new Random();
        System.out.println("AllHookedHS\tAllSubsequences\tallPivotMSS\tqueue.getMax\tGetInducedSubsequence\tHSMSS\tRefine_Time\tTime\tTotal_Time\tScore\t\tSum_of_Scores");
        while (patternsOver < stopNo && stopCount <= hardLimit) {
            begin = System.currentTimeMillis();
            stopCount++;
            switch (type) {
                case neigh:
                    nextPat = GenPattern.neighPattern(networkClone, changes);
                    break;
                case con:
                    nextPat = GenPattern.conPattern(networkClone);
                    break;
                case multcon:
                    nextPat = GenPattern.multConPattern(networkClone);
                    break;
                case multneigh:
                    nextPat = GenPattern.multNeighPattern(networkClone);
                    break;
                case hsmss:
                    if (stopCount == 1) {
                        //queueSeeds = new QueueEdgesTimeticks(network.size(),network.get(0).size());
                        before = System.currentTimeMillis();
                        nextPat = queueSeeds.GenerateHSMSS(networkClone, times);
                        after = System.currentTimeMillis();
                        hsmss = (after - before);
                    } else {
                        before = System.currentTimeMillis();
                        nextPat = queueSeeds.UpdateHSMSS(networkClone, patGraphLast, patSeqLast, times);
                        after = System.currentTimeMillis();
                        hsmss = (after - before);
                    }
                    break;
                case baseline:
                    System.out.print("Error: Baseline should not generate seeds");
                    break;
                case rand:
                    nextPat = GenPattern.randPattern(networkClone, rnd);
                    break;
                case max:
                    nextPat = GenPattern.maxPattern(networkClone);
                    break;
                case max2:
                    nextPat = GenPattern.max2Pattern(networkClone);
                    break;
                case multcon2:
                    nextPat = GenPattern.multCon2Pattern(networkClone);
                    break;
                case multcon3:
                    nextPat = GenPattern.multCon3Pattern(networkClone);
                    break;
                case hsmss2:
                    if (stopCount == 1)
                        nextPat = queueSeeds.GenerateHSMSS2(networkClone, times);
                    else
                        nextPat = queueSeeds.UpdateHSMSS2(networkClone, patGraphLast, patSeqLast, times);
                    break;
                case hsmss3:
                    if (stopCount == 1)
                        //queueSeeds = new QueueEdgesTimeticks(network.size(),network.get(0).size());
                        nextPat = queueSeeds.GenerateHSMSS3(networkClone, times);
                    else
                        nextPat = queueSeeds.UpdateHSMSS3(networkClone, patGraphLast, patSeqLast, times);
                    break;
                case maxed:
                    nextPat = GenPattern.maxPattern(networkClone);
                    break;
                default:
                    throw new Exception("invalid pattern type, see usage");
            }

            for (Subsequence s : nextPat.keySet())
                patSeq = s;
            patGraph = nextPat.get(patSeq);

            Ref<Subgraph> patGraphRef = new Ref<Subgraph>((Subgraph) patGraph.clone());
            Ref<Subsequence> patSeqRef = new Ref<Subsequence>((Subsequence) patSeq.clone());

            // refine

            before = System.currentTimeMillis();
            Refine(network, patGraphRef, patSeqRef, MSS);
            after = System.currentTimeMillis();
            refine = (after - before);

            Subgraph patGraphNew = patGraphRef.get();
            Subsequence patSeqNew = patSeqRef.get();

            float score = network.GetScore(patGraphNew, patSeqNew);

            // update structures
            if (score <= stopScore)
                patternsOver++;
            else {
                patternsOver = 0;
                cSum += score;
            }

            patsGraphOut.add(patGraphNew);
            patsSeqOut.add(patSeqNew);

            network.Erase(patGraphNew, patSeqNew);

            double jaccard = 0.0;
            if (compJaccard != null) {
                compJaccard.add(patGraphNew, patSeqNew);
                jaccard = compJaccard.getJaccard();
            }

            patGraphLast = patGraphNew;
            patSeqLast = patSeqNew;
            // measure time and print
//            System.out.println((System.currentTimeMillis()-start)/1000.0 + "\t" + score + "\t" + cSum + "\t(" + et.ed.source() + "," + et.ed.dest() + ")\t" + et.t + "\t" + value1 + "\t" + valueHS + "\t" + valueMSS);
            after = System.currentTimeMillis();
            System.out.println(refine + "\t\t" + (after - begin) + "\t" + (after - start) / 1000.0 + "\t\t" + score + "\t" + cSum + "\t");
            time.add((System.currentTimeMillis() - start) / 1000.0);
            scores.add((double) score);
            jaccards.add(jaccard);
        }
    }


    public static boolean Intersects(Subgraph g1, Subsequence s1, Subgraph g2, Subsequence s2) {
        boolean time = false;
        for (int i : s1)
            if (s2.contains(i)) {
                time = true;
                break;
            }
        if (time == false)
            return false;
        for (Edge e : g1)
            if (g2.contains(e))
                return true;
        return false;
    }

    public static void RefineAll(DynaGraph network, Subgraph patGraph, Subsequence patSeq, List<Subgraph> patsGraphOut, List<Subsequence> patsSeqOut, DynaGraph networkClone, MaximumScoreSubgraphMethod MSS, float coverage, int noSamples, float stopScore) throws IOException, Exception {
        DynaGraph oldnetwork = network;
        network = networkClone; //network.clone();

        patsGraphOut.clear();
        patsSeqOut.clear();

        PatternExtent ext = new PatternExtent(network, patGraph, patSeq);
        float energy = ext.GetEnergy();
        float residual = energy;

        // System.out.println("Subpattern\tScore\tEnergy\tConsumed\tResidual");
        //System.out.println("0\t\t\t\t" + energy); // initial energy

        boolean bExit = false;

        int i = 0;
        int zeroConsumed = 0;
        float area = 0f;
        ElapsedTime et;
        ArrayList<Float> sums = new ArrayList<Float>();
        ArrayList<Float> scores = new ArrayList<Float>();
        ArrayList<Float> areas = new ArrayList<Float>();
        float elapsedTime;
        float cSum = 0f;
        float score = 0f;
        float refineTime = 0f;
        stopScore *= noSamples;
        float curScore = Float.POSITIVE_INFINITY;
        int stopCount = 0;
        //while ((!ext.isEmpty()) && (residual > energy * (1.0-coverage)) && (bExit==false)) {
        while (curScore >= stopScore && stopCount <= 100) {
            stopCount++;
            System.out.println("%" + i);

            Ref<Subgraph> patGraphRef = new Ref<Subgraph>((Subgraph) patGraph.clone());
            Ref<Subsequence> patSeqRef = new Ref<Subsequence>((Subsequence) patSeq.clone());

            et = ElapsedTime.Start();
            Refine(network, patGraphRef, patSeqRef, MSS);
            elapsedTime = ElapsedTime.Stop(et);
            refineTime += elapsedTime;

            Subgraph patGraphNew = patGraphRef.get();
            Subsequence patSeqNew = patSeqRef.get();

//                patGraphNew.print();
//                patSeqNew.print();
//                System.out.println("Score: " + GetScore(patGraphNew, patSeqNew));

            patsGraphOut.add(patGraphNew);
            patsSeqOut.add(patSeqNew);

            float consumed = ext.GetEnergy(patGraphNew, patSeqNew);
            if (consumed == 0F)
                zeroConsumed++;
            else
                zeroConsumed = 0;

            score = oldnetwork.GetScore(patGraphNew, patSeqNew);
            cSum += score;
            area += cSum;

            scores.add(score);
            sums.add(cSum);
            areas.add(area);
            //System.out.print("" + i + "\t" + score + "\t" + network.GetEnergy(patGraphNew,patSeqNew) + "\t" + consumed + "\t");

            if (i == noSamples - 1) {
                curScore = 0;
                for (int s = 0; s <= i; s++)
                    curScore += scores.get(s);
            } else if (i >= noSamples)
                curScore = curScore + scores.get(i) - scores.get(i - noSamples);

            ext.Subtract(patGraphNew, patSeqNew);
            network.Erase(patGraphNew, patSeqNew);

            //System.out.println("" + ext.GetEnergy());

            if (patGraphNew.size() == 0) {
                System.out.println("# No more sub-patterns");
                bExit = true;
            }
            if (zeroConsumed >= 5) {
                System.out.println("# No more energy to consume");
                bExit = true;
            }
            i++;

        }

        System.out.println("# Refine_time: " + refineTime);

        System.out.println("# Subpattern\tScore\tSum_of_Scores\tArea");
        //Collections.sort(scores);
        for (i = 0; i < sums.size(); i++) {
            //System.out.println("" + i + "\t" + scores.get(scores.size() - i - 1) + "\t" + sums.get(i) + "\t" + areas.get(i));
            System.out.println("" + i + "\t" + scores.get(i) + "\t" + sums.get(i) + "\t" + areas.get(i));
        }
    }

    public static void RefineAll(DynaGraph network, Subgraph patGraph, Subsequence patSeq, List<Subgraph> patsGraphOut, List<Subsequence> patsSeqOut, DynaGraph networkClone, MaximumScoreSubgraphMethod MSS, float coverage) throws IOException, Exception {
        DynaGraph oldnetwork = network;
        network = networkClone; //network.clone();

        patsGraphOut.clear();
        patsSeqOut.clear();

        PatternExtent ext = new PatternExtent(network, patGraph, patSeq);
        float energy = ext.GetEnergy();
        float residual = energy;

        System.out.println("Subpattern\tScore\tEnergy\tConsumed\tResidual");
        System.out.println("0\t\t\t\t" + energy); // initial energy

        boolean bExit = false;

        int i = 1;
        int zeroConsumed = 0;
        while ((!ext.isEmpty()) && (residual > energy * (1.0 - coverage)) && (bExit == false)) {
            Ref<Subgraph> patGraphRef = new Ref<Subgraph>((Subgraph) patGraph.clone());
            Ref<Subsequence> patSeqRef = new Ref<Subsequence>((Subsequence) patSeq.clone());
            Refine(network, patGraphRef, patSeqRef, MSS);
            Subgraph patGraphNew = patGraphRef.get();
            Subsequence patSeqNew = patSeqRef.get();

//                patGraphNew.print();
//                patSeqNew.print();
//                System.out.println("Score: " + GetScore(patGraphNew, patSeqNew));

            patsGraphOut.add(patGraphNew);
            patsSeqOut.add(patSeqNew);

            float consumed = ext.GetEnergy(patGraphNew, patSeqNew);
            if (consumed == 0F)
                zeroConsumed++;
            else
                zeroConsumed = 0;
            System.out.print("" + i + "\t" + oldnetwork.GetScore(patGraphNew, patSeqNew) + "\t" + network.GetEnergy(patGraphNew, patSeqNew) + "\t" + consumed + "\t");

            ext.Subtract(patGraphNew, patSeqNew);
            network.Erase(patGraphNew, patSeqNew);

            System.out.println("" + ext.GetEnergy());

            if (patGraphNew.size() == 0) {
                System.out.println("No more sub-patterns");
                bExit = true;
            }
            if (zeroConsumed >= 5) {
                System.out.println("No more energy to consume");
                bExit = true;
            }
            i++;
        }
    }

    public static void TopKMeden(DynaGraph network, List<Subgraph> patsGraphOut, List<Subsequence> patsSeqOut, DynaGraph networkClone, float threshold, boolean exact, ArrayList<Double> times, ArrayList<Double> scores, int hardLimit) throws IOException, Exception {
        DynaGraph oldnetwork = network;
        network = networkClone; //network.clone();

        patsGraphOut.clear();
        patsSeqOut.clear();

        //System.out.println("Subpattern\tScore");
        System.out.println("Time\tScore\tAggScore");
        float sumScore = 0;
        long start = System.currentTimeMillis();
        Pattern pat = new Pattern();
        int i = 1;
        int count = 0;
        do {

            // TODO(Razvan,Misael): this is the call to a weighted a graph
            // top1 in Meden, if you need to run non-weighted experiments use the old function (commented below)
            //MaximumScoreDynamicSubgraph (network, pat, exact);
            count++;
            MaximumScoreDynamicSubgraphWeighted(network, pat, exact);
            Subgraph patGraphNew = pat.subgraph;
            Subsequence patSeqNew = pat.subsequence;

            //patGraphNew.print();
            //patSeqNew.print();

            patsGraphOut.add(patGraphNew);
            patsSeqOut.add(patSeqNew);

            if (pat.score >= threshold) sumScore += pat.score;
            double t = (System.currentTimeMillis() - start) / 1000.0;
            System.out.println("%" + i + "\t" + t + "\t" +
                    pat.score + "\t" +
                    sumScore);
            scores.add(pat.score);
            times.add(t);
            //System.out.println("" + i + "\t" + pat.score + "\t" + network.GetScore(patGraphNew,patSeqNew));

            network.Erase(patGraphNew, patSeqNew);
            i++;
        } while (pat.score > threshold && count < hardLimit);
    }

    public static void RefineAll(DynaGraph network, Subgraph patGraph, Subsequence patSeq, List<Subgraph> patsGraphOut, List<Subsequence> patsSeqOut) throws IOException {
        PatternExtent ext = new PatternExtent(network, patGraph, patSeq);

        while (!ext.isEmpty()) {
            PatternExtent.Entry entry = ext.PickEntry();

            Ref<Subgraph> patGraphRef = new Ref<Subgraph>((Subgraph) patGraph.clone());
            Ref<Subsequence> patSeqRef = new Ref<Subsequence>((Subsequence) patSeq.clone());
            Refine(network, patGraphRef, patSeqRef, entry.edge, entry.time);
            Subgraph patGraphNew = patGraphRef.get();
            Subsequence patSeqNew = patSeqRef.get();

            //patGraphNew.print();
            //patSeqNew.print();

            patsGraphOut.add(patGraphNew);
            patsSeqOut.add(patSeqNew);

            ext.Subtract(patGraphNew, patSeqNew);
            ext.remove(entry);

            if (patGraphNew.size() == 0)
                System.out.println("failed");
            System.out.println("Pattern extension: " + ext.size());
        }
    }

    public static void Refine(DynaGraph network, Subgraph patGraph, Subsequence patSeq, List<Subgraph> patsGraphOut, List<Subsequence> patsSeqOut) throws IOException {
        PatternExtent ext = new PatternExtent(network, patGraph, patSeq);

        while (!ext.isEmpty()) {
            PatternExtent.Entry entry = ext.PickEntry();

            Ref<Subgraph> patGraphRef = new Ref<Subgraph>((Subgraph) patGraph.clone());
            Ref<Subsequence> patSeqRef = new Ref<Subsequence>((Subsequence) patSeq.clone());
            Refine(network, patGraphRef, patSeqRef, entry.edge, entry.time);
            Subgraph patGraphNew = patGraphRef.get();
            Subsequence patSeqNew = patSeqRef.get();

//                patGraphNew.print();
//                patSeqNew.print();
//                System.out.println("Score: " + GetScore(patGraphNew, patSeqNew));

            patsGraphOut.add(patGraphNew);
            patsSeqOut.add(patSeqNew);

            ext.Subtract(patGraphNew, patSeqNew);
            ext.remove(entry);

            if (patGraphNew.size() == 0)
                System.out.println("failed");
            System.out.println("Pattern extension: " + ext.size());
        }
    }

    public static void Refine(DynaGraph network, Ref<Subgraph> patGraphRef, Ref<Subsequence> patSeqRef, Edge edgePivot, int timePivot) throws IOException {
        Subgraph patGraph = patGraphRef.get();
        Subsequence patSeq = patSeqRef.get();

        Sequence seq = network.GetInducedSubsequence(patGraph);
        patSeq = seq.MaximumScoreSubsequence(timePivot);
//        System.out.println("Score: " + GetScore(patGraph, patSeq));

        Subsequence patSeqOld;
        do {
            patSeqOld = patSeq;
            Graph g = network.GetInducedSubgraph(patSeq);
            patGraph = g.MaximumScoreSubgraph(edgePivot);
//            System.out.println("Score: " + GetScore(patGraph, patSeq));
            seq = network.GetInducedSubsequence(patGraph);
            patSeq = seq.MaximumScoreSubsequence(timePivot);
//            System.out.println("Score: " + GetScore(patGraph, patSeq));
        } while (network.GetScore(patGraph, patSeq) != network.GetScore(patGraph, patSeqOld));

        patGraphRef.set(patGraph);
        patSeqRef.set(patSeq);
    }

    public static void MaximumScoreDynamicSubgraphWeighted(DynaGraph network,
                                                           Pattern outPattern,
                                                           boolean exact) throws Exception {
        // first get the structure
        graph.Graph g = getMedenGraph(network.get(0));
        // now fill in the slices
        float[][] slices = new float[network.size()][g.getm()];
        for (int t = 0; t < slices.length; t++) {
            for (Edge e : network.get(t).keySet()) {
                slices[t][g.getEdgeIndex(e.source(), e.dest())] = network.get(t).get(e);
                slices[t][g.getEdgeIndex(e.dest(), e.source())] = network.get(t).get(e);
            }
        }
        WeightedTimeGraph tg = new WeightedTimeGraph(g, slices);
        // meas will collect statistics of the run (times, pruning, etc.)
        exp.RunMesasurements meas = new RunMesasurements();
        // the call to meden for the whole time span, alpha = 0.4, 10 top intervals for lowerbound
        graph.Pattern pt = exp.BoundUtility.runCoarse(tg, 0, network.size() - 1,
                network.size(), 0.4, 10, meas, exact);
        // now convert to the other pattern
        //outPattern = new Pattern();
        // score
        outPattern.score = pt.score;
        // interval
        outPattern.subsequence = new Subsequence();
        for (int i = pt.s; i <= pt.e; i++) outPattern.subsequence.add(i);
        outPattern.subgraph = new Subgraph();
        for (Edge e : network.get(0).keySet()) {
            if (pt.edges.get(tg.getEdgeIndex(e.source(), e.dest()))) {
                outPattern.subgraph.add(e);
            }
        }
    }

    public static void MaximumScoreDynamicSubgraph(DynaGraph network, Pattern outPattern, boolean exact) throws Exception {
        // first get the structure
        graph.Graph g = getMedenGraph(network.get(0));
        // now fill in the slices
        BitSet[] slices = new BitSet[network.size()];
        for (int t = 0; t < slices.length; t++) {
            slices[t] = new BitSet(g.getm());
            for (Edge e : network.get(t).keySet()) {
                if (network.get(t).get(e) > 0.0) {
                    slices[t].set(g.getEdgeIndex(e.source(), e.dest()));
                    slices[t].set(g.getEdgeIndex(e.dest(), e.source()));
                }
            }
        }
        TimeGraph tg = new TimeGraph(g, slices);
        // meas will collect statistics of the run (times, pruning, etc.)
        exp.RunMesasurements meas = new RunMesasurements();
        // the call to meden for the whole time span, alpha = 0.4, 10 top intervals for lowerbound
        graph.Pattern pt = exp.BoundUtility.runCoarse(tg, 0, network.size() - 1,
                network.size(), 0.4, 10, meas, exact);
        // now convert to the other pattern
        //outPattern = new Pattern();
        // score
        outPattern.score = pt.score;
        // interval
        outPattern.subsequence = new Subsequence();
        for (int i = pt.s; i <= pt.e; i++) outPattern.subsequence.add(i);
        outPattern.subgraph = new Subgraph();
        for (Edge e : network.get(0).keySet()) {
            if (pt.edges.get(tg.getEdgeIndex(e.source(), e.dest()))) {
                outPattern.subgraph.add(e);
            }
        }
    }

    // run Meden's TopDown
    public static Subgraph MaximumScoreSubgraphTopDown(Graph graph) {
        graph.Graph g = getMedenGraph(graph);
        BitSet sel = new BitSet(g.getm());
        index.IntervalBounds.topDownLBFast(g, sel);
        Subgraph sg = new Subgraph();
        for (Edge e : graph.keySet()) {
            if (sel.get(g.getEdgeIndex(e.source(), e.dest()))) {
                sg.add(e);
            }
        }
        return sg;
    }

    // Builds a meden graph from this
    private static graph.Graph getMedenGraph(Graph graph) {
        // First get number of nodes
        // Please, verify an if it is not the case, we have to discuss
        int n = 0;
        HashMap<Integer, Integer> deg = new HashMap<Integer, Integer>();
        for (Edge e : graph.keySet()) {
            n = Math.max(n, e.dest()); // this line is expecting that s < d, for all edges
            if (!deg.containsKey(e.source())) deg.put(e.source(), 0);
            deg.put(e.source(), deg.get(e.source()) + 1);
            if (!deg.containsKey(e.dest())) deg.put(e.dest(), 0);
            deg.put(e.dest(), deg.get(e.dest()) + 1);
        }
        // n is the biggest index we saw in the edges now increment to have the node count
        n = n + 1;
        graph.Graph g = new graph.Graph(n, graph.size() * 2);
        // fill in the degree
        g.ind[0] = 0;
        for (int i = 0; i < n; i++) {
            assert (deg.containsKey(i));
            g.ind[i + 1] = g.ind[i] + deg.get(i);
            g.names[i] = "" + i;
        }
        // fill in edges and weights
        Arrays.fill(g.endv, Integer.MIN_VALUE);
        Arrays.fill(g.we, Double.NEGATIVE_INFINITY);
        int eind = -1;
        for (Edge e : graph.keySet()) {
            // add e.d as a neighbor of e.s
            eind = g.ind[e.source()];
            while (g.endv[eind] > Integer.MIN_VALUE) eind++;
            assert (eind < g.ind[e.source() + 1]);
            g.endv[eind] = e.dest();
            g.we[eind] = graph.get(e);

            // add e.s as a neighbor of e.d
            eind = g.ind[e.dest()];
            while (g.endv[eind] > Integer.MIN_VALUE) eind++;
            assert (eind < g.ind[e.dest() + 1]);
            g.endv[eind] = e.source();
            g.we[eind] = graph.get(e);
        }
        return g;
    }

    // this is a unit test for running MaximumScoreSubgraphTopDown
    private static boolean testMedenTopDown() {
        // create a graph
        Graph g = new Graph();
        Edge e;
        e = new Edge(0, 1);
        g.put(e, (float) 2.0);
        e = new Edge(1, 2);
        g.put(e, (float) -1.0);
        e = new Edge(1, 3);
        g.put(e, (float) -1.0);
        e = new Edge(1, 4);
        g.put(e, (float) -8.0);
        e = new Edge(2, 3);
        g.put(e, (float) -2.0);
        e = new Edge(3, 4);
        g.put(e, (float) 2.0);
        e = new Edge(3, 5);
        g.put(e, (float) 1.0);
        e = new Edge(4, 5);
        g.put(e, (float) 3.0);

        // create the same graph Meden format
        graph.Graph gm = new graph.Graph(6, 16);
        String[] names = {"0", "1", "2", "3", "4", "5"};
        gm.names = names;
        int[] ind = {0, 1, 5, 7, 11, 14, 16};
        gm.ind = ind;
        int[] endv = {1, 0, 2, 3, 4, 1, 3, 1, 2, 4, 5, 1, 3, 5, 3, 4};
        gm.endv = endv;
        double[] we = {2.0, 2.0, -1.0, -1.0, -8.0, -1.0, -2.0, -1.0, -2.0, 2.0, 1.0, -8.0, 2.0, 3.0, 1.0, 3.0};
        gm.we = we;

        graph.Graph gm_conv = getMedenGraph(g);

        // Make sure gm_conv and g are the same
        assert (gm.getn() == gm_conv.getn()) : "Node counts do not match";
        assert (gm.getm() == gm_conv.getm()) : "Edge counts do not match";

        for (int i = 0; i < gm.getn(); i++) {
            for (int j = gm.ind[i]; j < gm.ind[i + 1]; j++) {
                int idx = gm_conv.getEdgeIndex(i, gm.endv[j]);
                assert (idx >= 0) : "Edge " + i + " - " + gm.endv[j] + " does not exist in constructed graph";
                assert (gm.we[j] == gm_conv.we[idx]) : "Edge " + i + " - " + gm.endv[j] + " has different weight in constructed graph";
            }
        }

        // verify the topDown invocation
        Subgraph sg = MaximumScoreSubgraphTopDown(g);
        assert (sg.size() == 5) : "top down did not find the right number of edges";
        for (Edge ed : sg) {
            assert ((ed.source() == 0 && ed.dest() == 1) ||
                    (ed.source() == 1 && ed.dest() == 3) ||
                    (ed.source() == 3 && ed.dest() == 4) ||
                    (ed.source() == 3 && ed.dest() == 5) ||
                    (ed.source() == 4 && ed.dest() == 5)
            );
        }
        return true;
    }

    public static Subgraph HookedHS(Graph graph, Edge e) {
        HSMSSTimes times = new HSMSSTimes();

        ElapsedTime et;
        float elapsedTime;

        GraphFast g = getFastGraph(graph);
        BitSet edges = new BitSet(g.getm());

        if (g.getn() == 1) return null;
        Arrays.fill(g.wn, 0);

        TreeSet<GraphFast.ComparableEdge> edg = new TreeSet<GraphFast.ComparableEdge>();
        // Order edges by increasing cost
        for (int i = 0; i < g.getn(); i++) {
            for (int j = g.ind[i]; j < g.ind[i + 1]; j++) {
                if (i < g.endv[j]) {
                    edg.add(new GraphFast.ComparableEdge(g.we[j], j, g.getReverseEdgeIndex(j), i, g.endv[j]));
                }
            }
        }

        et = ElapsedTime.Start();
        TreeFast t = computeMaxSTFast(g, edg);
        t.buildN2I();

        HashMap<Integer, HashSet<Edge>> aggEdges = new HashMap<Integer, HashSet<Edge>>();
        for (int i = 0; i < g.getn(); i++) {
            aggEdges.put(i, new HashSet<Edge>());
        }

        // distribute positive edge mass to adjacent nodes
        // the doubled edge score here is skipped
        for (GraphFast.ComparableEdge ce : edg.descendingSet()) {
            if (ce.w > 0) {
                int s = Integer.parseInt(g.names[ce.v1]);
                int d = Integer.parseInt(g.names[ce.v2]);
                if (t.n2i.containsKey(g.names[ce.v1])) {
                    aggEdges.get(s).add(new Edge(Math.min(s, d), Math.max(s, d)));
                    t.wn[t.n2i.get(g.names[ce.v1])] += ce.w;
                } else if (t.n2i.containsKey(g.names[ce.v2])) {
                    aggEdges.get(d).add(new Edge(Math.min(s, d), Math.max(s, d)));
                    t.wn[t.n2i.get(g.names[ce.v2])] += ce.w;
                }
            } else {
                break;
            }
        }

        // keep only negative edges in the tree,
        // take their absolute value
        for (int i = 0; i < t.getm(); i++) {
            if (t.we[i] >= 0) t.we[i] = 0;
            else t.we[i] = (-1) * t.we[i];
        }
        elapsedTime = ElapsedTime.Stop(et);
        times.timeMST += elapsedTime;

        HSTree tree = new HSTree();
        tree.build(t, aggEdges);

        et = ElapsedTime.Start();
        Subgraph subgraph = tree.computeHS(e);

        elapsedTime = ElapsedTime.Stop(et);
        times.timeHS += elapsedTime;


/*        // compute connected component of score == score of e
        Subgraph subgraph = new Subgraph();
        g.getLocalEdgeLists();
        HashMap<Edge, ArrayList<Edge>> graphList = g.localEdgeLists;
        BFS(g,graphList,e,Math.max(g.get(e),0),subgraph);
*/
        return subgraph;
    }

    public static Graph computeAllOptTreeScores2(TreeFast t, Graph graph, HashMap<Integer, HashSet<Edge>> aggEdges, HSMSSTimes times) {
        HSTree tree = new HSTree();
        tree.build(t, aggEdges);

        ElapsedTime et = ElapsedTime.Start();
        Graph gscore = tree.computeHSall(graph);
        float elapsedTime = ElapsedTime.Stop(et);
        times.timeHS += elapsedTime;

        return gscore;
    }

    private static void BFS(Graph graph, HashMap<Edge, ArrayList<Edge>> graphList, Edge e, double threshold, Subgraph exploredSubgraph) {
        if (exploredSubgraph.contains(e))
            return;
        exploredSubgraph.add(e);
        for (Edge n : graphList.get(e))
            if (graph.get(n) >= threshold)
                BFS(graph, graphList, n, threshold, exploredSubgraph);
    }

    // use this for edge-weighted fast version
    public static Graph AllHookedHS(Graph graph, HSMSSTimes times) {

        ElapsedTime et;
        float elapsedTime;
        long before, after, st, begin, opt, ed, cmp;
        begin = System.currentTimeMillis();

        GraphFast g = getFastGraph(graph);
        BitSet edges = new BitSet(g.getm());

        if (g.getn() == 1) return null;
        Arrays.fill(g.wn, 0);


        TreeSet<GraphFast.ComparableEdge> edg = new TreeSet<GraphFast.ComparableEdge>();

        before = System.currentTimeMillis();
        // Order edges by increasing cost
        for (int i = 0; i < g.getn(); i++) {
            for (int j = g.ind[i]; j < g.ind[i + 1]; j++) {
                if (i < g.endv[j]) {
                    edg.add(new GraphFast.ComparableEdge(g.we[j], j, g.getReverseEdgeIndex(j), i, g.endv[j]));
                }
            }
        }
        after = System.currentTimeMillis();
        ed = after - before;

        before = System.currentTimeMillis();
        et = ElapsedTime.Start();
        TreeFast t = computeMaxSTFast(g, edg);
        after = System.currentTimeMillis();
        st = after - before;

        t.buildN2I();

        HashMap<Integer, HashSet<Edge>> aggEdges = new HashMap<Integer, HashSet<Edge>>();
        for (int i = 0; i < g.getn(); i++) {
            aggEdges.put(i, new HashSet<Edge>());
        }


        // distribute positive edge mass to adjacent nodes
        // the doubled edge score here is skipped

        for (GraphFast.ComparableEdge ce : edg.descendingSet()) {
            if (ce.w > 0) {
                int s = Integer.parseInt(g.names[ce.v1]);
                int d = Integer.parseInt(g.names[ce.v2]);
                if (t.n2i.containsKey(g.names[ce.v1])) {
                    aggEdges.get(s).add(new Edge(Math.min(s, d), Math.max(s, d)));
                    t.wn[t.n2i.get(g.names[ce.v1])] += ce.w;
                } else if (t.n2i.containsKey(g.names[ce.v2])) {
                    aggEdges.get(d).add(new Edge(Math.min(s, d), Math.max(s, d)));
                    t.wn[t.n2i.get(g.names[ce.v2])] += ce.w;
                }
            } else {
                break;
            }
        }


        // keep only negative edges in the tree,
        // take their absolute value
        for (int i = 0; i < t.getm(); i++) {
            if (t.we[i] >= 0) t.we[i] = 0;
            else t.we[i] = (-1) * t.we[i];
        }


        elapsedTime = ElapsedTime.Stop(et);
        times.timeMST += elapsedTime;

        before = System.currentTimeMillis();
        Graph gnew = computeAllOptTreeScores2(t, graph, aggEdges, times);
        after = System.currentTimeMillis();
        opt = after - before;

        after = System.currentTimeMillis();

        //System.out.print((after - begin) + "-" + st +"-" + opt + "-" + ed + " " );
        return gnew;
    }

    // use this for edge-weighted fast version
    /*public static TreeFast topDownLBFast(GraphFast g, BitSet edges) {
        if(g.getn()==1) return null;
		Arrays.fill(g.wn, 0);

		long start = System.currentTimeMillis();
		TreeSet<ComparableEdge> edg = new TreeSet<ComparableEdge>();
		// Order edges by increasing cost
		for (int i = 0; i < g.getn(); i++) {
			for( int j = g.ind[i]; j < g.ind[i+1]; j++ ) {
				if (i < g.endv[j]) {
					edg.add(new ComparableEdge(g.we[j], j, g.getReverseEdgeIndex(j), i, g.endv[j]));
				}
			}
		}
		long ordert = System.currentTimeMillis() - start;

		start = System.currentTimeMillis();
		TreeFast t = ST.computeMaxSTFast(g, edg);t.buildN2I();
		long mstt = System.currentTimeMillis() - start;


		// distribute positive edge mass to adjacent nodes
		// the doubled edge score here is skipped
		for (ComparableEdge ce : edg.descendingSet()) {
			if (ce.w > 0) {
				if(t.n2i.containsKey(g.names[ce.v1])){
					t.wn[t.n2i.get(g.names[ce.v1])] += ce.w;
				}else if (t.n2i.containsKey(g.names[ce.v2])){
					t.wn[t.n2i.get(g.names[ce.v2])] += ce.w;
				}
			} else {
				break;
			}
		}

		// keep only negative edges in the tree,
		// take their absolute value
		for (int i=0; i < t.getm(); i++) {
			if(t.we[i] >= 0) t.we[i] = 0;
			else t.we[i] = (-1)*t.we[i];
		}


		start = System.currentTimeMillis();
		// TreeFast PCST
		TreeFast ot = computeOptTreeNW(t);
		long pcstt = System.currentTimeMillis() - start;

		start = System.currentTimeMillis();
		if (null != edges) {
			// first add connector edges from the spanning tree
			ArrayList<Integer> el = g.getSubgraphEdgesStructureOnly(ot);
			for (int e:el) edges.set(e);
			// now add all positive edges adjacent to included nodes
			HashSet<String> nds = new HashSet<String>();
			for(String n:ot.names) nds.add(n);
			for (int i=0; i < g.getn(); i++){
				if(!nds.contains(g.names[i])) continue;
				for (int j=g.ind[i]; j < g.ind[i+1]; j++) {
					if (g.we[j] >=0)
						edges.set(j);
				}
			}
		}
		long reconstructt = System.currentTimeMillis() - start;


		//ExtToolsDriver.showInGoogleEarth(ot.names);
		start = System.currentTimeMillis();
		double s =  ot.getScore();
		long scoret = System.currentTimeMillis() - start;

//		System.err.print("" + s + "\tmst: " + mstt  +
//						 "ms\tpcst: " + pcstt +
//						 "ms\tscore: " + scoret + "ms\n");
		return ot;
	}*/

    public static TreeFast computeMaxSTFast(GraphFast g, TreeSet<GraphFast.ComparableEdge> edges) {
        long start = System.currentTimeMillis();
        TreeSet<GraphFast.ComparableEdge> tmpedges = new TreeSet<GraphFast.ComparableEdge>();
        int[] comp = new int[g.getn()];
        for (int i = 0; i < g.getn(); i++) comp[i] = i;

        start = System.currentTimeMillis();
        // Get the MST
        int c2 = -1;
        double sum = 0;
        for (GraphFast.ComparableEdge e : edges.descendingSet()) {
            if (comp[e.v1] != comp[e.v2]) {
                // add the edge indices
                tmpedges.add(e);
                c2 = comp[e.v2];
                for (int i = 0; i < comp.length; i++) {
                    if (comp[i] == c2) {
                        comp[i] = comp[e.v1];
                    }
                }
                if (tmpedges.size() >= (g.getn() - 1)) break;
            }
        }
        long kruskalt = System.currentTimeMillis() - start;

        start = System.currentTimeMillis();
        TreeFast t = new TreeFast(g.subgraph(tmpedges));
        long constt = System.currentTimeMillis() - start;

        //System.err.print(sum + "( ord: " + ordert + " + krus: " + kruskalt + " + constr: " + constt + ")" );
        return t;
    }

    public static GraphFast getFastGraph(Graph graph) {
        // First get number of nodes
        // Please, verify an if it is not the case, we have to discuss
        int n = 0;
        HashMap<Integer, Integer> deg = new HashMap<Integer, Integer>();
        for (Edge e : graph.keySet()) {
            n = Math.max(n, e.dest()); // this line is expecting that s < d, for all edges
            assert e.source() < e.dest();
            if (!deg.containsKey(e.source())) deg.put(e.source(), 0);
            deg.put(e.source(), deg.get(e.source()) + 1);
            if (!deg.containsKey(e.dest())) deg.put(e.dest(), 0);
            deg.put(e.dest(), deg.get(e.dest()) + 1);
        }
        // n is the biggest index we saw in the edges now increment to have the node count
        n = n + 1;
        GraphFast g = new GraphFast(n, graph.size() * 2);

        // Some nodes might be isolated add them with degree 0
        for (int i = 0; i < n; i++) {
            if (!deg.containsKey(i)) deg.put(i, 0);
        }

        // fill in the degree
        g.ind[0] = 0;
        for (int i = 0; i < n; i++) {
            assert (deg.containsKey(i));
            g.ind[i + 1] = g.ind[i] + deg.get(i);
            g.names[i] = "" + i;
        }
        // fill in edges and weights
        Arrays.fill(g.endv, Integer.MIN_VALUE);
        Arrays.fill(g.we, Double.NEGATIVE_INFINITY);
        int eind = -1;
        for (Edge e : graph.keySet()) {
            // add e.d as a neighbor of e.s
            eind = g.ind[e.source()];
            while (g.endv[eind] > Integer.MIN_VALUE) eind++;
            assert (eind < g.ind[e.source() + 1]);
            g.endv[eind] = e.dest();
            g.we[eind] = graph.get(e);

            // add e.s as a neighbor of e.d
            eind = g.ind[e.dest()];
            while (g.endv[eind] > Integer.MIN_VALUE) eind++;
            assert (eind < g.ind[e.dest() + 1]);
            g.endv[eind] = e.source();
            g.we[eind] = graph.get(e);
        }
        return g;
    }


    public static Graph computeAllOptTreeScores(TreeFast t, Graph graph) {
        HashMap<Integer, Double> vals = new HashMap<Integer, Double>();
        HashMap<Integer, ArrayList<Integer>> e_children = new HashMap<Integer, ArrayList<Integer>>();

        // Aggregate values up
        for (int l = 0; l <= t.height; l++) {
            for (int i = 0; i < t.ind.length - 1; i++) {
                if (t.levels[i] == l) {
                    vals.put(i, t.wn[i]);
                    e_children.put(i, new ArrayList<Integer>());
                    for (int j = t.ind[i]; j < t.ind[i + 1]; j++) {
                        if (t.levels[t.endv[j]] == l - 1 &&
                                vals.get(t.endv[j]) - t.we[j] >= 0) {
                            vals.put(i, vals.get(i) + vals.get(t.endv[j]) - t.we[j]);
                            e_children.get(i).add(t.endv[j]);
                        }
                    }
                }
            }
        }
        int newcomp = 0;
        HashMap<Integer, Integer> comp = new HashMap<Integer, Integer>();

        // Update values down
        for (int l = t.height; l > 0; l--) {

            for (int i = 0; i < t.ind.length - 1; i++) {
                if (t.levels[i] == l) {
                    if (l == t.height)
                        comp.put(i, newcomp++);
                    for (int j = t.ind[i]; j < t.ind[i + 1]; j++) {
                        if (t.levels[t.endv[j]] == l - 1) {
                            double up = vals.get(i);
                            if (e_children.get(i).contains(t.endv[j]))
                                up = up - vals.get(t.endv[j]) + t.we[j];
                            if (up - t.we[j] > 0) {
                                vals.put(t.endv[j], vals.get(t.endv[j]) + up - t.we[j]);
                                comp.put(t.endv[j], comp.get(i));
                            } else
                                comp.put(t.endv[j], newcomp++);
                        }
                    }
                }
            }
        }

        Graph gnew = graph.clone();
        for (Edge e : graph.keySet()) {
            double score = Math.max(vals.get(t.n2i.get("" + e.source())), vals.get(t.n2i.get("" + e.dest())));
            if (graph.get(e) < 0 && !e_children.get(t.n2i.get("" + e.source())).contains(t.n2i.get("" + e.dest())) && !e_children.get(t.n2i.get("" + e.dest())).contains(t.n2i.get("" + e.source()))) {
                score += graph.get(e);
            }
            if (comp.get(t.n2i.get("" + e.source())) != comp.get(t.n2i.get("" + e.dest()))) {
                score += Math.min(vals.get(t.n2i.get("" + e.source())), vals.get(t.n2i.get("" + e.dest())));
            }
            gnew.put(e, (float) score);
        }
        return gnew;
    }
/*
    // TODO(petko): add the Steiner nodes along an expanded path
	// to consider in the next steps of the MST on the G_c
	// It basically selects edges and then does a subgraph of t
	public static TreeFast computeOptTreeNW(TreeFast t) {
		HashMap<Integer,Double> vals = new HashMap<Integer, Double>();
		HashMap<Integer,ArrayList<Integer>> e_children = new HashMap<Integer, ArrayList<Integer>>();

		// Aggregate values up
		for (int l=0; l<=t.height; l++) {
			for (int i = 0; i < t.ind.length-1; i++){
				if (t.levels[i]==l) {
					vals.put(i, t.wn[i]);
					e_children.put(i, new ArrayList<Integer>());
					for (int j=t.ind[i]; j < t.ind[i+1]; j++) {
						if(t.levels[t.endv[j]] == l-1 &&
						   vals.get(t.endv[j]) - t.we[j] >= 0) {
						   vals.put(i, vals.get(i) + vals.get(t.endv[j]) - t.we[j]);
						   e_children.get(i).add(j);
						}
					}
				}
			}
		}

		// Find Max
		int  ind = -1;
		double max = -1;
		for (int i:vals.keySet()) {
			if (vals.get(i) > max) {
				max = vals.get(i);
				ind = i;
			}
		}


		PriorityQueue<Integer> q = new PriorityQueue<Integer>();
		ArrayList<Integer> edges = new ArrayList<Integer>();
		q.add(ind);
		// collect the best subtree of ind
		Integer v = null;
		while (q.size() > 0) {
			v = q.poll();
			for(int edge:e_children.get(v)){
				edges.add(edge);
				edges.add(t.getEdgeIndex(t.endv[edge], v));
				assert(!q.contains(t.endv[edge]));
				q.add(t.endv[edge]);
			}
		}
		if(edges.size()==0) {
			return new TreeFast(new GraphFast(t.names[ind],t.wn[ind]));
		}

		return new TreeFast(t.subgraph(edges));
	}
*/

    public static void debugRefine(Subgraph patGraph, Subsequence patSeq) {
        System.out.print(patGraph.size() + " " + patSeq.size() + "\n");
        for (int rr = 0; rr < Math.min(10, patGraph.size()); rr++)
            System.out.print(patGraph.toArray()[rr] + " ");
        System.out.println();
        for (int rr = 0; rr < Math.min(10, patSeq.size()); rr++)
            System.out.print(patSeq.get(rr) + " ");
        System.out.println();

    }

    public static void main(String[] args) throws IOException {
        testMedenTopDown();
        DynaGraph dg = new DynaGraph();
        dg.read(args[0]);
        Pattern out = null;
        try {
            MaximumScoreDynamicSubgraph(dg, out, false);
        } catch (Exception ex) {
            System.err.print(ex.toString());
        }
    }
}
