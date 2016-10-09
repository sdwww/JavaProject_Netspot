import java.util.Vector;


public class Sequence extends Vector<Float> {
    public Subsequence MaximumScoreSubsequence() {
        int i = 0;
        float s = 0;
        float sMax = Float.NEGATIVE_INFINITY;
        int iMax = 0, jMax = -1;
        for (int j = 0; j < this.size(); j++) {
            float w = this.get(j);
            s += w;
            if (w > s) {
                i = j;
                s = w;
            }
            if (s > sMax) {
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

    public Subsequence MaximumScoreSubsequence(int pivot) {
        // the pivot is constrained to belong to the Maximum Score Subsequence
        Subsequence subseq = new Subsequence();

        int i = 0;
        float s = 0;
        float sMax = Float.NEGATIVE_INFINITY;
        int iMax = 0, jMax = -1;
        for (int j = 0; (j < this.size()) && (i <= pivot); j++) {
            float w = this.get(j);
            s += w;
            if ((w > s) && (j <= pivot)) {
                i = j;
                s = w;
            }
            if ((s > sMax) && (i <= pivot) && (pivot <= j)) {
                sMax = s;
                iMax = i;
                jMax = j;
            }
        }

        for (int t = iMax; t <= jMax; t++) {
            subseq.add(t);
        }
        return subseq;
    }

    float GetTotalWeight() {
        float sc = 0F;
        for (int i = 0; i < this.size(); i++) {
            sc += this.get(i);

        }
        return sc;
    }
}
