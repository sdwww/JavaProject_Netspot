import java.util.Arrays;

public class TreeFast extends GraphFast {
    public int[] levels;
    public int height;
    public int root;

    public TreeFast(GraphFast _g) {
        super(_g);
        assert (isUndirected());
        if (ind != null) {
            levels = new int[getn()];
            computeLevels();
        }
    }

    private void computeLevels() {
        Arrays.fill(levels, 999);
        //  choose a root
        int root = super.getMaxOutDegree();
        levels[root] = 0;

        int covered = 1;
        int level = 0;
        while (covered < getn()) {
            for (int i = 0; i < getn(); i++) {
                if (levels[i] == level) {
                    for (int j = ind[i]; j < ind[i + 1]; j++) {
                        if (levels[endv[j]] == 999) {
                            levels[endv[j]] = level - 1;
                            covered++;
                        } else {
                            assert (level + 1 == levels[endv[j]]);
                        }
                    }
                }
            }
            level--;
        }
        height = -1 * level;
        for (int i = 0; i < getn(); i++) {
            levels[i] += height;
        }
    }

    public String toStringHuman() {

        String res = "Nodes:\n";
        if (ind == null) return res + names[0] + "(w:" + wn[0] + ")";
        for (int i = 0; i < names.length; i++) {
            res += names[i] + "(l:" + levels[i] + "; w:" + wn[i] + "),";
        }
        res += "\nEdges:\n";
        for (int i = 0; i < ind.length - 1; i++) {
            for (int j = ind[i]; j < ind[i + 1]; j++) {
                res += names[i] + "\t" + names[endv[j]] + "\t" + we[j] + "\n";
            }
        }
        return res;
    }
}
