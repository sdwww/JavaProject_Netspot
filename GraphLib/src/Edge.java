import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;


public class Edge {
    private int source;
    private int dest;

    public int source() {
        return source;
    }

    public int dest() {
        return dest;
    }

    public Edge(int ls, int ld) {
        source = Math.min(ls, ld);
        dest = Math.max(ls, ld);
    }

    @Override
    public boolean equals(Object obj) {
        Edge other = (Edge) obj;
        return (source == other.source) && (dest == other.dest);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public int hashCode() {
        return 123456791 * source + dest;
    }

    @Override
    public String toString() {
        return "(" + source + ", " + dest + ")";

    }

    public void save(DataOutputStream out) throws IOException {
        out.writeInt(source);
        out.writeInt(dest);
    }

    public void load(DataInputStream in) throws IOException {
        source = in.readInt();
        dest = in.readInt();
    }
}
