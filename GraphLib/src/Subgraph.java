import java.io.*;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Scanner;


public class Subgraph extends HashSet<Edge> {

    public void read(String fFileName) throws IOException {
        Scanner scanner = new Scanner(new FileInputStream(fFileName), Config.fEncoding);
        LinkedList<Edge> l = new LinkedList<Edge>();
        try {
            while (scanner.hasNextLine()) {
                Edge e = processLine(scanner.nextLine());
                if (e == null)
                    continue;
                this.add(new Edge(e.source(), e.dest()));
            }
        } finally {
            scanner.close();
        }
    }

    static public Edge processLine(String aLine) throws IOException {
        //use a second Scanner to parse the content of each line
        Scanner scanner = new Scanner(aLine);
        scanner.useDelimiter(",");
        int s, d;
        if (scanner.hasNext()) {
            String token = scanner.next();
            s = new Integer(token);

            if (!scanner.hasNext())
                throw new IOException("Incomplete row");
            token = scanner.next();
            d = new Integer(token);
            if (s > d) {
                int dest = s;
                s = d;
                d = dest;
            } //swap

            return new Edge(s, d);
        }
        return null;
    }

    public void print() {
        System.out.print("Graph: ");
        for (Edge e : this) {
            System.out.print("<" + e.source() + "," + e.dest() + "> ");
        }
        System.out.println("");
    }

    public String getPrintShort() {
        HashSet<Integer> nodes = new HashSet<Integer>();
        for (Edge e : this) {
            nodes.add(e.source());
            nodes.add(e.dest());
        }
        String ret = "";
        for (int i : nodes)
            ret += i + "-";
        return ret;
    }

    public void printCytoscape(String fFileName) throws IOException {
        Writer out = new OutputStreamWriter(new FileOutputStream(fFileName), Config.fEncoding);
        try {
            for (Edge e : this) {
                out.write(String.valueOf(e.source()) + " (pp) " + String.valueOf(e.dest()) + "\t1\n");
            }
        } finally {
            out.close();
        }
    }

    public void write(String filename) throws IOException {
        Writer writer = null;
        try {
            File file = new File(filename);
            writer = new BufferedWriter(new FileWriter(file));
            for (Edge e : this) {
                writer.write("" + e.source() + "," + e.dest() + "\n");
            }
        } finally {
            writer.close();
        }
    }

    public void write(Writer writer) throws IOException {
        try {
            for (Edge e : this) {
                writer.write("" + e.source() + "," + e.dest() + "\n");
            }
        } finally {
            writer.close();
        }
    }

    @Override
    public Object clone() {
        Subgraph newObj = new Subgraph();
        newObj.addAll(this);
        return newObj;
    }

    public boolean equals(Object ot) {
        Subgraph other = (Subgraph) ot;
        if (size() != other.size())
            return false;
        for (Edge e : this) {
            if (!other.contains(e))
                return false;
        }
        return true;
    }

    public int hashCode() {
        int hc = 0;
        for (Edge e : this)
            hc += e.hashCode();
        return hc;
    }

    public void save(DataOutputStream out) throws IOException {
        out.writeInt(size());
        for (Edge e : this) {
            e.save(out);
        }
    }

    public void load(DataInputStream in) throws IOException {
        int n = in.readInt();
        for (int i = 0; i < n; i++) {
            Edge e = new Edge(0, 0);
            e.load(in);
            add(e);
        }
    }
}
