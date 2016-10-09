import java.io.*;
import java.util.Scanner;
import java.util.Vector;


public class Subsequence extends Vector<Integer> {

    public void read(String fFileName) throws IOException {
        Scanner scanner = new Scanner(new FileInputStream(fFileName), Config.fEncoding);
        try {
            while (scanner.hasNextLine()) {
                String token = scanner.nextLine();
                Integer t;
                try {
                    t = new Integer(token);
                } catch (Exception e) {
                    continue;
                }
                if (t == null)
                    continue;
                this.add(t);
            }
        } finally {
            scanner.close();
        }
    }

    public String getPrintShort() {
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (int t : this) {
            min = Math.min(t, min);
            max = Math.max(t, max);
        }
        return "[" + min + "," + max + "]";
    }

    public void print() {
        System.out.print("Sequence: ");
        for (int t : this) {
            System.out.print("" + t + ",");
        }
        System.out.println("");
    }

    public void write(String filename) throws IOException {
        Writer writer = null;
        try {
            File file = new File(filename);
            writer = new BufferedWriter(new FileWriter(file));
            for (int i = 0; i < size(); i++) {
                writer.write("" + get(i) + "\n");
            }
        } finally {
            writer.close();
        }
    }

    public void write(Writer writer) throws IOException {
        try {
            for (int i = 0; i < size(); i++) {
                writer.write("" + get(i) + "\n");
            }
        } finally {
            writer.close();
        }
    }

    @Override
    public Object clone() {
        Subsequence newObj = new Subsequence();
        newObj.addAll(this);
        return newObj;

    }

    public boolean equals(Object ot) {
        Subsequence other = (Subsequence) ot;
        if (size() != other.size())
            return false;
        for (int i = 0; i < size(); i++) {
            if (get(i) != other.get(i))
                return false;
        }
        return true;
    }

    public int hashCode() {
        int hc = 0;
        int next = 1337;
        for (int i : this) {
            hc += i * next;
            next = next * 1337 % 12345;
        }
        return hc;
    }

    public void save(DataOutputStream out) throws IOException {
        out.writeInt(size());
        for (int i : this) {
            out.writeInt(i);
        }
    }

    public void load(DataInputStream in) throws IOException {
        int n = in.readInt();
        for (int i = 0; i < n; i++) {
            add(in.readInt());
        }
    }
}
