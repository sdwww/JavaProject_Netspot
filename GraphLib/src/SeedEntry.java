import java.util.*;

public class SeedEntry implements Map.Entry<Seed, Float>, Comparable<SeedEntry> {

    private final Seed key;
    private Float value;

    public SeedEntry(Seed k, Float v) {
        key = k;
        value = v;
    }

    public Seed getKey() {
        return key;
    }

    public Float getValue() {
        return value;
    }

    public Float setValue(Float v) {
        Float old = value;
        value = v;
        return old;
    }

    public int compareTo(SeedEntry that) {
        return getValue().compareTo(that.getValue()) * -1;
    }

    public boolean equals(Object that) {
        if (this == that)
            return true;

        if (!(that instanceof SeedEntry))
            return false;

        SeedEntry thatS = (SeedEntry) that;

        return getKey().equals(thatS.getKey());
    }

    public int hashCode() {
        return getKey().hashCode();
    }
}
