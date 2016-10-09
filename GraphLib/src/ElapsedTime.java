

public class ElapsedTime {
    public static ElapsedTime Start() {
        return new ElapsedTime();
    }

    public static float GetPartial(ElapsedTime handle) {
        return handle.GetElapsedTime();
    }

    public static float Stop(ElapsedTime handle) {
        float time = handle.GetElapsedTime();
        try {
            handle.finalize();
        } catch (Throwable e) {
            System.out.println("Error: object not finalized");
        }
        return time;
    }

    private long start;

    private ElapsedTime() {
        start = System.currentTimeMillis();
    }

    private float GetElapsedTime() {
        long elapsedTimeMillis = System.currentTimeMillis() - start;
        float elapsedTimeSec = elapsedTimeMillis / 1000F;
        return elapsedTimeSec;
    }
}
