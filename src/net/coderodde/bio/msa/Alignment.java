package net.coderodde.bio.msa;

public class Alignment {

    private final String[] alignment;
    private final int cost;

    Alignment(String[] alignment, int cost) {
        this.alignment = alignment;
        this.cost = cost;
    }

    public String[] getAlignemnt() {
        return alignment;
    }

    public int getCost() {
        return cost;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        String separator = "";

        for (String row : alignment) {
            sb.append(separator).append(row);
            separator = "\n";
        }

        sb.append("\nCost: ").append(cost);
        return sb.toString();
    }
}
