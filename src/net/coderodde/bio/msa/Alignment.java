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
        sb.append("Cost: ").append(cost).append("\n");
        String separator = "";
        
        for (String row : alignment) {
            sb.append(separator).append(row);
            separator = "\n";
        }
        
        return sb.toString();
    }
}
