package net.coderodde.bio.msa;

final class App {

    private static final String[] SEQUENCES = {
        "ACGH",
        "CFG",
        "EAC",
    };
    
    public static void main(String[] args) {
        MultipleSequenceAlignmentInstance instance = 
                new MultipleSequenceAlignmentInstance(
                        PAM250CostMatrix.getPAM250CostMatrix(),
                        4, 
                        SEQUENCES);
        System.out.println(instance.align());
        System.out.println(instance.alignBrute());
    }
}
