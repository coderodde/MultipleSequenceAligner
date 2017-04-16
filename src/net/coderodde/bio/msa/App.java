package net.coderodde.bio.msa;

final class App {

    private static final String[] SEQUENCES = {
        "ACGHKGMNPFQEKKFKLMNRW",
        "CFGPQWYRTLMEKKFKNR",
        "EACLMNRPQWTR",
        "TIMWAYHTMGIEKKFK"
    };
    
    public static void main(String[] args) {
        MultipleSequenceAlignmentInstance instance = 
                new MultipleSequenceAlignmentInstance(
                        PAM250CostMatrix.getPAM250CostMatrix(),
                        8, 
                        SEQUENCES);
        
        long start = System.currentTimeMillis();
        Alignment alignment1 = instance.align();
        long end = System.currentTimeMillis();
        
        System.out.println(alignment1);
        System.out.println(end - start + " ms.");
        System.out.println();

        start = System.currentTimeMillis();
        Alignment alignment2 = instance.alignBrute();
        end = System.currentTimeMillis();
        
        System.out.println(alignment2);
        System.out.println(end - start + " ms.");
    }
}
