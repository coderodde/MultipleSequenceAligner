package net.coderodde.bio.msa;

final class App {

    private static final String[] SEQUENCES = {
        "ACGHKGMNPFQW",
        "CFGPQWYRT",
        "EACLMNRPQWTR",
//        "TIMWAYH"
    };
    
    public static void main(String[] args) {
        MultipleSequenceAlignmentInstance instance = 
                new MultipleSequenceAlignmentInstance(
                        PAM250CostMatrix.getPAM250CostMatrix(),
                        4, 
                        SEQUENCES);
//        System.out.println(instance.align());
//        System.out.println(instance.alignBrute());

        long start = System.currentTimeMillis();
        Alignment alignment1 = instance.align();
        long end = System.currentTimeMillis();
        
        System.out.println(alignment1);
        System.out.println(end - start + " ms.");

        start = System.currentTimeMillis();
        Alignment alignment2 = instance.alignBrute();
        end = System.currentTimeMillis();
        
        System.out.println(alignment2);
        System.out.println(end - start + " ms.");
    }
}
