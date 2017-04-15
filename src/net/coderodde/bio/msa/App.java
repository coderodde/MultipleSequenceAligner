package net.coderodde.bio.msa;

import net.coderodde.bio.msa.support.PAM250CostMatrix;

final class App {

    private static final String[] SEQUENCES = {
        "ACCD",
        "ADC",
        "CDA",
    };
    
    public static void main(String[] args) {
        MultipleSequenceAlignmentInstance instance = 
                new MultipleSequenceAlignmentInstance(
                        PAM250CostMatrix.getPAM250CostMatrix(),
                        8, 
                        SEQUENCES);
    }
}
