package net.coderodde.bio.msa;

import java.util.Objects;
import java.util.Set;

public final class MultipleSequenceAlignmentInstance {

    /**
     * The penalty for pairs in which there is one valid character and one gap.
     */
    private final int gapPenalty;
    
    /**
     * The character pairs cost matrix.
     */
    private final CostMatrix<Integer> costMatrix;
    
    /**
     * The array of sequences to be aligned.
     */
    private final String[] sequenceArray;
    
    public MultipleSequenceAlignmentInstance(CostMatrix<Integer> costMatrix,
                                             int gapPenalty,
                                             String... sequenceArray) {
        this.costMatrix = Objects.requireNonNull(costMatrix,
                                                 "Cost matrix is null");
        this.gapPenalty = gapPenalty;
        this.sequenceArray = sequenceArray.clone();
        
        for (int i = 0; i != sequenceArray.length; ++i) {
            this.sequenceArray[i] = 
                    checkIsValidGenomicSequence(sequenceArray[i]);
        }
    }

    /**
     * Create the source node for the sequence alignment task, i.e., a lattice 
     * node with all coordinates set to zero.
     * 
     * @return the source node.
     */
    LatticeNode getSourceNode() {
        int[] sourceCoordinates = new int[sequenceArray.length];
        return new LatticeNode(this, sourceCoordinates);
    }
    
    /**
     * Create the target node for the sequence alignment task, i.e., a lattice
     * node with all coordinates set to the respective sequence lengths.
     * 
     * @return the target node.
     */
    LatticeNode getTargetNode() {
        int[] targetCoordinates = new int[sequenceArray.length];
        
        for (int i = 0; i != sequenceArray.length; ++i) {
            targetCoordinates[i] = sequenceArray[i].length();
        }
        
        return new LatticeNode(this, targetCoordinates);
    }
    
    String[] getSequenceArray() {
        return sequenceArray;
    }
    
    private String checkIsValidGenomicSequence(String string) {
        Set<Character> characterSet = 
                AminoAcidAlphabet.getAminoAcidAlphabet()
                                 .getCharacterSet();
        
        for (char c : string.toCharArray()) {
            if (!characterSet.contains(c)) {
                throw new IllegalArgumentException("Unknown amino acid: " + c);
            }
        }
        
        return string;
    }
}
