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
    
    // A small speed optimization:
    private final char[] column;
    
    public MultipleSequenceAlignmentInstance(CostMatrix<Integer> costMatrix,
                                             int gapPenalty,
                                             String... sequenceArray) {
        this.costMatrix = Objects.requireNonNull(costMatrix,
                                                 "Cost matrix is null");
        this.gapPenalty = gapPenalty;
        this.sequenceArray = sequenceArray.clone();
        this.column = new char[sequenceArray.length];
        
        for (int i = 0; i != sequenceArray.length; ++i) {
            this.sequenceArray[i] = 
                    checkIsValidGenomicSequence(sequenceArray[i]);
        }
    }
    
    public Alignment align() {
        return null;
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
    
    Integer getWeight(LatticeNode tail, LatticeNode head) {
        // Extract the column represented by taking a single hop from 'tail' to
        // 'head':
        int[] tailCoordinates = tail.getCoordinates();
        int[] headCoordinates = head.getCoordinates();
        
        for (int i = 0; i < sequenceArray.length; ++i) {
            if (tailCoordinates[i] + 1 == headCoordinates[i]) {
                column[i] = sequenceArray[i].charAt(tailCoordinates[i]);
            } else {
                column[i] = AminoAcidAlphabet.GAP_CHARACTER;
            }
        }
        
        // Compute the hop cost as the sum of pairwise hops in any plane:
        int cost = 0;
        
        for (int i = 0; i < column.length; ++i) {
            char character1 = column[i];
            
            for (int j = i + 1; j < column.length; ++j) {
                char character2 = column[j];
                
                if (character1 != AminoAcidAlphabet.GAP_CHARACTER) {
                    if (character2 != AminoAcidAlphabet.GAP_CHARACTER) {
                        cost += costMatrix.getCost(character1, character2);
                    } else {
                        cost += gapPenalty;
                    }
                } else {
                    // character1 IS the gap character:
                    if (character2 != AminoAcidAlphabet.GAP_CHARACTER) {
                        cost += gapPenalty;
                    } else {
                        // Do nothing since we have a pair (gap, gap).
                    }
                }
            }
        }
        
        return cost;
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
