package net.coderodde.bio.msa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.PriorityQueue;
import java.util.Queue;
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
        HeuristicFunction hf = new HeuristicFunctionComputer()
                .computeHeuristicFunction(this);
        Queue<LatticeNodeHolder> open = new PriorityQueue<>();
        Map<LatticeNode, Integer> distance = new HashMap<>();
        Map<LatticeNode, LatticeNode> parents = new HashMap<>();
        
        LatticeNode sourceNode = getSourceNode();
        LatticeNode targetNode = getTargetNode();
        
        open.add(new LatticeNodeHolder(sourceNode, 0));
        distance.put(sourceNode, 0);
        parents.put(sourceNode, null);
        
        Set<LatticeNode> closed = new HashSet<>();
        
        while (true) {
            LatticeNode currentNode = open.remove().getNode();
            
            if (currentNode.equals(targetNode)) {
                return tracebackPath(parents, distance.get(currentNode));
            }
            
            if (closed.contains(currentNode)) {
                continue;
            }
            
            closed.add(currentNode);
            
            for (LatticeNode childNode : currentNode.getChildren()) {
                if (closed.contains(childNode)) {
                    continue;
                }
                
                int tentativeCost = distance.get(currentNode) +
                                    getWeight(currentNode, childNode);
                Integer currentCost = distance.get(childNode);
                
                if (currentCost == null || currentCost > tentativeCost) {
                    open.add(new LatticeNodeHolder(childNode, 
                                                   tentativeCost + 
                                                           hf.get(childNode)));
                    distance.put(childNode, tentativeCost);
                    parents.put(childNode, currentNode);
                }
            }
        }
    }
    
    public Alignment alignBrute() {
        Queue<LatticeNodeHolder> open = new PriorityQueue<>();
        Map<LatticeNode, Integer> distance = new HashMap<>();
        Map<LatticeNode, LatticeNode> parents = new HashMap<>();
        
        LatticeNode sourceNode = getSourceNode();
        LatticeNode targetNode = getTargetNode();
        
        open.add(new LatticeNodeHolder(sourceNode, 0));
        distance.put(sourceNode, 0);
        parents.put(sourceNode, null);
        
        while (true) {
            LatticeNode currentNode = open.remove().getNode();
            
            if (currentNode.equals(targetNode)) {
                return tracebackPath(parents, distance.get(currentNode));
            }
            
            for (LatticeNode childNode : currentNode.getChildren()) {
                int tentativeCost = distance.get(currentNode) +
                                    getWeight(currentNode, childNode);
                Integer currentCost = distance.get(childNode);
                
                if (currentCost == null || currentCost > tentativeCost) {
                    open.add(new LatticeNodeHolder(childNode, tentativeCost));
                    distance.put(childNode, tentativeCost);
                    parents.put(childNode, currentNode);
                }
            }
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
    
    Integer getWeight(LatticeNode tail,
                      LatticeNode head, 
                      int dimension1,
                      int dimension2) {
        int[] tailCoordinates = tail.getCoordinates();
        int[] headCoordinates = head.getCoordinates();
        
        if (tailCoordinates[dimension1] == headCoordinates[dimension1]) {
            //System.out.println("1");
            return gapPenalty;
        } else if (tailCoordinates[dimension2] == headCoordinates[dimension2]) {
            //System.out.println("2");
            return gapPenalty;
        } else {
//            System.out.println("3");
            char character1 = sequenceArray[dimension1]
                    .charAt(tailCoordinates[dimension1]);
            
            char character2 = sequenceArray[dimension2]
                    .charAt(tailCoordinates[dimension2]);
            return costMatrix.getCost(character1, character2);
        }
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
    
    private Alignment tracebackPath(Map<LatticeNode, LatticeNode> parents,
                                    Integer cost) {
        List<LatticeNode> path = new ArrayList<>();
        LatticeNode node = getTargetNode();
        
        while (node != null) {
            path.add(node);
            node = parents.get(node);
        }
        
        Collections.<LatticeNode>reverse(path);
        
        String[] strings = new String[getSequenceArray().length];
        StringBuilder[] stringBuilders = new StringBuilder[strings.length];
        
        for (int i = 0; i < stringBuilders.length; ++i) {
            stringBuilders[i] = new StringBuilder();
        }
        
        for (int i = 1; i < path.size(); ++i) {
            LatticeNode tail = path.get(i - 1);
            LatticeNode head = path.get(i);
            
            int[] tailCoordinates = tail.getCoordinates();
            int[] headCoordinates = head.getCoordinates();
            
            for (int j = 0; j < tailCoordinates.length; ++j) {
                if (tailCoordinates[j] != headCoordinates[j]) {
                    stringBuilders[j].
                            append(sequenceArray[j].charAt(tailCoordinates[j]));
                } else {
                    stringBuilders[j].append(AminoAcidAlphabet.GAP_CHARACTER);
                }
            }
        }
        
        for (int i = 0; i < strings.length; ++i) {
            strings[i] = stringBuilders[i].toString();
        }
        
        return new Alignment(strings, cost);
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
    
    private static final class LatticeNodeHolder
            implements Comparable<LatticeNodeHolder> {

        private final Integer fScore;
        private final LatticeNode node;
        
        LatticeNodeHolder(LatticeNode node, Integer fScore) {
            this.fScore = fScore;
            this.node = node;
        }
        
        LatticeNode getNode() {
            return node;
        }
        
        @Override
        public int compareTo(LatticeNodeHolder o) {
            return Integer.compare(fScore, o.fScore);
        }
    }
}
