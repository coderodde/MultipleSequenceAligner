package net.coderodde.bio.msa;

import java.awt.Point;
import java.util.PriorityQueue;
import java.util.Queue;

final class HeuristicFunctionComputer {

    HeuristicFunction computeHeuristicFunction(
            MultipleSequenceAlignmentInstance instance) {
        String[] sequenceArray = instance.getSequenceArray();
        HeuristicFunction heuristicFunction = new HeuristicFunction(instance);
        
        for (int dimension1 = 0; 
                dimension1 < sequenceArray.length; 
                dimension1++) {
            for (int dimension2 = dimension1 + 1; 
                    dimension2 < sequenceArray.length; 
                    dimension2++) {
                loadPartialHeuristicFunction(heuristicFunction,
                                             dimension1,
                                             dimension2,
                                             instance);
            }
        }
        
        return heuristicFunction;
    }
    
    private void loadPartialHeuristicFunction(
            HeuristicFunction heuristicFunction,
            int dimension1,
            int dimension2,
            MultipleSequenceAlignmentInstance instance) {
        LatticeNode target = instance.getTargetNode();
        Queue<LatticeNodeHolder> open = new PriorityQueue<>();
        heuristicFunction.put(dimension1,
                              dimension2,
                              extractPoint(target, 
                                           dimension1, 
                                           dimension2),
                              0);
        
        while (!open.isEmpty()) {
            LatticeNodeHolder currentNodeHolder = open.remove();
            LatticeNode currentNode = currentNodeHolder.getLatticeNode();
            Integer currentNodeCost = currentNodeHolder.getCost();
            
            for (LatticeNode parent : currentNode.getParents()) {
                int tentativeCost = 
                        currentNodeCost + instance.getWeight(parent, 
                                                             currentNode);
                Point parentPoint = extractPoint(parent, 
                                                 dimension1, 
                                                 dimension2);
                
                if (!heuristicFunction.containsPartial(dimension1,
                                                       dimension2, 
                                                       parentPoint)
                        || heuristicFunction.getPartial(dimension1,
                                                        dimension2,
                                                        parentPoint)
                        > tentativeCost) {
                    open.add(new LatticeNodeHolder(parent, tentativeCost));
                }
            }
        }
    }
    
    private static Point extractPoint(LatticeNode latticeNode, int i, int j) {
        return new Point(latticeNode.getCoordinates()[i],
                         latticeNode.getCoordinates()[j]);
    }
    
    private static final class LatticeNodeHolder 
         implements Comparable<LatticeNodeHolder> {

        private final LatticeNode node;
        private final Integer cost;
        
        LatticeNodeHolder(LatticeNode node, Integer cost) {
            this.node = node;
            this.cost = cost;
        }
        
        LatticeNode getLatticeNode() {
            return node;
        }
        
        Integer getCost() {
            return cost;
        }
        
        @Override
        public int compareTo(LatticeNodeHolder o) {
            return Integer.compare(cost, o.cost);
        }
    }
}
