package net.coderodde.bio.msa;

import java.awt.Point;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;

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
    
    // Basically, this method runs Dijkstra backwards from the target node until
    // the entire 2D-grid is explored.
    private void loadPartialHeuristicFunction(
            HeuristicFunction heuristicFunction,
            int dimension1,
            int dimension2,
            MultipleSequenceAlignmentInstance instance) {
        LatticeNode target = instance.getTargetNode();
        Queue<LatticeNodeHolder> open = new PriorityQueue<>();
        Set<LatticeNode> closed = new HashSet<>();
        open.add(new LatticeNodeHolder(target, 0));
        Point point = new Point();
        
        while (!open.isEmpty()) {
            LatticeNodeHolder currentNodeHolder = open.remove();
            LatticeNode currentNode = currentNodeHolder.getLatticeNode();
            
            if (closed.contains(currentNode)) {
                continue;
            }
            
            closed.add(currentNode);
            
            Integer currentNodeCost = currentNodeHolder.getCost();
            extractPoint(point, currentNode, dimension1, dimension2);
            heuristicFunction.put(dimension1,
                                  dimension2, 
                                  new Point(point),
                                  currentNodeCost);
            
            for (LatticeNode parent : currentNode.getParents()) {
                if (closed.contains(parent)) {
                    continue;
                }
                
                int tentativeCost =
                        currentNodeCost + instance.getWeight(parent, 
                                                             currentNode, 
                                                             dimension1, 
                                                             dimension2);
                extractPoint(point, parent, dimension1, dimension2);
                
                if (!heuristicFunction.containsPartial(dimension1,
                                                       dimension2, 
                                                       point)
                        || heuristicFunction.getPartial(dimension1,
                                                        dimension2,
                                                        point)
                        > tentativeCost) {
                    open.add(new LatticeNodeHolder(parent, tentativeCost));
                }
            }
        }
    }
    
    private static void extractPoint(Point point,
                                     LatticeNode latticeNode, 
                                     int i, 
                                     int j) {
        int[] coordinates = latticeNode.getCoordinates();
        point.x = coordinates[i];
        point.y = coordinates[j];
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
