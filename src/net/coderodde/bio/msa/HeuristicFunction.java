package net.coderodde.bio.msa;

import java.awt.Point;
import java.util.HashMap;
import java.util.Map;

final class HeuristicFunction {

    private final Map<Point, Integer>[][] matrix;
    
    HeuristicFunction(MultipleSequenceAlignmentInstance instance) {
        int sequences = instance.getSequenceArray().length;
        this.matrix = new Map[sequences][];
        
        for (int i = 0; i < sequences; ++i) {
            for (int j = i + 1; j < sequences; ++j) {
                this.matrix[i][j] = new HashMap<>();
            }
        }
    }
    
    void put(int dimension1, int dimension2, Point point, Integer cost) {
        matrix[dimension1][dimension2].put(point, cost);
    }
        
    boolean containsPartial(int dimension1, int dimension2, Point point) {
        return matrix[dimension1][dimension2].containsKey(point);
    }
    
    Integer getPartial(int dimension1, int dimension2, Point point) {
        return matrix[dimension1][dimension2].get(point);
    }
    
    Integer get(LatticeNode node) {
        int cost = 0;
        int[] coordinates = node.getCoordinates();
        Point point = new Point();
        
        for (int dimension1 = 0; 
                dimension1 < coordinates.length; 
                dimension1++) {
            point.x = dimension1;
            
            for (int dimension2 = dimension1 + 1;
                    dimension2 < coordinates.length; 
                    dimension2++) {
                point.y = dimension2;
                cost += matrix[dimension1][dimension2].get(point);
            }
        }
        
        return cost;
    }
    
    private static Point extractPoint(int[] coordinates, int i, int j) {
        return new Point(coordinates[i], coordinates[j]);
    }
}
