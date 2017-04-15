package net.coderodde.bio.msa;

import java.awt.Point;
import java.util.HashMap;
import java.util.Map;

final class HeuristicFunction {

    private final Map<Point, Integer> map = new HashMap<>();
    
    void put(Point point, Integer cost) {
        map.put(point, cost);
    }
        
    boolean contains(Point point) {
        return map.containsKey(point);
    }
    
    Integer get(Point point) {
        return map.get(point);
    }
}
