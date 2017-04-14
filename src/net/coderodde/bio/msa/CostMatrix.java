package net.coderodde.bio.msa;

public interface CostMatrix<C> {
    
    public C getCost(Character aminoAcidChar1, Character aminoAcidChar2);
}
