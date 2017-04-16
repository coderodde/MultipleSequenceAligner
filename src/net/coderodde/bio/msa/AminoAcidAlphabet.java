package net.coderodde.bio.msa;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public final class AminoAcidAlphabet {

    public static final Character GAP_CHARACTER = '-';
    private static AminoAcidAlphabet instance;
    
    private final Set<Character> alphabet = new HashSet<>();
    
    public static AminoAcidAlphabet getAminoAcidAlphabet() {
        if (instance == null) {
            instance = new AminoAcidAlphabet();
        }
        
        return instance;
    }
    
    private AminoAcidAlphabet() {
        String aminoAcids = "ACDEF" +
                            "GHIKL" +
                            "MNPQR" +
                            "STVWY";
        
        for (char c : aminoAcids.toCharArray()) {
            alphabet.add(c);
        }
    }
    
    public Set<Character> getCharacterSet() {
        return Collections.<Character>unmodifiableSet(alphabet);
    }
}
