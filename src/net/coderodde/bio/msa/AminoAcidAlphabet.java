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
        alphabet.add('A');
        alphabet.add('C');
        alphabet.add('D');
        alphabet.add('E');
        alphabet.add('F');
        
        alphabet.add('G');
        alphabet.add('H');
        alphabet.add('I');
        alphabet.add('K');
        alphabet.add('L');
        
        alphabet.add('M');
        alphabet.add('N');
        alphabet.add('P');
        alphabet.add('Q');
        alphabet.add('R');
        
        alphabet.add('S');
        alphabet.add('T');
        alphabet.add('V');
        alphabet.add('W');
        alphabet.add('Y');
    }
    
    public Set<Character> getAlphabet() {
        return Collections.<Character>unmodifiableSet(alphabet);
    }
}
