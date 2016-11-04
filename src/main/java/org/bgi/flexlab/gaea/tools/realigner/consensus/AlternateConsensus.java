package org.bgi.flexlab.gaea.tools.realigner.consensus;

import htsjdk.samtools.Cigar;

import java.util.ArrayList;
import java.util.Arrays;

import org.bgi.flexlab.gaea.util.Pair;

public class AlternateConsensus {
	private final byte[] str;
	private final ArrayList<Pair<Integer, Integer>> readIndexes;
	private final int positionOnReference;
	private int mismatchSum;
	private Cigar cigar;

    public AlternateConsensus(byte[] str, Cigar cigar, int positionOnReference) {
        this.str = str;
        this.cigar = cigar;
        this.positionOnReference = positionOnReference;
        mismatchSum = 0;
        readIndexes = new ArrayList<Pair<Integer, Integer>>();
    }

    @Override
    public boolean equals(Object o) {
    	if(o instanceof AlternateConsensus){
    		AlternateConsensus other = (AlternateConsensus)o;
    		return ( this == other || Arrays.equals(this.str,other.getSequence()) ) ;
    	}
    	return false;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(this.str);
    }
    
    public byte[] getSequence(){
    	return this.str;
    }
    
    public ArrayList<Pair<Integer, Integer>> getReadIndexes(){
    	return this.readIndexes;
    }
    
    public int getPositionOnReference(){
    	return this.positionOnReference;
    }
    
    public int getMismatch(){
    	return this.mismatchSum;
    }
    
    public Cigar getCigar(){
    	return this.cigar;
    }
}
