package org.cisrg.hyperstructures;

import java.util.ArrayList;


public abstract class GAPlugins {
	protected GAPlugins() { }
	
	public GAPlugins( double ub ) {
		this.upperBound = ub;
	}
	
	
	/**
	 * 
	 * @param chromosomes
	 * @param cindex
	 * @param probability	probability that a mutation is performed in a chromosome
	 * @param upperBound	the maximum index
	 * 
	 * Over-ride this at your own risk!  This standard mutation operator simply evaluates for each chromosome.  If 
	 * the random number generator falls below the specified probability for a given chromosome of the specified 
	 * set of chromosomes, then one of the integers is replaced with a random other integer between 0, and the
	 * upper bound.
	 */
	/*
	public ArrayList<Integer> mutate( ArrayList<Integer> chromosome, double probability, double upperBound  ) {
		//ArrayList<Integer> c2 = new ArrayList<Integer>(chromosomes.get(cindex));
		ArrayList<Integer> c2 = new ArrayList<Integer>( chromosome );
		
		for( int n = 0; n < c2.size(); n++ ) {

			if( Math.random() < probability ) {
				//c2.set(n, (int) Math.round( (Math.random() * upperBound)) );
				mutationType( c2, n );
			}
		}

		return c2;
	}
	*/
	
	
	/**
	 * Override this function with something that evaluates a chromosome, returning a value to represent how "good" it is.  
	 * The higher the value, the fitter.  For standardisation, it's recommended that this value is between 0 & 1.
	 * 
	 * @param chromosome
	 * @return
	 */
	public abstract double fitness(ArrayList<Integer> chromosome);
	
	public int mutationType( int n ) {
		return (int) Math.round( (Math.random() * upperBound ) );
	}
	
	
	public double upperBound;
}

