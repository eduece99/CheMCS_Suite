package org.cisrg.mapping;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.cisrg.BitSetExtended;
import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.GenerateCompatibilityGraphEdges.EdgeType;
import org.openscience.cdk.interfaces.IAtomContainer;

public class GenerateCompatibilityGraphEdgesEdgeTypes extends
		GenerateCompatibilityGraphEdges {

	
	
	protected List<Collection<Integer>> cEdges;  // c-edge - the modular product nodes are both joined to each other in each molecule
	protected List<Collection<Integer>> dEdges;  // d-edge - the modular product nodes are both NOT joined to each other in each molecule
	protected List<Collection<Integer>> rEdges;  // red edges - complement graph of modular product (so basically, an absence of edge)
	
	
	
	
	
	/**
     * Generates a compatibility graph between two molecules
     * 
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @throws java.io.IOException
     */
    public GenerateCompatibilityGraphEdgesEdgeTypes(IAtomContainer source,
            IAtomContainer target,
            boolean shouldMatchBonds, 
            boolean useRaymondHeuristics,
            boolean useRingHeuristics,
            int topoDistance,
            String atomTypeProperty,
            boolean buildEdges
         ) throws IOException {
    	
    	super( source, target, shouldMatchBonds, useRaymondHeuristics, useRingHeuristics, topoDistance, atomTypeProperty, buildEdges );
    	
    	// red edge handling
    	for (int a = 0; a < edges.size(); a++ ) { 
    		edges.get(a);
    		
    		BitSetExtended<Integer> inverted = new BitSetExtended<>( edges.get(a) );
    		inverted.getBitSet().flip(0, edges.size() );
    		//inverted.remove( a );
    		rEdges.add( inverted );
    	}
    	
    }
	
    
    
    @Override
    protected void createEdgeDataStructures(int nodesCount) {
    	
    	
    	edges = new ArrayList<Collection<Integer>>( nodesCount );
    	cEdges = new ArrayList<Collection<Integer>>( nodesCount );
    	dEdges = new ArrayList<Collection<Integer>>( nodesCount );
    	rEdges = new ArrayList<Collection<Integer>>( nodesCount );
    	
    	for (int a = 0; a < nodesCount; a++ ) {
        	edges.add( ConvenienceTools.createCollection( nodesCount ) );
        	cEdges.add( ConvenienceTools.createCollection( nodesCount ) );
        	dEdges.add( ConvenienceTools.createCollection( nodesCount ) );
        	//rEdges.add( ConvenienceTools.createCollection( nodesCount ) );
        }
    }
    
    
    @Override
	protected void addEdge( int a, int b, EdgeType type ) {
    	edges.get(a).add(b);  
        edges.get(b).add(a);
        //System.err.println("test");
        
       if( type == EdgeType.CEDGE ) {
        	cEdges.get(a).add(b);  
            cEdges.get(b).add(a);
        } else if( type == EdgeType.DEDGE ) {
        	dEdges.get(a).add(b);  
            dEdges.get(b).add(a);
        } else if( type == EdgeType.REDEDGE ) {
        	rEdges.get(a).add(b);  
            rEdges.get(b).add(a);
        }
        
        numberOfEdges += 2;
    }
    
    
    /**
     * specify a row in the modular products adjacency list/matrix.  Returns the intersection of the supplied set, with said row.
     * 
     * @param index
     * @param set
     * @return
     */
    public Collection<Integer> intersectNodeNeighbours( int index, Collection<Integer> set, EdgeType eType ) {
    	
    	if( set.isEmpty() )
    		return new ArrayList<Integer>(0);
    	
    	Collection<Integer> neighbours = null;
    	
    	switch( eType ) {
    		case CEDGE:
    			neighbours = cEdges.get(index);
    		case DEDGE:
    			neighbours = dEdges.get(index);
    		default:
    			neighbours = edges.get(index);
    	}
    	
    

    	Collection<Integer> oSet = null;
    	oSet = ConvenienceTools.createCollection(set);
   
    	
    	//System.out.println( "class of collection is " + oSet.getClass().toString() );
    	oSet.retainAll(neighbours);
    	
    	return oSet;
    }
    
    
    
    public List<Collection<Integer>> getCEdges() {
    	return cEdges;
    }
    
    public List<Collection<Integer>> getREdges() {
    	return rEdges;
    }
    
}
