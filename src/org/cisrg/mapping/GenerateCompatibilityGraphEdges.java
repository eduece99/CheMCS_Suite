/*
 *
 * Copyright (C) 2006-2010  Syed Asad Rahman <asad@ebi.ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received iIndex copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.cisrg.mapping;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.hyperstructures.DefaultHyperstructureAtomMatcher;
import org.openscience.cdk.CDKConstants;
//import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.ConnectedComponents;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.cdk.smsd.algorithm.matchers.DefaultBondMatcher;
import org.openscience.cdk.tools.manipulator.BondManipulator;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also markes edges in the compatibility graph as c-edges or d-edges.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
//@TestClass("org.openscience.cdk.smsd.SMSDBondSensitiveTest")
public class GenerateCompatibilityGraphEdges {

	protected List<int[]> graphNodes = null;  // each node is a 2-long array containing references to the source and target molecule bonds respectively
	protected List<Collection<Integer>> edges = null;  // contains indices of the graphNodes List.  Collection type for BitString extensions, Adj List etc...
	
	protected IAtomContainer source = null;
	protected IAtomContainer target = null;
	protected boolean shouldMatchBonds = false;
	protected String atomTypeProperty = null;
	protected boolean verbose = false;
    
    
    // heuristic variables
	protected boolean raymondHeuristics = true;
	protected boolean ringHeuristics = false;
	protected boolean useTopologicalDistance = false;
    protected int topologicalDistanceLimit = 500;
    protected int flexTopoDistLim = 1;
    public int[][] pathDistancesHsMol;
    public int[][] pathDistancesQMol;
    List<List<Integer>> sourceAdjList;
    List<List<Integer>> targetAdjList;
	List<List<Integer>> sourceRings; 
	List<List<Integer>> targetRings;
	List<List<Integer>> sourceBondRings; 
	List<List<Integer>> targetBondRings;
	List<Integer> sourceCCInducedSubgraphBondIDs; 
	List<Integer> targetCCInducedSubgraphBondIDs;
	List<Integer> sourceRigidInducedSubgraphBondIDs; 
	List<Integer> targetRigidInducedSubgraphBondIDs;
	List<Integer> sourceFlexibleInducedSubgraphBondIDs; 
	List<Integer> targetFlexibleInducedSubgraphBondIDs;
	int[] sourceFlexibleInducedSubgraphPathDistances;
	int[] targetFlexibleInducedSubgraphPathDistances;
	List<List<Integer>> sourceRingBonds;
	List<List<Integer>> targetRingBonds;
	
	protected SmilesGenerator sg;
	
	
    public int numberOfEdges = 0;
    
    
    public enum inducedSubgraphType {
    	CC,
    	RIGID,
    	NONRING
    }
    
    public enum EdgeType {
    	LOOP,
    	REDEDGE,
    	CEDGE,
    	DEDGE
    }
    
    
    // ascending order of first element
 	protected Comparator<int[]> pairComparator = new Comparator<int[]>() {
 					 public int compare(int[] o1, int[] o2) {

 				            int dComp;
 				            // sort in descending order of degree first
 				            if( o1[0] == o2[0] ) {
 				            	dComp = 0;
 				            } else if( o1[0] < o2[0] ) {
 				            	dComp = -1;
 				            	//Collections.swap( sortedVertices2, o1, o2 );
 				            } else {
 				            	dComp = 1;
 				            	//Collections.swap( sortedVertices2, o2, o1 );
 				            }
 				            
 				            return dComp;
 					 } 
 				};
 				
 				

    /**
    * Default constructor added 
    */
    public GenerateCompatibilityGraphEdges(){
        
    }

    /**
     * Generates a compatibility graph between two molecules
     * 
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @throws java.io.IOException
     */
    public GenerateCompatibilityGraphEdges(IAtomContainer source,
            IAtomContainer target,
            boolean shouldMatchBonds, 
            boolean useRaymondHeuristics,
            boolean useRingHeuristics,
            int topoDistance,
            String atomTypeProperty,
            boolean buildEdges
         ) throws IOException {
    	
        setMatchBond(shouldMatchBonds);
        this.source = source;
        this.target = target;
        this.atomTypeProperty = atomTypeProperty;
        this.raymondHeuristics = useRaymondHeuristics;
        this.ringHeuristics = useRingHeuristics;
        this.topologicalDistanceLimit = topoDistance;
        
        sg = new SmilesGenerator().aromatic();
        
        graphNodes = new ArrayList<int[]>( source.getBondCount() * target.getBondCount() );
        
        
        
        //adjMatrix = new IBond[source.getAtomCount()][target.getAtomCount()];  // null = no bond/edge
        
        
        
        useTopologicalDistance = ( topoDistance >= 0 && topoDistance < Math.max( source.getBondCount(), target.getBondCount() ) );
        
        if( raymondHeuristics || useTopologicalDistance ) {
        	
	        pathDistancesHsMol = PathTools.computeFloydAPSP( ConvenienceTools.bondAdjacencyMatrix(source) );
	        //pathDistancesHsMol = ConvenienceTools.bfsShortestPathLengths( ConvenienceTools.createBondAdjacencyList(source) );
			pathDistancesQMol = PathTools.computeFloydAPSP( ConvenienceTools.bondAdjacencyMatrix(target) );
	        //pathDistancesQMol = ConvenienceTools.bfsShortestPathLengths( ConvenienceTools.createBondAdjacencyList(target) );
        }
        
		if( raymondHeuristics ) {
			/*
			ConvenienceTools.countRings(source);
			ConvenienceTools.countRings(target);
			
			
			try {
				ConvenienceTools.calculateAromaticity( source );
				ConvenienceTools.calculateAromaticity( target );
			} catch (CDKException e) {
				e.printStackTrace();
			}
			*/
			
			
			if( verbose ) {
				try {
					System.out.println( "source SMILES - " + sg.create( source ) );
					System.out.println( "target SMILES - " + sg.create( target ) );
				} catch (CDKException e1) {
					e1.printStackTrace();
				}
			}
			
			sourceAdjList = ConvenienceTools.createBondAdjacencyList(source);
			targetAdjList = ConvenienceTools.createBondAdjacencyList(target);
			
			
			sourceRings = (List<List<Integer>>) source.getProperty( CDKConstants.RELEVANT_RINGS ); 
			targetRings = (List<List<Integer>>) target.getProperty( CDKConstants.RELEVANT_RINGS );
			
			
			if( buildEdges ) {
			
				sourceBondRings = new ArrayList<List<Integer>>( source.getBondCount() );
				targetBondRings = new ArrayList<List<Integer>>( target.getBondCount() );
				
				
				
				for( int n=0; n < source.getBondCount(); n++ ) {
					sourceBondRings.add( (List<Integer>) source.getBond( n ).getProperty( CDKConstants.RING_CONNECTIONS ) );
				}
				
				for( int n=0; n < target.getBondCount(); n++ ) {
					targetBondRings.add( (List<Integer>) target.getBond( n ).getProperty( CDKConstants.RING_CONNECTIONS ) );
				}
				
				sourceCCInducedSubgraphBondIDs = findInducedSubgraphBonds( source, inducedSubgraphType.CC );
				targetCCInducedSubgraphBondIDs = findInducedSubgraphBonds( target, inducedSubgraphType.CC );
				
				sourceRigidInducedSubgraphBondIDs = findInducedSubgraphBonds( source, inducedSubgraphType.RIGID );
				targetRigidInducedSubgraphBondIDs = findInducedSubgraphBonds( target, inducedSubgraphType.RIGID );
				
				sourceFlexibleInducedSubgraphBondIDs = findInducedSubgraphBonds( source, inducedSubgraphType.NONRING );
				targetFlexibleInducedSubgraphBondIDs = findInducedSubgraphBonds( target, inducedSubgraphType.NONRING );
				
				sourceFlexibleInducedSubgraphPathDistances = new int[ sourceFlexibleInducedSubgraphBondIDs.size() / 2 + 1 ];
				targetFlexibleInducedSubgraphPathDistances = new int[ targetFlexibleInducedSubgraphBondIDs.size() / 2 + 1 ];
				
				
				
				// XXX  find max path length in chain-based subgraphs, of either molecule
		    	for( int b1 = 0; b1 < sourceFlexibleInducedSubgraphBondIDs.size(); b1++ ) {
		    		int b1sgId = sourceFlexibleInducedSubgraphBondIDs.get(b1);
		    		
		    		if( b1sgId < 0 )
		    			continue;
		    		
		    		for( int b2 = 0; b2 < sourceFlexibleInducedSubgraphBondIDs.size(); b2++ ) {
		    			if( sourceFlexibleInducedSubgraphBondIDs.get(b1) == sourceFlexibleInducedSubgraphBondIDs.get(b2) ) {
		    				if( sourceFlexibleInducedSubgraphPathDistances[ b1sgId ] < pathDistancesHsMol[b1][b2] )
		    					sourceFlexibleInducedSubgraphPathDistances[ b1sgId ] = pathDistancesHsMol[b1][b2];
		    				
		    				if( flexTopoDistLim < pathDistancesHsMol[b1][b2] )
		    					flexTopoDistLim = pathDistancesHsMol[b1][b2];
		    			}
		    		}
		    	}
		    	
		    	for( int b1 = 0; b1 < targetFlexibleInducedSubgraphBondIDs.size(); b1++ ) {
		    		int b1sgId = targetFlexibleInducedSubgraphBondIDs.get(b1);
		    		
		    		if( b1sgId < 0 )
		    			continue;
		    		
		    		for( int b2 = 0; b2 < targetFlexibleInducedSubgraphBondIDs.size(); b2++ ) {
		    			if( targetFlexibleInducedSubgraphBondIDs.get(b1) == targetFlexibleInducedSubgraphBondIDs.get(b2) ) {
		    				if( targetFlexibleInducedSubgraphPathDistances[ b1sgId ] < pathDistancesQMol[b1][b2] )
		    					targetFlexibleInducedSubgraphPathDistances[ b1sgId ] = pathDistancesQMol[b1][b2];
		    				
		    				if( flexTopoDistLim < pathDistancesQMol[b1][b2] )
		    					flexTopoDistLim = pathDistancesQMol[b1][b2];
		    			}
		    		}
		    	}
			}
			
			//compatibilityGraphNodes( );
			compatibilityGraphNodesRaymondHeuristics();
			
			
			if( buildEdges ) {
				createEdgeDataStructures( graphNodes.size() );
				
				int maxSRingIndex = 0;
				int maxTRingIndex = 0;
				
				for( int b=0; b < sourceBondRings.size(); b++ ) {
					if( sourceBondRings.get(b) == null )
						continue;
					
					int currentSize = Collections.max(sourceBondRings.get(b));
					if( currentSize > maxSRingIndex ) 
						maxSRingIndex = currentSize;
				}
				
				for( int b=0; b < targetBondRings.size(); b++ ) {
					if( targetBondRings.get(b) == null )
						continue;
					
					int currentSize = Collections.max(targetBondRings.get(b));
					if( currentSize > maxTRingIndex ) 
						maxTRingIndex = currentSize;
				}
			
				sourceRingBonds = new ArrayList<List<Integer>>( maxSRingIndex+1 );
				targetRingBonds = new ArrayList<List<Integer>>( maxTRingIndex+1 );
				
				// initial autoboxing
				for( int b=0; b <= maxSRingIndex+1; b++ ) {
					sourceRingBonds.add( new ArrayList<Integer>() );
				}
				for( int b=0; b <= maxTRingIndex+1; b++ ) {
					targetRingBonds.add( new ArrayList<Integer>() );
				}
				
				
				
				for( int b=0; b < sourceBondRings.size(); b++ ) {
					List<Integer> bonds = sourceBondRings.get(b);
					
					if( bonds == null )
						continue;
					
					for( Integer bIndex : bonds ) {
						sourceRingBonds.get( bIndex ).add(b);
					}
				}
				
				//ConvenienceTools.printList(sourceRingBonds);
				
				for( int b=0; b < targetBondRings.size(); b++ ) {
					List<Integer> bonds = targetBondRings.get(b);
					
					if( bonds == null )
						continue;
					
					for( Integer bIndex : bonds ) {
						targetRingBonds.get( bIndex ).add(b);
					}
				}
				
				if( verbose ) {
					System.out.println( "sourceRingBonds " + sourceRingBonds );
					System.out.println( "targetRingBonds " + targetRingBonds );
				}
				compatibilityGraphEdgesRaymondHeuristics();
			}
			
			//System.out.println( sourceCCInducedSubgraphBondIDs );
		} else {
	        compatibilityGraphNodes( );
	        
	        if( buildEdges ) {
	        	createEdgeDataStructures( graphNodes.size() );
	        	compatibilityGraphEdges( );
	        }
		}
		//boolean finished = false;
		/*        
        for( int[] ea : getNodes() ) {
			IBond bondA = source.getBond( ea[0] );
			IBond bondB = target.getBond( ea[1] );
			System.out.println(  " [ " + ea[0] + "," + ea[1] + "]  " + 
					bondA.getOrder() + " " + bondA.getAtom(0).getSymbol() + " " +  bondA.getAtom(1).getSymbol() + " " +
					bondB.getOrder() + " " + bondB.getAtom(0).getSymbol() + " " +  bondB.getAtom(1).getSymbol() + " " );
		}*/
    }

    
    
    
    /*
     * Edge modular product:
     * 
     * Go through all edges of graphs
     * Create a node to represent 2 edges, if they match (same bond type and 2 atoms are of compatible labels) - node is the 2 edges
     * An edge is if the 2 edges are: non-identical (that is, both vertices are not the same); adjacent; share an end-vertex of the same label     
     */
    
    /**
     * Generate Compatibility Graph Nodes, where a node represents an edge from each of the 2 molecules
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraphNodes() throws IOException {

    	graphNodes.clear();
        //compGraphNodes.clear();

        IAtomContainer reactant = source;
        IAtomContainer product = target;
        
        for( int a = 0; a < reactant.getBondCount(); a++ ) {
        	
        	for( int b = 0; b < product.getBondCount(); b++ ) {
        		
        		// build aromatic (substituent) pairs - only aromatics can match other aromatics so such pairs can be removed from consideration 
        		// after being used.
        		
        		IBond rBond = reactant.getBond(a);
        		IBond pBond = product.getBond(b);
        		if( nodesMatch( reactant, rBond, product, pBond, shouldMatchBonds, atomTypeProperty ) 	
        				) {
        			
        			graphNodes.add( new int[]{ a, b } );
        		}
        		
        	}
        }

        return 0;
    }
    
    
    
    /**
     * Generate Compatibility Graph Nodes, where a node represents an edge from each of the 2 molecules
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraphNodesRaymondHeuristics() throws IOException {

    	graphNodes.clear();
        //compGraphNodes.clear();

        IAtomContainer reactant = source;
        IAtomContainer product = target;
        
        ringNodeMatching();
        
        for( int a = 0; a < reactant.getBondCount(); a++ ) {
        	
        	for( int b = 0; b < product.getBondCount(); b++ ) {
        		
        		// build aromatic (substituent) pairs - only aromatics can match other aromatics so such pairs can be removed from consideration 
        		// after being used.
        		
        		IBond rBond = reactant.getBond(a);
        		IBond pBond = product.getBond(b);
        		if( nodesMatch( reactant, rBond, product, pBond, shouldMatchBonds, atomTypeProperty ) 
        			&& ! ( rBond.getFlag( CDKConstants.ISAROMATIC ) & pBond.getFlag( CDKConstants.ISAROMATIC ) )	
        				) {
        			
        			graphNodes.add( new int[]{ a, b } );
        		}
        		
        	}
        }

        
        // need to sort nodes on ascending order of first element
        // TODO  Only needed for RASCAL AFAIK.  Needs a workaround
        Collections.sort( graphNodes, pairComparator );
     		
        return 0;
    }

    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraphEdges() throws IOException {

    
        // initialise with dummies
        //int largestGraphSize = Math.max( source.getAtomCount(), target.getAtomCount() );
    	int nodesCount = graphNodes.size();
        
        for (int a = 0; a < nodesCount; a++ ) {
        	
        	for (int b = a + 1; b < nodesCount; b++ ) {
        		
        		EdgeType et = edgesMatch(
                		source, 
                		target, 
                		source.getBond( graphNodes.get(a)[0] ), 
                		target.getBond( graphNodes.get(a)[1] ), 
                		source.getBond( graphNodes.get(b)[0] ), 
                		target.getBond( graphNodes.get(b)[1] ),
                		true
                	);
        		boolean matches = (et == EdgeType.CEDGE || et == EdgeType.DEDGE);
                
                
                // label/weight matching
                if( ! matches )
                	continue;
                
                
                if( useTopologicalDistance ) {
	                int topoDifference = Math.abs( pathDistancesHsMol[ graphNodes.get(a)[0] ][ graphNodes.get(b)[0] ] - pathDistancesQMol[ graphNodes.get(a)[1] ][ graphNodes.get(b)[1] ] );
	                
	        		// special case for topological distance constraint thing (in general)
	                if ( topoDifference > topologicalDistanceLimit ) {
	                	continue;
	                }
                }
                
                
                
                addEdge( a, b, et );
        	}
        }
        
        
        // connected common subgraph constraint goes here
      /*  for (int a = 0; a < edges.size(); a++ ) {
        	
        	boolean violates = true;
        	IBond bondA1 = source.getBond( graphNodes.get(a)[0] );
        	IBond bondA2 = target.getBond( graphNodes.get(a)[1] );
        	
        	
        	for (int b = a + 1; b < edges.size(); b++ ) {
        		
        		
        		IBond bondB1 = source.getBond( graphNodes.get(b)[0] );
            	IBond bondB2 = target.getBond( graphNodes.get(b)[1] );
  
        		// valid if at least one pair of connected bonds has a matching pair of connected bonds
        		if( ! (bondA1.isConnectedTo(bondB1) && bondA2.isConnectedTo(bondB2)) ) {
        			violates = false;
        			edges.get(a).remove( new Integer(b) );
        			//break;
        		}
        	}
        	
        	if( violates ) {
        		//edges.get(a).clear();
        		//System.err.println("violates");
        		//edges.get(b).clear();
        		for (int b = a + 1; b < graphNodes.size(); b++ ) {
	        		 edges.get(a).remove(b);  
	                 edges.get(b).remove(a);
        		}
        	}
        }*/
        
        return 0;
    }
    
    
    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraphEdgesRaymondHeuristics() throws IOException {
        
    	int nodesCount = graphNodes.size();
    	
    	/*int arCount = 0;
    	for( int n=0; n < target.getBondCount(); n++ ) {
    		if( ConvenienceTools.isAromatic( target.getBond(n) ) ) {
    			arCount++;
    		}
    	}
    	System.out.println( "arCount - " + arCount );*/
    	
        // initialise with dummies
        //int largestGraphSize = Math.max( source.getAtomCount(), target.getAtomCount() );
        for (int a = 0; a < graphNodes.size(); a++ ) {
        	edges.add( ConvenienceTools.createCollection( nodesCount ) );
        }
        
        for (int a = 0; a < graphNodes.size(); a++ ) {
        	
        	for (int b = a + 1; b < graphNodes.size(); b++ ) {
        		
        		/*
        		IBond reactantBond = null;
                IBond productBond = null;

                reactantBond = source.getBond(source.getAtom( graphNodes.get(a)[0] ), source.getAtom( graphNodes.get(b)[0] ));
                productBond = target.getBond(target.getAtom( graphNodes.get(a)[1] ), target.getAtom( graphNodes.get(b)[1] ));
                */
        		EdgeType et = edgesMatch(
                		source, 
                		target, 
                		source.getBond( graphNodes.get(a)[0] ), 
                		target.getBond( graphNodes.get(a)[1] ), 
                		source.getBond( graphNodes.get(b)[0] ), 
                		target.getBond( graphNodes.get(b)[1] ),
                		true
                	);

        		boolean matches = (et == EdgeType.CEDGE || et == EdgeType.DEDGE);
                
                // label/weight matching
                if( ! matches )
                	continue;
                
                
                int topoDifference = Math.abs( pathDistancesHsMol[ graphNodes.get(a)[0] ][ graphNodes.get(b)[0] ] - pathDistancesQMol[ graphNodes.get(a)[1] ][ graphNodes.get(b)[1] ] );
                
                // special case for topological distance constraint thing (in general)
        		if( useTopologicalDistance ) {
	                if ( topoDifference > topologicalDistanceLimit ) {
	                    continue;
	                }
        		}
        		
        		int maxRingTopoDistance = 0;
        		
        		List<Integer> source1Rings = sourceBondRings.get( graphNodes.get(a)[0] );
            	List<Integer> target1Rings = targetBondRings.get( graphNodes.get(a)[1] );
            	List<Integer> source2Rings = sourceBondRings.get( graphNodes.get(b)[0] );
            	List<Integer> target2Rings = targetBondRings.get( graphNodes.get(b)[1] );
            	
            
                
            	/*int minSRingSize = Math.min( sourceRingPaths.get(source1Rings.get(0)).size(), sourceRingPaths.get(source1Rings.get(1)).size() );
            	int minTRingSize = Math.min( targetRingPaths.get(target1Rings.get(0)).size(), targetRingPaths.get(target1Rings.get(1)).size() );
                */
        		
        		
        		
        		
        		
        		
        	
               
        		

        		/*
        		 * Divide molecules into rigid and flexible fragments
        		 * 
        		 * Rigid fragment heuristic - for a pair of mod prod nodes, if both edges land on the same rigid fragment per molecule,
        		 * the topological distance difference must < threshold
        		 * 
        		 * threshold is defined as the difference in size between the (single) rings being compared.  With fused rings
        		 * it's the two single rings of the fused ring being compared.  This should mean that same ring sizes prove to 
        		 * have the strictest constraint
        		 * 
        		 */
        		
        		
        		
        		
        		/* XXX 
        		 * "flexible" induced subgraph (topological distance) heuristic
        		 * 
        		 * Similar to CC induced subgraph heuristic - if the edge pairs are in the same chain in both molecules, apply
        		 * a topological distance constraint.  Instead of same path length however (topological difference of 0), the max topological
        		 * difference is the size difference between the 2 chains.
        		 */
        		Integer subgraphFIdS1 = sourceFlexibleInducedSubgraphBondIDs.get( graphNodes.get(a)[0] );
        		Integer subgraphFIdT1 = targetFlexibleInducedSubgraphBondIDs.get( graphNodes.get(a)[1] );
        		Integer subgraphFIdS2 = sourceFlexibleInducedSubgraphBondIDs.get( graphNodes.get(b)[0] );
        		Integer subgraphFIdT2 = targetFlexibleInducedSubgraphBondIDs.get( graphNodes.get(b)[1] );
        		
        		// all 4 bonds are in a "flexible" induced subgraph (only test source and target to see if they exist - non existence is -1)
        		if( subgraphFIdS1 >= 0 && subgraphFIdT1 >= 0  ) {
        			if( 
        					subgraphFIdS1 == subgraphFIdS2 && 
        					subgraphFIdT1 == subgraphFIdT2
        			) {
        				int minChainSize = Math.min( sourceFlexibleInducedSubgraphPathDistances[subgraphFIdS1] , targetFlexibleInducedSubgraphPathDistances[subgraphFIdT1] ) - 1;
        				int flexTopoLimit = Math.min( minChainSize, Math.abs(sourceFlexibleInducedSubgraphPathDistances[subgraphFIdS1] - targetFlexibleInducedSubgraphPathDistances[subgraphFIdT1]) );
        				if( topoDifference > flexTopoLimit ) {
        					//System.out.println("Flexible induced subgraph test failed " + flexTopoLimit);
        					continue;
        				}
        			}
        		}
        		
        		
        		
                // weak ring heuristic - source & target edges must be in same ring, with "or" rule on being fused bonds
                if( 
                		source1Rings != null
                		&& target1Rings != null
                		&& source2Rings != null
                		&& target2Rings != null
                    ) {
                	
                	
                	 /*XXX
                	 * additional constraint for strong ring heuristic - path length shouldn't be 0 if the rings are of different sizes.
                	 * 
                	 * Accounts for fused bonds too - looks for largest difference in size between the rings on the bond.
                	 * 
                	 */
                	for( int sr1=0; sr1 < source1Rings.size(); sr1++ ) {
                		for( int sr2=0; sr2 < source2Rings.size(); sr2++ ) {
                			for( int tr1=0; tr1 < target1Rings.size(); tr1++ ) {
                				for( int tr2=0; tr2 < target2Rings.size(); tr2++ ) {
                					
                					int rtd = Math.abs(
                                			Math.min( sourceRingBonds.get( source1Rings.get(sr1) ).size(), sourceRingBonds.get( source2Rings.get(sr2) ).size() ) -
                                			Math.min( targetRingBonds.get( target1Rings.get(tr1) ).size(), targetRingBonds.get( target2Rings.get(tr2) ).size() )
                                	);
                					
                					if( rtd > maxRingTopoDistance )
                						maxRingTopoDistance = rtd;
                				}
                			}
                		}
                	
                	}
                	
                	if( ringHeuristics ) {
	                	List<Integer> sourceRingsCommon = new ArrayList<Integer>( source1Rings );  
	                	sourceRingsCommon.retainAll( source2Rings );  // rings common between the two bonds in the molecule (i.e. check if in same ring)
	                	List<Integer> targetRingsCommon = new ArrayList<Integer>( target1Rings );
	                	targetRingsCommon.retainAll( target2Rings );
	                	
	                	boolean sourceSameRing = (sourceRingsCommon.size() > 0);
	                	boolean targetSameRing = (targetRingsCommon.size() > 0);
	                	
	                	// at least 1 bond must not be a fusion junction
	                	boolean sourceUnfusedCondition = (source1Rings.size() == 1) || (source2Rings.size() == 1);  
	                	boolean targetUnfusedCondition = (target1Rings.size() == 1) || (target2Rings.size() == 1);  
	                	
	                	
	                	if( sourceSameRing && sourceUnfusedCondition && targetUnfusedCondition ) {
	                		if( ! targetSameRing ) {
	                			//System.out.println( "weak ring test failed");
	                			continue;  // don't add bond
	                		} else {
	                			
	                			// passed weak ring heuristic
	                    		// so apply strong ring heuristic - aromatic ring bonds must have same path length
	                    		if(
	                    				ConvenienceTools.isAromatic( source.getBond( graphNodes.get(a)[0] ) )
	                            		&& ConvenienceTools.isAromatic( target.getBond( graphNodes.get(a)[1] ) )
	                            		&& ConvenienceTools.isAromatic( source.getBond( graphNodes.get(b)[0] ) )
	                            		&& ConvenienceTools.isAromatic( target.getBond( graphNodes.get(b)[1] ) )
	                    		) {
	                    			if( topoDifference > maxRingTopoDistance ) {
	                    				//System.out.println("strong ring test failed");
	                    				continue;
	                    			}
	                    		}
	                		}
	                	}
	
	                } 
                }
                
                
                /* XXX
                 *  "rigid" induced subgraph (topological distance) heuristic
                 *  
                 *  Acts on non-rotatable bonds.  Works in same way as CC-induced heuristic
                 */
        		Integer subgraphRIdS1 = sourceRigidInducedSubgraphBondIDs.get( graphNodes.get(a)[0] );
        		Integer subgraphRIdT1 = targetRigidInducedSubgraphBondIDs.get( graphNodes.get(a)[1] );
        		Integer subgraphRIdS2 = sourceRigidInducedSubgraphBondIDs.get( graphNodes.get(b)[0] );
        		Integer subgraphRIdT2 = targetRigidInducedSubgraphBondIDs.get( graphNodes.get(b)[1] );
        		
        		// all 4 bonds are in a "rigid" induced subgraph
        		if( subgraphRIdS1 >= 0 && subgraphRIdT1 >= 0 ) {
        			if( 
        					subgraphRIdS1.equals(subgraphRIdS2) && 
        					subgraphRIdT1.equals(subgraphRIdT2)
        			) {
        				if( topoDifference > maxRingTopoDistance ) {
        				     /*if( 
        		                		source1Rings != null
        		                		&& target1Rings != null
        		                		&& source2Rings != null
        		                		&& target2Rings != null
        		                    ) {
        					System.out.println("Rigid induced subgraph test failed " + maxRingTopoDistance + " - " + 
        							source1Rings.get(0) + " " + source2Rings.get(0) + " " + target1Rings.get(0) + " " + target2Rings.get(0) + " "  );
        				     }*/
        					continue;
        				}
        			}
        		}
                
             	/* XXX 
             	 * C-C induced subgraph (topological distance) heuristic.
             	 * 
             	 * Note that the topological difference limit is no longer 0 like in the original paper.  This only differs from 0 when
             	 * 2 rings of different sizes are being compared. 
             	 */
        		Integer subgraphIdS1 = sourceCCInducedSubgraphBondIDs.get( graphNodes.get(a)[0] );
        		Integer subgraphIdT1 = targetCCInducedSubgraphBondIDs.get( graphNodes.get(a)[1] );
        		Integer subgraphIdS2 = sourceCCInducedSubgraphBondIDs.get( graphNodes.get(b)[0] );
        		Integer subgraphIdT2 = targetCCInducedSubgraphBondIDs.get( graphNodes.get(b)[1] );
        		
        		// all 4 bonds are in a C-C induced subgraph
        		if( subgraphIdS1 >= 0 && subgraphIdT1 >= 0 ) {
        			if( 
        					subgraphIdS1.equals(subgraphIdS2) && 
        					subgraphIdT1.equals(subgraphIdT2)
        			) {
        				if( topoDifference > maxRingTopoDistance ) {
        					//System.out.println("C-C induced subgraph test failed");
        					continue;
        				}
        			}
        		}
                
        		
        		/* XXX
        		 * Additional hack when comparing molecules with several flexible chains.  
        		 * 
        		 * If the edge pairs aren't in rings, then the topological distance limit is equal to the largest chain length
        		 * in either molecule.
        		 */
                if( source1Rings == null
                		&& target1Rings == null
                		&& source2Rings == null
                		&& target2Rings == null ) {
	                if( topoDifference > flexTopoDistLim ) {
	        				//System.out.println("non-ring topological constraint failed " + flexTopoDistLim );
	        				continue;
	        		}
                }
                
                
        		
                
                addEdge( a, b, et );
        		
        		
        		
        	}
        }
        
        //System.out.println("non-ring topological constraint failed " + flexTopoDistLim );
        
        return 0;
    }


    
    
    
    /**
     * Remember that a node in this modular product represents a pair of molecule bonds, NOT atoms!
     * 
     * Thus, the bond types of the matched nodes must be compatible, as well as the atom labels
     * 
     * @param ac1
     * @param bondA1	the bond from the first molecule
     * @param ac2
     * @param bondA2	the bond from the other molecule
     * @param shouldMatchBonds
     * @return
     */
    public static boolean nodesMatch(IAtomContainer ac1,
            IBond bondA1,
            IAtomContainer ac2,
            IBond bondA2,
            boolean shouldMatchBonds,
            String atomTypeProperty ) {
    	
    	
    	/* 
         * aromatic bond hack - aromatic bonds should not match non aromatic bonds
         * 
         * I don't know why the native CDK bond matcher doesn't account for this because the source
         * code seems to indicate that it does?
         */

    	
        if( shouldMatchBonds && 
        	!(bondA1 instanceof IQueryBond) && 
        	ConvenienceTools.isAromatic(bondA1) ^ ConvenienceTools.isAromatic(bondA2) )
        	return false;

        //Bond Matcher
    	DefaultBondMatcher bondMatcher = null;
        
        if( bondA1 instanceof IQueryBond && shouldMatchBonds )
        	bondMatcher = new DefaultBondMatcher( (IQueryBond) bondA1 );
        else 
        	bondMatcher = new DefaultBondMatcher( ac1, bondA1, shouldMatchBonds );
        
        //bondMatcher.setBondMatchFlag(shouldMatchBonds);
        
        if( shouldMatchBonds )
	        if( ! bondMatcher.matches(ac2, bondA2) )
	        	return false;
        
        
        
        // testing of custom atom type definitions
        if( atomTypeProperty != null ) {
        	
        	Object typeA1_0 = bondA1.getAtom(0).getProperty( atomTypeProperty );
        	Object typeA2_0 = bondA2.getAtom(0).getProperty( atomTypeProperty );
        	Object typeA1_1 = bondA1.getAtom(1).getProperty( atomTypeProperty );
        	Object typeA2_1 = bondA2.getAtom(1).getProperty( atomTypeProperty );
        	
        	if( 
        			typeA1_0.equals(typeA2_0) &&
        			typeA1_1.equals(typeA2_1) 
        	)
        		return true;
        	
        	if( 
        			typeA1_0.equals(typeA2_1) &&
        			typeA1_1.equals(typeA2_0) 
        	)
        		return true;
        	
        	
        	return false;
        }
        
        //Atom Matchers
        AtomMatcher atomMatcher1 = null;
        AtomMatcher atomMatcher2 = null;
        
        // What matches the first atom in the first bond? 
        if( bondA1.getAtom(0) instanceof IQueryAtom )
        	atomMatcher1 = new DefaultHyperstructureAtomMatcher( (IQueryAtom) bondA1.getAtom(0), (IQueryAtomContainer) ac1 );
        else
        	atomMatcher1 = new DefaultHyperstructureAtomMatcher(ac1, bondA1.getAtom(0), false); 
        
        // What matches the second atom in the first bond?
        if( bondA1.getAtom(1) instanceof IQueryAtom )
        	atomMatcher2 = new DefaultHyperstructureAtomMatcher( (IQueryAtom) bondA1.getAtom(1), (IQueryAtomContainer) ac1 );
        else
        	atomMatcher2 = new DefaultHyperstructureAtomMatcher(ac1, bondA1.getAtom(1), false);
        
        if( atomMatcher1.matches(ac2, bondA2.getAtom(0)) && atomMatcher2.matches(ac2, bondA2.getAtom(1)) )
        	return true;
        
        if( atomMatcher1.matches(ac2, bondA2.getAtom(1)) && atomMatcher2.matches(ac2, bondA2.getAtom(0)) )
    		return true;
        
        /*
        if( 
        	atomMatcher1.matches(ac2, bondA2.getAtom(0)) && atomMatcher2.matches(ac2, bondA2.getAtom(1))  // atom matching is bugged
        	//bondA1.getAtom(0).getSymbol().equals( bondA2.getAtom(0).getSymbol() ) && bondA1.getAtom(1).getSymbol().equals( bondA2.getAtom(1).getSymbol() )
        		)
        	return true;
        */
        return false;
    }
    
    
    /**
     * 
     * @param ac1
     * @param ac2
     * @param bondA1	edge 1 from the first pair (node in the modular product) - first molecule
     * @param bondA2	edge 2 from the first pair - second molecule
     * @param bondB1	edge 1 from the second pair - first molecule
     * @param bondB2
     * @return
     */
    public static EdgeType edgesMatch( 
    		IAtomContainer ac1,
    		IAtomContainer ac2,
    		IBond bondA1,
    		IBond bondA2,
    		IBond bondB1,
    		IBond bondB2,
    		boolean allowDEdges
    ) {
    	
    	// identity is redundant - we're not gonna connect any edges to the same vertex!
    	if( bondA1 == bondB1 || bondA2 == bondB2 )
    		return EdgeType.LOOP;
    	
    	// c-edges: true if a1 & a2 are adjacent, and b1 & b2 are also adjacent 
    	// (and that the common vertex is the same label as the other common vertex)
    	
    	// find common vertex in the first molecule
    	IAtom commonVertex1 = null;
    	if( bondA1.contains( bondB1.getAtom(0) ) ) {
    		commonVertex1 = bondB1.getAtom(0);
    	} else if( bondA1.contains( bondB1.getAtom(1) ) ) {
    		commonVertex1 = bondB1.getAtom(1);
    	}
    	
    	IAtom commonVertex2 = null;
    	if( bondA2.contains( bondB2.getAtom(0) ) ) {
    		commonVertex2 = bondB2.getAtom(0);
    	} else if( bondA2.contains( bondB2.getAtom(1) ) ) {
    		commonVertex2 = bondB2.getAtom(1);
    	}
    	
    	// connected common substructure constraint
    	//if( ! (bondA1.isConnectedTo(bondB1) && bondA2.isConnectedTo(bondB2)) )
    		//return false;
    	
    	if( commonVertex1 != null && commonVertex2 != null ) {
    		
	    	AtomMatcher atomMatcher = null;
	    	if( commonVertex1 instanceof IQueryAtom )
	        	atomMatcher = new DefaultHyperstructureAtomMatcher( (IQueryAtom) commonVertex1, (IQueryAtomContainer) ac1 );
	        else
	        	atomMatcher = new DefaultHyperstructureAtomMatcher( ac1, commonVertex1, true);
	    	
	    	if( atomMatcher.matches(ac2, commonVertex2) )
	    		return EdgeType.CEDGE;  
	    	
    	} else if(  (commonVertex1 == null && commonVertex2 == null) ) {
    		return EdgeType.DEDGE; 
    	}
    	
    	
    	return EdgeType.REDEDGE;
    }
    
    
    
    protected void createEdgeDataStructures(int nodesCount) {
    	edges = new ArrayList<Collection<Integer>>( graphNodes.size() );
    	
        for (int a = 0; a < nodesCount; a++ ) {
        	//edges.add( new BitSetExtended<Integer>( largestGraphSize ) );
        	edges.add( ConvenienceTools.createCollection( nodesCount ) );
        }
    }
    
    protected void addEdge( int a, int b, EdgeType type ) {
    	
    	//if( (type == EdgeType.CEDGE || type == EdgeType.DEDGE) ) {
	    	edges.get(a).add(b);  
	        edges.get(b).add(a);
	        
	        numberOfEdges += 2;
    	//}
    }
    
    
    
    /*
     * Identify all benzene rings in both molecules, note their bond indices
     * search for common substituents between all aromatic ring pairs
     * 
     * go thru all benzene rings
     * - at least one is benzene (the other can be fused or non-benzoid):
     * 	- if any have common substituents, map an arbitrary path for each common substituent
     * 	- if not, just make a 1:1 mapping (if rings are of same size)
     * - otherwise (both are non benzene), map via normal mod prod rules
     * 
     * Start from each common substituent mapping, then keep adjacent pairs (topo distance of 0)
     * 
     * After this, benzenes with non-common substituents can be 1:1 matched
     * 
     * 
     * Then can do other bond matches (though if both are in aromatic rings, ignore them), outside this method
     * 
     */
    private void ringNodeMatching() {
    	

    	List<List<Integer>> sourceAromatics = new ArrayList<List<Integer>>( sourceRings.size() );
    	List<List<Integer>> targetAromatics = new ArrayList<List<Integer>>( targetRings.size() );
    	
    	List<Boolean> sourceBigFusedBonds = new ArrayList<Boolean>( sourceRings.size() );
    	List<Boolean> targetBigFusedBonds = new ArrayList<Boolean>( targetRings.size() );
    	
    	
    	/*
    	SSSRFinder sssrSource = new SSSRFinder(source);
    	for( IAtomContainer ring : sssrSource.findRelevantRings().atomContainers() ) {
    		ArrayList<Integer> ringBonds = new ArrayList<Integer>( 6 );
    		
    		for( IBond b : ring.bonds() ) {
    			//System.out.print( source.getBondNumber(b) + "," );
    			ringBonds.add( source.getBondNumber(b) );
    		}
    		//sourceRings.add( ringBonds );
    		//System.out.println();
    	}
    	
    	SSSRFinder sssrTarget = new SSSRFinder(target);
    	for( IAtomContainer ring : sssrTarget.findRelevantRings().atomContainers() ) {
    		ArrayList<Integer> ringBonds = new ArrayList<Integer>( 6 );
    		
    		for( IBond b : ring.bonds() ) {
    			//System.out.print( source.getBondNumber(b) + "," );
    			ringBonds.add( target.getBondNumber(b) );
    		}
    		//targetRings.add( ringBonds );
    		//System.out.println();
    	}*/
    	
    	//List<List<Integer>> sourceSubstituentBonds = new ArrayList<List<Integer>>( sourceAromatics.size() );
    	//List<List<Integer>> targetSubstituentBonds = new ArrayList<List<Integer>>( targetAromatics.size() );
    	ArrayList<Map<Integer, int[]>> sourceSubstituentsToRings = new ArrayList<Map<Integer, int[]>>( sourceAromatics.size() );
    	ArrayList<Map<Integer, int[]>> targetSubstituentsToRings = new ArrayList<Map<Integer, int[]>>( targetAromatics.size() );
    	
    	//sourceRings.remove(0);
    	//targetRings.remove(0);
    	
    	ArrayList<int[]> ringPairs = new ArrayList<int[]>();
    	
    	/* keep all aromatic ring bonds (some fused rings will have aromatic & non-aromatic rings within them)
    	 */
    	for( List<Integer> bondIndices : sourceRings ) {
    		
    		int arCount = 0;
    		
    		List<Integer> arRing = new ArrayList<Integer>( bondIndices.size() );
    		
    		for( Integer rBond : bondIndices ) {
    			if( source.getBond( rBond ).getFlag( CDKConstants.ISAROMATIC ) ) {
    				arRing.add( rBond );
    				++arCount;
    			}
    		}
    		
    		/*
    		 * Due to a bug in aromatic perception of the Daylight (and others as of 4/2/2015) method in CDK, some bonds in fused
    		 * rings are not marked as aromatic.  Thus, this accounts for that (i.e includes the whole ring when a certain threshold
    		 * of aromatic bonds is found)
    		 */
    		if( arCount >= bondIndices.size() - 1 )
    			arRing = bondIndices;
    		
    		
    		if( arRing.size() > 0 ) {
    			sourceAromatics.add( arRing );
    			
    			/*int[] f1 = new int[ arRing.size() ];
    	    	for( int n = 0; n < f1.length; n++ ) { f1[n] = n; }
    	    	sourceFusedBonds.add( findFusedBonds( source, f1 ) );*/
    		}
    		
    	
    	}
    	
    	for( List<Integer> bondIndices : targetRings ) {
    		
    		int arCount = 0;
    		
    		List<Integer> arRing = new ArrayList<Integer>( bondIndices.size() );
    		
    		for( Integer rBond : bondIndices ) {
    			if( target.getBond( rBond ).getFlag( CDKConstants.ISAROMATIC ) ) {
    				arRing.add( rBond );
    				++arCount;
    			}
    		}
    		
    		if( arCount >= bondIndices.size() - 1 )
    			arRing = bondIndices;
    		
    		if( arRing.size() > 0 ) {
    			targetAromatics.add( arRing );
    			
    			/*int[] f2 = new int[ arRing.size() ];
    	    	for( int n = 0; n < f2.length; n++ ) { f2[n] = n; }
    	    	targetFusedBonds.add( findFusedBonds( target, f2 ) );*/
    		}
    	}
    	
    	
    	// find "big" fused rings
    	for( List<Integer> bondIndices : sourceAromatics ) {
    		for( int index : bondIndices ) {
    			if( ((ArrayList) source.getBond(index).getProperty( CDKConstants.RING_CONNECTIONS )).size() > 1 ) {
    				sourceBigFusedBonds.add( true );
    				break;
    			}
    		}
    		
			sourceBigFusedBonds.add( false );
    	}
    	
    	// find "big" fused rings in target
    	for( List<Integer> bondIndices : targetAromatics ) {
    		for( int index : bondIndices ) {
    			
    			if(verbose)
    				System.out.println( "ring members " + target.getBond(index).getProperty( CDKConstants.RING_CONNECTIONS ) );
    			
    			if( ((ArrayList) target.getBond(index).getProperty( CDKConstants.RING_CONNECTIONS )).size() > 1 ) {
    				
    				targetBigFusedBonds.add(true);
    				break;
    			}
    		}
    		
			targetBigFusedBonds.add( false );
    	}
    	
    	
    	
    	// find ring paths (bond indices)
    	List<List<int[]>> sourceRingPaths = generateRingPaths( source, sourceAromatics );
    	List<List<int[]>> targetRingPaths = generateRingPaths( target, targetAromatics );
		

    	
    	// catalogue substituents around each ring (a substituent bond must not be in another ring)
    	for( int ri = 0; ri < sourceAromatics.size(); ri++ ) {
    		
    		//Map<Integer, int[]> subToRingBonds = getRingSubstituents( sourceAromatics.get(ri), source );
    		Map<Integer, int[]> subToRingBonds = getRingSubstituents( sourceRingPaths.get(ri), source );
    		
    		sourceSubstituentsToRings.add( subToRingBonds );
    	}
    	
    	
    	// catalogue substituents around each ring (a substituent bond must not be in another ring)
    	for( int ri = 0; ri < targetAromatics.size(); ri++ ) {
    		//sourceSubstituentsToRings.add( ri, new HashMap<Integer, int[]>() );
    		//Map<Integer, int[]> subToRingBonds = getRingSubstituents( targetAromatics.get(ri), target );
    		Map<Integer, int[]> subToRingBonds = getRingSubstituents( targetRingPaths.get(ri), target );
    		
    		targetSubstituentsToRings.add( subToRingBonds );
    	}
    	

    	
    	
    	
    	
    	if( verbose ) {
	    	System.out.println( "source rings - " +  sourceRings );
	    	System.out.println( "target rings - " + targetRings );
	    	ConvenienceTools.printList( sourceRingPaths );
	    	ConvenienceTools.printList( targetRingPaths );
	    	System.out.println( "source aromatic rings - " + sourceAromatics );
	    	System.out.println( "target aromatic rings - " + targetAromatics );
	    	System.out.println( "source fused - " + sourceBigFusedBonds );
	    	System.out.println( "target fused - " + targetBigFusedBonds );
	    	System.out.println( "s substituents " + sourceSubstituentsToRings );
	    	System.out.println( "t substituents " + targetSubstituentsToRings );
    	}
    	
    	
    	// now to go through all the rings and start performing matches
    	for( int sInd = 0 ; sInd < sourceAromatics.size(); sInd++ ) {
    		for( int tInd = 0 ; tInd < targetAromatics.size(); tInd++ ) {
    			
    			
    			
    			ArrayList<int[]> arbitraryPairs = new ArrayList<int[]>();
    			ArrayList<int[]> nonSubPairs = new ArrayList<int[]>( sourceAromatics.get(sInd).size() );
				
				int maxSubstituentRingSize = 0;
				ArrayList<int[]> substituentPairs = new ArrayList<int[]>( sourceAromatics.get(sInd).size() );
				
				if( verbose ) {
					StringBuilder sourceString = new StringBuilder(40);
					for( Integer n : sourceAromatics.get(sInd) ) {
						IBond b = source.getBond(n);
						sourceString.append( b.getAtom(0).getSymbol() );
						//sourceString.append( b.getOrder().toString() );
						sourceString.append( b.getAtom(1).getSymbol() );
						sourceString.append(" ");
					}
					
					StringBuilder targetString = new StringBuilder(40);
					for( Integer n : targetAromatics.get(tInd) ) {
						IBond b = target.getBond(n);
						targetString.append( b.getAtom(0).getSymbol() );
						//sourceString.append( b.getOrder().toString() );
						targetString.append( b.getAtom(1).getSymbol() );
						targetString.append(" ");
					}
					
					System.out.println("source " + sInd + " " + sourceString + " | target " + tInd + " " + targetString   );
					System.out.println("fused - " + (sourceBigFusedBonds.get(sInd) & targetBigFusedBonds.get(tInd)) );
				}
				
				/* if both rings are "non-benzoid", avoid using heuristics
				 * 
				 * non-benzoid - either a ring with an odd number of bonds, or a fused ring
				 */
				boolean sourceNonBenzoid = ( sourceAromatics.get(sInd).size() % 2 != 0 || sourceBigFusedBonds.get(sInd) );
				boolean targetNonBenzoid = ( targetAromatics.get(tInd).size() % 2 != 0 || targetBigFusedBonds.get(tInd) );
				boolean bothNonBenzoid = sourceNonBenzoid && targetNonBenzoid;
				
				if( verbose )
					System.out.println( "ne " + bothNonBenzoid + " " + sourceAromatics.get(sInd).size() + " " + sourceBigFusedBonds.get(sInd)
						+ " " + targetAromatics.get(tInd).size() + " " + targetBigFusedBonds.get(tInd) );
				
    			// if both are fused/non-benzoid, skip
    			if( ! bothNonBenzoid && (sourceBigFusedBonds.get(sInd) & targetBigFusedBonds.get(tInd)) == false ) {
    				
    				
    				 // flag if either of the rings is symmetrical - if 1 is, it doesn't matter if
    				  // 2 is so don't bother repeating the check
    				  boolean srSymmetric = false; 
    				  boolean trSymmetric = false;
    				  List<int[]> sParaSubs = getParaSubs( pathDistancesHsMol, sourceSubstituentsToRings.get(sInd), sourceAromatics.get(sInd).size() );
    				  List<int[]> tParaSubs = getParaSubs( pathDistancesQMol, targetSubstituentsToRings.get(tInd), targetAromatics.get(tInd).size() );
    				  
    				  if( sourceSubstituentsToRings.get(sInd).size() == 1 ||
    				      ( sourceSubstituentsToRings.get(sInd).size() == 2 && sParaSubs != null ) ) {
    					  srSymmetric = true;
    				  }
    				  if( ! srSymmetric && ( targetSubstituentsToRings.get(tInd).size() == 1 ||
    				                       ( targetSubstituentsToRings.get(tInd).size() == 2 && tParaSubs != null ) ) ) {
    					  trSymmetric = true;
    				  }

    			
	    			// see if there're any (common) substituents
	    			if( sourceSubstituentsToRings.get(sInd).size() > 0 && targetSubstituentsToRings.get(tInd).size() > 0 ) {
	    				
	    				/* construct a list of source-target substitution pairs
	    				 * 
	    				 * Importance of this will become apparent when we start dealing with para-subbed rings 
	    				 */
	    				List<int[]> substituentPairList = new ArrayList<int[]>( sourceSubstituentsToRings.size() * targetSubstituentsToRings.size() );
	    				
	    				/*
    					 * Deal with para-substituted ring substituents first - idea that if both rings have common
    					 * para-substituents, then all other substituents can be ignored
    					 */
    					if(  sParaSubs != null && tParaSubs != null && (sourceSubstituentsToRings.get(sInd).size() == 2 || targetSubstituentsToRings.get(tInd).size() == 2 ) ) {
    						if( verbose ) {
	    						System.out.println( "sps " +  sParaSubs );
	    						System.out.println( "tps " +  tParaSubs );
    						}
    						
    						for( int[] spSubPair : sParaSubs ) {
    							
    								IBond spSubBond = source.getBond(spSubPair[0]);
    								
	    							for( int[] tpSubPair : tParaSubs ) {
	    								
	    									IBond tpSubBond = target.getBond(tpSubPair[0]);
	    									
	    									if( nodesMatch(source, spSubBond, target, tpSubBond, shouldMatchBonds, atomTypeProperty) ) {
	    										substituentPairList.add( new int[]{ spSubPair[0], tpSubPair[0] } );
	    										
	    										if( verbose )
	    											System.out.println("sub pair - " + spSubPair[0] + "," + tpSubPair[0] + "  " + spSubBond.getAtom(0).getSymbol() + spSubBond.getAtom(1).getSymbol() );
	    									}
	    								
	    								
	    							}
    							
    							
    						}
    					} else {
	    				
		    				for( Integer sSubInd : sourceSubstituentsToRings.get(sInd).keySet() ) {
		    					for( Integer tSubInd : targetSubstituentsToRings.get(tInd).keySet() ) {
		    						substituentPairList.add( new int[]{ sSubInd, tSubInd } );
		    					}
		    				}
	    				
    					}
	    				
	    				for( int[] sPair : substituentPairList ) {
	    					int sSubInd = sPair[0];
	    					int tSubInd = sPair[1];
	    					
	    					IBond sSub = source.getBond( sSubInd );
	    					
	    					
	    						
	    						HashSet<Integer> currSubstituentRingSize = new HashSet<Integer>( sourceAromatics.get(sInd).size() );
	    						
	    						IBond tSub = target.getBond( tSubInd );
	    						
	    						// perform ring matching here if substituents are common
	    		    			if( nodesMatch(source, sSub, target, tSub,  shouldMatchBonds, atomTypeProperty) ) {
	    		    				
	    		    				/*
	    		    				 * map all compatible paths - based on the substituents (as starting points for path enumeration)
	    		    				 */
	    		    				
	    		    				if( verbose )
	    		    					System.out.println( sInd + " " + tInd + " common substitution " + sSubInd + " vs " + tSubInd  );
	    		    				
	    		    				List<int[]> currentPairs = null;
	    		    				
	    		    				
	    		    				
	    		    					
	    		    					currentPairs = new ArrayList<int[]>( sourceRingPaths.get(sInd).get(0).length * 2 );
	    		    					
	    		    					
	    		    					
	    		    					
	    		    					/* as we're mapping single rings to fused rings, fused rings may have more than
	    		    					 * one path, hence the nested for loop.  This' the idea that we map a single ring
	    		    					 * to one of the ring components of the fused ring
	    		    					 */
	    		    					for( int sPathIndex = 0; sPathIndex < sourceRingPaths.get(sInd).size(); sPathIndex++ ) {
	    		    						for( int tPathIndex = 0; tPathIndex < targetRingPaths.get(tInd).size(); tPathIndex++ ) {
			    		    					int[] srPath = sourceRingPaths.get(sInd).get(sPathIndex);
			    		    					int[] trPath = targetRingPaths.get(tInd).get(tPathIndex);
			    		    					
			    		    					int sStart = -1; 
			    		    					int tStart = -1;
			    		    					
			    		    					/*
		    		    						 * Tests for aligning rings correctly via paths:
		    		    						 * 
		    		    						 * take first ring, find sStart (a bond on the path that's adjacent to the substitution).
		    		    						 * Find the path distance (ring-wise) from sStart to the other bond adjacent to substitution
		    		    						 * 
		    		    						 * now for the next ring, we need the tStart whose distance to other sub bond is identical.
		    		    						 * 4 tests to find identical match:
		    		    						 * - first bond found, loop forward along path to other bond
		    		    						 * - first bond found, loop reverse direction to other bond
		    		    						 * - second bond found, loop forward to 1st bond
		    		    						 * - second bond found, loop reverse to 1st bond
		    		    						 * 
		    		    						 * if no identical solution found, choose arbitrary start from one of 2 bonds
		    		    						 * 
		    		    						 * return int[]{ startBond, distance }
		    		    						 * 
		    		    						 */
			    		    					
			    		    					int sDistance = 0;
			    		    					//int tDistance = 0;
			    		    					int sStart2 = -1;
			    		    					int tStart2 = -1;
			    		    					boolean reversePaths = false;
			    		    					
			    		    					//int sStart = 0, tStart = 0;
			    		    					for( int s = 0; s < srPath.length; s++ ) {
	
			    		    						if( sourceSubstituentsToRings.get(sInd).get( sSubInd )[0] == srPath[s]  ) {
			    		    							sStart = s;
			    		    						}
			    		    						
			    		    						if( sourceSubstituentsToRings.get(sInd).get( sSubInd )[1] == srPath[s]  ) {
			    		    							sStart2 = s;
			    		    						}
			    		    						
			    		    						if( sStart >= 0 && sStart2 >= 0 ) {
			    		    							sDistance = ringPathDistance(srPath, sStart, sStart2, false);
			    		    							
			    		    							if( verbose )
			    		    								System.out.println( sPathIndex + " s start1 " + sStart + " " + sDistance );
			    		    							
			    		    							break;
			    		    						}
			    		    					}
			    		    					
			    		    					
			    		    					
			    		    					
			    		    					for( int t = 0; t < trPath.length; t++ ) {
	
			    		    						if( targetSubstituentsToRings.get(tInd).get( tSubInd )[0] == trPath[t]  ) {
			    		    							tStart = t;
			    		    							
			    		    							if( srPath.length != trPath.length )
			    		    								break;
			    		    						}
			    		    						
			    		    						if( targetSubstituentsToRings.get(tInd).get( tSubInd )[1] == trPath[t]  ) {
			    		    							tStart2 = t;
			    		    						}
			    		    						
			    		    						if( tStart >= 0 && tStart2 >= 0 ) {
			    		    							int tDistance1 = ringPathDistance(trPath, tStart, tStart2, false);
			    		    							int tDistance2 = ringPathDistance(trPath, tStart, tStart2, true);
			    		    							
			    		    							if( tDistance1 == sDistance )
			    		    								reversePaths = false;
			    		    							else if( tDistance2 == sDistance ) {
			    		    								reversePaths = true;
			    		    							}
			    		    							
			    		    							if( verbose )
			    		    								System.out.println( tPathIndex + " t start1 " + tStart + " " + tDistance1 + " " + tDistance2 + " " + reversePaths );
			    		    							
			    		    							break;
			    		    						}
			    		    					}
			    		    		
			    		    					

			    		    					
			    		    					/* if one is -1, then there're no common substituents
			    		    					 * 
			    		    					 * However, we still attempt to match the rings anyway
			    		    					 */
			    		    					if( sStart == -1 )
			    		    						sStart = 0;
			    		    					
			    		    					if( tStart == -1 )
			    		    						tStart = 0;
			    		    					
			    		    					
				    		    				currentPairs.addAll( matchRingsViaPaths( srPath, trPath, sStart, tStart, reversePaths ) );
				    		    				//currentPairs.addAll( matchRingsViaPaths( srPath, trPath, sStart2, tStart2, !reversePaths ) );	
				    		    				
				    		    				// need reverse mapping too, due to ambiguity from starting points
				    		    				if( ! srSymmetric && ! trSymmetric ) {
				    		    					currentPairs.addAll( matchRingsViaPaths( srPath, trPath, sStart, tStart, ! reversePaths ) );
				    		    					//currentPairs.addAll( matchRingsViaPaths( srPath, trPath, sStart2, tStart2,  reversePaths ) );
				    		    				}
			    		    					
	    		    						}
	    		    					}
	    		    					//currentPairs.add( new int[]{ sSubInd, tSubInd } );
	    		    				
	    		    				
		    		    			for( int[] pair : currentPairs ) {
		    		    				currSubstituentRingSize.add( pair[0] );
		    		    				substituentPairs.add( pair );
		    		    			}
		    		    			
		    		    			
		    		    			if( currSubstituentRingSize.size() > maxSubstituentRingSize ) {
		    		    				maxSubstituentRingSize = currSubstituentRingSize.size();
		    		    			}
	    					}
	    				}
	    			}
    			
	    			
    				
    				/* In addition to common substituents we attempt to match hetero-atoms (provided both rings are not fused)
    				 * 
    				 * heteroatom bonds go first - see if there're common heteroatoms.  If so, match from the first common.  Traverse one ring in both directions
    				 * 
    				 * if no common heteroatoms, start from a non-common pair (but don't match them), one ring in BOTH directions
    				 * if no heteroatoms just choose arbitrary start point.  One direction only
    				*/

	    			if( (sourceBigFusedBonds.get(sInd) & targetBigFusedBonds.get(tInd)) == false ) {
    			
	    				// catalogue heteroatom bonds (a non CC or CH bond)
	    				ArrayList<Integer> sHets = new ArrayList<Integer>(6);  // indices of heteroAtoms in the ring
	    				ArrayList<Integer> tHets = new ArrayList<Integer>(6);
	    				
	    				
						int[] srPath = sourceRingPaths.get(sInd).get(0);
    					int[] trPath = targetRingPaths.get(tInd).get(0);
	    				
	    				for( int i = 0; i < srPath.length; i += 2 ) {
	    					IBond b = source.getBond( srPath[i] );
	    					
	    					for( IAtom a : b.atoms() ) {
	    						if( a.getAtomicNumber() != 6 ) {
	    							sHets.add( i );
	    							break;
	    						}
	    					}
	    					
	    					
	    				}
	    				
	    				for( int i = 0; i < trPath.length; i++ ) {
	    					IBond b = target.getBond( trPath[i] );
	    					
	    					for( IAtom a : b.atoms() ) {
	    						if( a.getAtomicNumber() != 6 ) {
	    							tHets.add( i );
	    							break;
	    						}
	    					}

	    				}
	    				
	    				if( verbose ) {
	    					System.out.println( sourceAromatics.get(sInd).size() + " " + targetAromatics.get(tInd).size() + " " + sourceBigFusedBonds.get(sInd) + " " + targetBigFusedBonds.get(tInd) + " " + (sourceBigFusedBonds.get(sInd) & targetBigFusedBonds.get(tInd)) );
		    				System.out.println( "sHets - " + sHets);
		    				System.out.println( "tHets - " + tHets);
	    				}
	    				
	    				ArrayList<int[]> hetMatchPairs = new ArrayList<int[]>( sourceAromatics.get(sInd).size() );
	
	    				// find common heteroatom bonds
	    				for( Integer sHet : sHets ) {
	    					//IBond sHetBond = source.getBond(sHet);
	    					for( Integer tHet : tHets ) {
	    						//IBond tHetBond = target.getBond(tHet);
	    						
	    						
	    						
	    						if( verbose )
	    							System.out.println( sInd + " " + tInd + " Heteroatom matching " + sHet + " " + tHet );
	    						
	    						
	    						
								/*hetMatchPairs = matchRingsViaRefPoint(
										sourceAromatics.get(sInd), targetAromatics.get(tInd), 
										sourceAromatics.get(sInd).get(sHet), targetAromatics.get(tInd).get(tHet) 
								);*/
								
		    					hetMatchPairs.addAll( matchRingsViaPaths( srPath, trPath, sHet, tHet, false ) );
			    				//currentPairs.addAll( matchRingsViaPaths( srPath, trPath, sStart2, tStart2, !reversePaths ) );	
			    				
			    				/*// need reverse mapping too, due to ambiguity from starting points
			    				if( ! srSymmetric && ! trSymmetric ) {
			    					currentPairs.addAll( matchRingsViaPaths( srPath, trPath, sStart, tStart, ! reversePaths ) );
			    					//currentPairs.addAll( matchRingsViaPaths( srPath, trPath, sStart2, tStart2,  reversePaths ) );
			    				}*/
								
			    				
			    				if( hetMatchPairs.size() > nonSubPairs.size() )
			    					nonSubPairs = hetMatchPairs;
	
	    					}
	    				}
	    			}
    				 
	    			// perform arbitrary match if the other 2 stages fail
    				if( nonSubPairs.size() == 0 && substituentPairs.size() == 0 ) {
	    				arbitraryPairs = new ArrayList<int[]>( sourceAromatics.get(sInd).size() );
	    				
	    				if( verbose )
	    					System.out.println( sInd + " " + tInd + " non fused arbitrary matching");
	    				
	    				int sStart = 0; 
	    				int tStart = 0;
	    				
	    				// now to perform the arbitrary match
	    				for( int bondInd = 0; bondInd < Math.min( sourceAromatics.get(sInd).size(), targetAromatics.get(tInd).size() ); bondInd++ ) {
	    					IBond sBond = source.getBond( sourceAromatics.get(sInd).get(sStart) );
	    					IBond tBond = target.getBond( targetAromatics.get(tInd).get(tStart) );
	    					
	    					if( nodesMatch(source, sBond, target, tBond, shouldMatchBonds, atomTypeProperty) ) {
	    						arbitraryPairs.add( new int[]{ sourceAromatics.get(sInd).get(sStart), targetAromatics.get(tInd).get(tStart) } );
	    					}
	    					
	    					if( sStart + 1 >= sourceAromatics.get(sInd).size() )
	    						sStart = 0;
	    					else
	    						sStart++;
	    					
	    					if( tStart + 1 >= targetAromatics.get(tInd).size() )
	    						tStart = 0;
	    					else
	    						tStart++;
	    				}
    				}
    			} else {
    				// seeing as the symmetry conditions have been violated we won't perform any node deletions
    				// i.e. standard modular product node matching
    				
    				if( verbose )
    					System.out.println("both non-benzoid - performing arbitrary match");
    				
					for( Integer a : sourceAromatics.get(sInd) ) {
			        	for( Integer b : targetAromatics.get(tInd) ) {
			        		
			        		IBond rBond = source.getBond(a);
			        		IBond pBond = target.getBond(b);
			        		
			        		if( nodesMatch( source, rBond, target, pBond, shouldMatchBonds, atomTypeProperty ) ) {
			        			arbitraryPairs.add( new int[]{ a, b } );
			        		}
			        	}
			        }
    			}
    				
    				
    				
    			if( arbitraryPairs.size() > nonSubPairs.size() )
    				nonSubPairs = arbitraryPairs;
    			
    			if( nonSubPairs.size() > maxSubstituentRingSize ) {
    				ringPairs.addAll(nonSubPairs);
    			} else {
    				ringPairs.addAll(substituentPairs);
    			}
    				//ringPairs.addAll(nonSubPairs);
    			
    			
    			if( verbose )
    				System.out.println("-------------------------");
    		}
    	}
    		
    	
    	
    	
    	// uniquification via string
    	HashMap<String, int[]> uniquePairs = new HashMap<String, int[]>();
    	
    	for( int[] sPair : ringPairs ) {
    		IBond sBond = source.getBond(sPair[0]);
    		IBond tBond = target.getBond(sPair[1]);
    		
    		if( verbose )
    			System.out.println( "ring pair: " + sPair[0] + ", " + sPair[1] + " : " +  sBond.getAtom(0).getSymbol() + sBond.getAtom(1).getSymbol() + " - " + tBond.getAtom(0).getSymbol() + tBond.getAtom(1).getSymbol() );
    	
    		uniquePairs.put( sPair[0] + "," + sPair[1] , sPair );
    	}
    	
    	if( verbose ) {
    		System.out.println( "number of ring pairs - " + ringPairs.size() );
    		System.out.println( "number of unique ring pairs - " + uniquePairs.size() );
    	}
    	
    	for( int[] pair : uniquePairs.values() ) {
    		graphNodes.add(pair);
    	}
    	
    	//graphNodes.addAll( ringPairs );
    }
    
    
    
    private int ringPathDistance( int[] path, int start, int end, boolean reverse ) {
    	int distance = 0;
    	int pos = start;
    	
    	while( pos != end ) {
    		if( reverse ) {
    			--pos;
    			
    			if( pos < 0 )
    				pos = path.length - 1;
    		} else {
    			++pos;
    			
    			if( pos >= path.length )
    				pos = 0;
    		}
    		
    		++distance;
    	}
    	
    	return distance;
    }
    
    
    /**
     * Calculate this based on path-distances
     * 
     * in a benzene ring for example, para-substituted bonds should have a topological distance of 3 bonds from each other
     * 8-membered ring - 4 bonds, etc
     * 
     * return false if ring has an odd-numbered size
     * 
     * This also assumes that being topological distance-based, that no bridges will occur
     * 
     * @param ringPath
     * @param subInfo
     * @return
     */
    private List<int[]> getParaSubs( int[][] pathDistances, Map<Integer, int[]> subInfo, int ringSize ) {
    
    	
    	
    	if( ringSize % 2 != 0)
    		return null;
		
    	List<int[]> paraSubs = new ArrayList<int[]>( ringSize );
    	
		// int sStart = 0, tStart = 0;
		for( Integer s1 : subInfo.keySet() ) {
			for( Integer s2 : subInfo.keySet() ) {
				if( pathDistances[s1][s2] == (ringSize / 2) + 1 ) {
					paraSubs.add( new int[]{ s1, s2 } );
				}
			}
		}
		
		if( paraSubs.size() == 0 )
			return null;
			
		return paraSubs;
    }
    
    
    
    private List<List<int[]>> generateRingPaths( IAtomContainer molecule, List<List<Integer>> ringLists ) {
    	
    	//int[][] ringPaths = new int[ ringLists.size() ][];
    	List<List<int[]>> allRingPaths = new ArrayList<List<int[]>>();
    	int[][] adjList = GraphUtil.toAdjList(molecule);
    	
    	CycleFinder cf = Cycles.mcb();
    	
		
		for( int ringInd = 0 ; ringInd < ringLists.size(); ringInd++ ) {
			HashSet<Integer> vertices = new HashSet<Integer>();
			List<Integer> ring = ringLists.get( ringInd );
			List<int[]> ringPaths = new ArrayList<int[]>();
			
			for( Integer bInd : ring ) {
				vertices.add( molecule.getAtomNumber( molecule.getBond(bInd).getAtom(0) ) );  
				vertices.add( molecule.getAtomNumber( molecule.getBond(bInd).getAtom(1) ) );  
			}
			
			int[] vs = new int[ vertices.size() ];
			int v = 0;
			for( Integer vInd : vertices ) {
				vs[ v++ ] = vInd;
			}
			
			int[][] subAdjList = GraphUtil.subgraph(adjList, vs);
			try {
				  int[][] paths = cf.find(molecule, subAdjList, 22).paths();
				  
				  if( verbose )
					  System.out.println( "paths - " + paths.length );
				  
				  for( int pInd = 0; pInd < paths.length; pInd++ ) {
				  
					  int[] path = paths[ pInd ];
					  int[] bondPathIndices = new int[ path.length - 1 ];
					  
					  IBond[] bondPath = new IBond[path.length - 1];
		    	        // map each bond back to the parent
		    	      for (int i = 0, last = path.length - 1; i < last; i++) {
		    	            bondPath[i] = molecule.getBond( 
		    	            		molecule.getAtom( vs[ path[i  ] ] ), 
		    	            		molecule.getAtom( vs[ path[i+1] ] ) 
		    	            );
		    	            Integer bondIndex = molecule.getBondNumber(bondPath[i]);
 
		    	            
		    	            bondPathIndices[i] = bondIndex;
		    	            
		    	            if( verbose )
		    	            	System.out.println( pInd + ": " + paths[pInd][i] + "," + bondIndex + " " );
		    	      }
		    	        
		    	        ringPaths.add( bondPathIndices );
				  }
				  
				
				  
			} catch (Intractable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			allRingPaths.add(ringPaths);
		}
    	
		return allRingPaths;
    }
    
    
    /**
     * Traverses the bonds along the path of a ring.  
     * 
     * If said bond and it's next bond have a common vertex with a substituent (degree > 2) that is non-aromatic, store the indices
     * of the 2 bonds, and map to the substituent bond's index
     * 
     * @param bondIndices
     * @param molecule
     * @return
     */
    
    private Map<Integer, int[]> getRingSubstituents( List<int[]> ringPaths, IAtomContainer mol ) { 
    	
    	Map<Integer, int[]> subToRingBonds = new HashMap<Integer, int[]>();
    	
    	
    	for( int[] ringPath : ringPaths ) {
			for( int b1 = 0, b2 = 1; b1 < ringPath.length; b1++, b2++ ) {
	
				if( b2 == ringPath.length )
					b2 = 0;
				
				IBond bond1 = mol.getBond( ringPath[b1] );
				IBond bond2 = mol.getBond( ringPath[b2] );
				
				IAtom commonVertex1 = null;
		    	if( bond1.contains( bond2.getAtom(0) ) ) {
		    		commonVertex1 = bond2.getAtom(0);
		    	} else if( bond1.contains( bond2.getAtom(1) ) ) {
		    		commonVertex1 = bond2.getAtom(1);
		    	}
		    	
		    	
		    	List<IBond> cvNbs = mol.getConnectedBondsList(commonVertex1);  // slow, adjacency list would be better
		    	
		    	for( IBond nb : cvNbs ) {
		    		if( nb != bond1 && nb != bond2 ) {
		    			int subIndex = mol.getBondNumber(nb);
		    			
		    			if( ! ConvenienceTools.isAromatic( nb ) )  // not in an aromatic ring
		    				subToRingBonds.put(subIndex, new int[]{ ringPath[b1], ringPath[b2] } );
		    				
		    			break;
		    		}
		    	}
				
				
								
			}    		
		}
    	
    	return subToRingBonds;
    	
    }
    
    
    
    /**
     * Thanks to John May for this code, for finding rings in fused ring systems
     * 
     * @param molecule
     * @return
     */
    private ArrayList<Integer> findFusedBonds( IAtomContainer molecule, int[] fused ) {
    	
    	
    	ArrayList<Integer> bondRingCounts = new ArrayList<Integer>( molecule.getBondCount() );
    	for( int b = 0; b < molecule.getBondCount(); b++ ) {
    		bondRingCounts.add( 0 );
    	}
    	
    	EdgeToBondMap bondMap = EdgeToBondMap.withSpaceFor(molecule);
    	int[][] adjList = GraphUtil.toAdjList(molecule, bondMap);

    	RingSearch ringSearch = new RingSearch(molecule, adjList);

    	CycleFinder cf = Cycles.mcb();
    	int maxPathLen = molecule.getBondCount();
    	                       
    	// for each fused ring system
    	//for (int[] fused : ringSearch.fused()) {
    	    
    	    // take the subgraph and find the MCB
    	    int[][] subAdjList = GraphUtil.subgraph(adjList, fused);
    	    int[][] paths = null;
			try {
				paths = cf.find(molecule, subAdjList, maxPathLen).paths();
			} catch (Intractable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
    	                                         
    	    // for each ring in this system...
    	    for (int[] path : paths) {
    	        IBond[] bondPath = new IBond[path.length - 1];
    	        // map each bond back to the parent
    	        for (int i = 0, last = path.length - 1; i < last; i++) {
    	            bondPath[i] = bondMap.get(fused[path[i]], fused[path[i+1]]);
    	            Integer bondIndex = molecule.getBondNumber(bondPath[i]);
    	            bondRingCounts.set( bondIndex, bondRingCounts.get(bondIndex) + 1 );
    	            System.out.print( bondIndex + ", " +  bondRingCounts.get(bondIndex) + "| " );
    	        }
    	        
    	        if( verbose )
    	        	System.out.println(bondPath.length + " ");
    	    }
    	    
    	    if( verbose )
    	    	System.out.println();
    	//}
    	
    	ArrayList<Integer> fusedBonds = new ArrayList<Integer>( molecule.getBondCount() );
    	for( int b = 0; b < bondRingCounts.size(); b++ ) {
    		if( bondRingCounts.get(b) > 1 )
    			fusedBonds.add(b);
    	}
    	
    	return fusedBonds;
    }
    
    
    
    private ArrayList<int[]> matchRingsViaPaths( int[] sPath, int[] tPath, int sStart, int tStart, boolean reverse ) { 
    	
    	ArrayList<int[]> ringPairs = new ArrayList<int[]>( sPath.length );	
    	
    	for( int n = 0; n < sPath.length; n++ ) {
    		
    		IBond sBond = source.getBond( sPath[ sStart ] );
			IBond tBond = target.getBond( tPath[ tStart ] );
			
			if( nodesMatch(source, sBond, target, tBond, shouldMatchBonds, atomTypeProperty) ) {
				ringPairs.add( new int[]{ sPath[ sStart ], tPath[ tStart ] } );
				
				//currSubstituentRingSize.add( sIndex );
				//substituentsPerformed = true;
			}
			

				
				++sStart;
				
				if( sStart >= sPath.length )
					sStart = 0;
				
				if( reverse ) {
					--tStart;
					
					if( tStart < 0 )
						tStart = tPath.length - 1;
				} else {
					++tStart;
					
					if( tStart >= tPath.length )
						tStart = 0;
				}
				
				
			
			
    		
    	}
    	
    	return ringPairs;
    }
    
    
    private ArrayList<int[]> matchRingsViaRefPoint( List<Integer> sRing, List<Integer> tRing, int sRefPoint, int tRefPoint ) {
    	

		ArrayList<int[]> ringPairs = new ArrayList<int[]>( sRing.size() );	
    	
			/*
			 * enumerate bonds of relevant ring of mol1 & mol2
			 * map path distance from substituent, to the bond index
			 * Then create pairs from the maps, with the substituent
			 */
			HashMap<Integer, ArrayList<Integer>> sSubDistances = new HashMap<Integer, ArrayList<Integer>>();
			HashMap<Integer, ArrayList<Integer>> tSubDistances = new HashMap<Integer, ArrayList<Integer>>();
			
			for( Integer bIndex : sRing ) {
				
				int pDist = pathDistancesHsMol[ sRefPoint ][ bIndex ];
				
				if( ! sSubDistances.containsKey( pDist ) ) {
					sSubDistances.put( pDist, new ArrayList<Integer>(6) );
				}
					
				sSubDistances.get( pDist ).add( bIndex );
			}
			
			for( Integer bIndex : tRing ) {
				
				int pDist = pathDistancesQMol[ tRefPoint ][ bIndex ];
				
				if( ! tSubDistances.containsKey( pDist ) ) {
					tSubDistances.put( pDist, new ArrayList<Integer>(6) );
				}
					
				tSubDistances.get( pDist ).add( bIndex );
			}
			
			if( verbose ) {
			System.out.println( sSubDistances );
			System.out.println( tSubDistances );
			}
			
			// now to construct pairs
			if( sRing.size() == tRing.size() && sRing.size() <= 7 ) {  // simple rings, match via paths to preserve topological matching
				
				
				for( Integer pDist : sSubDistances.keySet() ) {
					
					List<Integer> sDists = sSubDistances.get(pDist);
					List<Integer> tDists = tSubDistances.get(pDist);
					
					for( int n = 0; n < sDists.size(); n++ ) {
						
						
							IBond sBond = source.getBond( sDists.get(n) );
							IBond tBond = target.getBond( tDists.get(n) );
							
							if( nodesMatch(source, sBond, target, tBond, shouldMatchBonds, atomTypeProperty) ) {

								ringPairs.add( new int[]{ sDists.get(n), tDists.get(n) } );
								
								
								//currSubstituentRingSize.add( sIndex );
								//substituentsPerformed = true;
							}
						
					}
				}
				
			} else {
			
				for( Integer pDist : sSubDistances.keySet() ) {
					for( Integer sIndex : sSubDistances.get(pDist) ) {
						
						if( ! tSubDistances.containsKey(pDist) )
							continue;
						
						for( Integer tIndex : tSubDistances.get(pDist) ) {
							IBond sBond = source.getBond(sIndex);
							IBond tBond = target.getBond(tIndex);
							
							if( nodesMatch(source, sBond, target, tBond, shouldMatchBonds, atomTypeProperty) ) {
								 
								ringPairs.add( new int[]{ sIndex, tIndex } );
								
								
								//currSubstituentRingSize.add( sIndex );
								//substituentsPerformed = true;
							}
						}
					}
				}
				
			}
			
			return ringPairs;
    	
    }
    
    
    /**
     * Induced subgraph heuristic - C-C bonds
     * 
     * @param molecule
     * @return a map for the bonds of a molecule, with a number identifying which induced subgraph it is in (or -1 if it isn't in one)
     */
    private ArrayList<Integer> findInducedSubgraphBonds( IAtomContainer molecule, inducedSubgraphType type ) {
    
    	
    	
    	ArrayList<Integer> desiredBonds = new ArrayList<Integer>( molecule.getBondCount() );
    	ArrayList<Integer> subgraphToOrig = new ArrayList<Integer>( molecule.getBondCount() );
    	ArrayList<Integer> bondComponentMap = new ArrayList<Integer>( molecule.getBondCount() );
    	 
    	if( type == inducedSubgraphType.RIGID ) {
    		for( int b = 0; b < molecule.getBondCount(); b++ ) {
	    		IBond bond = molecule.getBond(b);
	    		if( ! isRotatableBond(molecule, bond) )
	    			desiredBonds.add( b );
	    	}
    	} else if( type == inducedSubgraphType.NONRING ) {
    		for( int b = 0; b < molecule.getBondCount(); b++ ) {
	    		IBond bond = molecule.getBond(b);
	    		if( /*isRotatableBond(molecule, bond) */ ! ConvenienceTools.isRingBond(bond) )
	    			desiredBonds.add( b );
	    	}
    	} else {
	    	// note all C-C bonds
	    	for( int b = 0; b < molecule.getBondCount(); b++ ) {
	    		IBond bond = molecule.getBond(b);
	    		if( isMatchBond() ) {
		    		if( bond.getAtom(0).getAtomicNumber() == 6 && bond.getAtom(1).getAtomicNumber() == 6 && 
		    				bond.getOrder() == IBond.Order.SINGLE && ! ConvenienceTools.isAromatic(bond) )
		    			desiredBonds.add( b );
	    		} else {
	    			if( bond.getAtom(0).getAtomicNumber() == 6 && bond.getAtom(1).getAtomicNumber() == 6 )
		    			desiredBonds.add( b );
	    		}
	    	}
    	}
    	//int[] singleCarbonBondsArray = new int[ singleCarbonBonds.size() ];
    	//for( int n = 0; n < singleCarbonBonds.size(); n++ ) { singleCarbonBondsArray[n] = singleCarbonBonds.get(n); }
    	
    	if( verbose )
    		System.out.println( type + " requested carbon bonds - " + desiredBonds );
    	
    	
    	IAtomContainer subgraph = new org.openscience.cdk.AtomContainer();
    	for( int b = 0; b < desiredBonds.size(); b++ ) {
    		Integer subgraphInd = desiredBonds.get(b);
    		subgraphToOrig.add( subgraphInd );
    		
    		IBond bond = molecule.getBond(subgraphInd);
    		subgraph.addBond( bond );
    		
    		if( ! subgraph.contains( bond.getAtom(0) ))
    			subgraph.addAtom( bond.getAtom(0) );
    		
    		if( ! subgraph.contains( bond.getAtom(1) ))
    			subgraph.addAtom( bond.getAtom(1) );
    	}
    	
    	// remove subgraph bonds with no neighbours
    /*	for( int bInd = 0, bs = subgraph.getBondCount(); bInd < bs; bInd++ ) {
    		IBond b = subgraph.getBond(bInd);
    		if( b == null )
    			continue;
    		
    		if( subgraph.getConnectedAtomsCount( b.getAtom(0) ) <= 1 && subgraph.getConnectedAtomsCount( b.getAtom(1) ) <= 1 ) {
    			subgraph.removeAtom(b.getAtom(0));
    			subgraph.removeAtom(b.getAtom(1));
    			subgraph.removeBond(b);
    		}
    	}*/
    	
    	if( verbose ) {
	    	try {
				System.out.println( type + " subgraph - " + sg.create( subgraph ) );
			} catch (CDKException e1) {
				e1.printStackTrace();
			}
    	}
    	
    	// this stage is used to identify the connected components in the created (currently disconnected) subgraph
    	// namely - a map showing which bond is in which subgraph?
    	int[][] g  = GraphUtil.toAdjList( subgraph );
    	ConnectedComponents cc = new ConnectedComponents(g);
    	int[]  components = cc.components();
    	HashMap<Integer, Integer> atomToComponent = new HashMap<Integer, Integer>();
    	
    	for (int v = 0; v < g.length; v++) {
    		atomToComponent.put( v, components[v] );
    	}
    	
    	
    	// autoboxing
    	for( int b = 0; b < molecule.getBondCount(); b++ ) {
    		bondComponentMap.add( -1 );  // -1 means that a given bond is not in an induced subgraph
    	}
    	
    	// translates a bond ID to a particular identified induced subgraph
    	for( int b = 0; b < subgraph.getBondCount(); b++ ) {
    		int subgraphID = atomToComponent.get( 
    				subgraph.getAtomNumber( subgraph.getBond(b).getAtom(0) ) 
    		);
    		
    		bondComponentMap.set( subgraphToOrig.get(b), subgraphID );
    	}
    	
    	

    	
    	
    	
    	return bondComponentMap;
    }
    
    
    
    
    /**
	 *  The method determines if a bond is rotatable or not
	 *  
	 *  Based on the methodology used in SmartRotatableBondsCountDescriptor (CDK KNIME implementation)
	 *
	 *@param  ac                AtomContainer
	 *@return                   number of rotatable bonds
	 */
    public boolean isRotatableBond(IAtomContainer ac, IBond bond) {


			IAtom atom0 = bond.getAtom(0);
			IAtom atom1 = bond.getAtom(1);
            if (atom0.getSymbol().equals("H") || atom1.getSymbol().equals("H")) 
            	return false;
            
            if ( ! ConvenienceTools.isRingBond(bond) && bond.getOrder() == CDKConstants.BONDORDER_SINGLE) {
            	
            	if (isAmide(atom0, atom1, ac) || isAmide(atom1, atom0, ac)) {
            		return false;
            	}
            	
            	// this thing reduces accuracy
            	// TODO we count one-bond ring bridges as being non-rotatable (though in reality they can be)
            	// only however, if they're terminal atoms next to an atom
            	if( atom0.getFlag(CDKConstants.ISINRING) || atom1.getFlag(CDKConstants.ISINRING) ) {
            		// either is a terminal atom
            		if( ac.getConnectedAtomsCount(atom0) == 1 || ac.getConnectedAtomsCount(atom1) == 1 )
            			return false;
            		
            		// single bond link between rings where at least one of the rings is aromatic
            		if( (atom0.getFlag(CDKConstants.ISINRING) & atom1.getFlag(CDKConstants.ISINRING)) &&
            			(atom0.getFlag(CDKConstants.ISAROMATIC) | atom1.getFlag(CDKConstants.ISAROMATIC))	)
            			return false;
            	}
            	
            	/*// anything immediate aromatic ring substitution we also fix
            	if( ConvenienceTools.isAromatic(atom0) | ConvenienceTools.isAromatic(atom1) )
            		return false;*/
            		
            	
            	
            	// things across triple bonds are not rotatable
				if ((BondManipulator.isLowerOrder(ac.getMaximumBondOrder(atom0), IBond.Order.TRIPLE)) && 
					(BondManipulator.isLowerOrder(ac.getMaximumBondOrder(atom1), IBond.Order.TRIPLE))) {
					
                     return true;
				}
			}
		

            return false;
	}
	
	private boolean isAmide(IAtom atom0, IAtom atom1, IAtomContainer ac) {
		
		// ignores
    	// (a) 1,3-H tautomers where HO-C=N
    	// (b) charge states where (-)O-C=N
    	// (c) anything else that does not work but should
    	if (atom0.getSymbol().equals("C") && atom1.getSymbol().equals("N")) {
    		for (IAtom neighbor : ac.getConnectedAtomsList(atom0)) {
    			if (neighbor.getSymbol().equals("O") && ac.getBond(atom0, neighbor).getOrder() == CDKConstants.BONDORDER_DOUBLE) {
    				return true;
    			}
    		}
    	}
    	return false;
	}
    
	
	
    /**
     * specify a row in the modular products adjacency list/matrix.  Returns the intersection of the supplied set, with said row.
     * 
     * @param index
     * @param set
     * @return
     */
    public Collection<Integer> intersectNodeNeighbours( int index, Collection<Integer> set ) {
    	
    	if( set.isEmpty() )
    		return new ArrayList<Integer>(0);
    	
    	Collection<Integer> neighbours = edges.get(index);
    	//Collection<Integer> oSet = new HashSet<Integer>( set );
    	Collection<Integer> oSet = null;
    	oSet = ConvenienceTools.createCollection(set);
    	
    	/*try {
			Constructor constructor = set.getClass().getConstructor(new Class[]{ Collection.class });

			oSet = (Collection<Integer>) constructor.newInstance( set );
		} catch (InstantiationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
    	
    	
    	//System.out.println( "class of collection is " + oSet.getClass().toString() );
    	oSet.retainAll(neighbours);
    	
    	return oSet;
    }
    
  
    public void clear() {
    	graphNodes = null;
    	edges = null;
    }
    
    public List<int[]> getNodes() {
    	return graphNodes;
    }
    
    public List<Collection<Integer>> getAdjacencyList() {
    	return edges;
    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isMatchBond() {
        return shouldMatchBonds;
    }

    /**
     * @param shouldMatchBonds the shouldMatchBonds to set
     */
    public void setMatchBond(boolean shouldMatchBonds) {
        this.shouldMatchBonds = shouldMatchBonds;
    }
}
