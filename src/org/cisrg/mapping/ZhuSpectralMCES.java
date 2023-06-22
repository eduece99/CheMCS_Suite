package org.cisrg.mapping;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.zip.DataFormatException;
import java.util.Map.Entry;

//import org.apache.xmlbeans.impl.piccolo.io.FileFormatException;
import org.cisrg.utilities.HungarianAlgorithm;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.matrix.AdjacencyMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.AtomContainer;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.doublealgo.Transform;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.jet.math.Functions;
//import chemaxon.calculations.topology.SimilarityChecker;

public class ZhuSpectralMCES extends MCSMethods {

	
	
	
	// descending order of bond score
		protected static Comparator<double[]> similarityComparator = new Comparator<double[]>() {
			 public int compare(double[] o1, double[] o2) {
				 	
		            if( o1[2] == o2[2] ) {
		            	return 0;
		            } 
		            
		            return o1[2] > o2[2] ? -1 : 1;
		           
			 } 
		};
		
		
		// pair object
		protected class SimPair implements Comparable<SimPair> {
			
			protected SimPair( int[] p, double s ) {
				pair = p;
				similarity = s;
			}
			
			@Override
			public int hashCode() {
				return pair[0] * pair[1];
			}
			
			@Override
			public int compareTo(SimPair o) {
				 return Double.compare( similarity, o.getSimilarity() ) * -1;
			}
			
			@Override
			public boolean equals( Object o ) {
				
				if( o instanceof SimPair ) {
					int[] oPair = ((SimPair) o).getPair();
					
					if( pair[0] == oPair[0] && pair[1] == oPair[1] )
						return true;
				}
				
				return false;
			}
			
			@Override
			public String toString() {
				return "[" + pair[0] + "," + pair[1] + "] " + similarity;
			}
			
			
			public int[] getPair() {
				return pair;
			}
			
			public double getSimilarity() {
				return similarity;
			}
			
			
			private int[] pair;
			private double similarity;
			
		}
	
	
	public static IAtomContainer createSubgraph( IAtomContainer sourceMol, Collection<Integer> atoms, GraphUtil.EdgeToBondMap etbMap ) {
		
		IAtomContainer nMol = new AtomContainer();
		for( Integer i : atoms ) {
			nMol.addAtom( sourceMol.getAtom(i) );
			for( Integer j : atoms ) {
				
				IBond b = etbMap.get(i, j);
				
				if( b != null && ! nMol.contains(b) ) {
					nMol.addBond( b );
				}
			}
		}
		
		return nMol;
	}
	
	
	
	
	protected DoubleMatrix2D globalSimilarity( IAtomContainer g1, IAtomContainer g2 ) {
		
		// adjacency matrices
		DenseDoubleMatrix2D g1AdjMat = new DenseDoubleMatrix2D( ConvenienceTools.int2DtoDouble2D( AdjacencyMatrix.getMatrix(g1) ) );
		DenseDoubleMatrix2D g2AdjMat = new DenseDoubleMatrix2D( ConvenienceTools.int2DtoDouble2D( AdjacencyMatrix.getMatrix(g2) ) );
		/*DenseDoubleMatrix2D g1AdjMat = new DenseDoubleMatrix2D( ConvenienceTools.int2DtoDouble2D( ConvenienceTools.bondAdjacencyMatrix(g1) ) );
		DenseDoubleMatrix2D g2AdjMat = new DenseDoubleMatrix2D( ConvenienceTools.int2DtoDouble2D( ConvenienceTools.bondAdjacencyMatrix(g2) ) );*/
		
		
		
		// get diagonal matrices of degrees per graph
				double[] degreesN1 = new double[ adj1.length ];
				double[] degreesN2 = new double[ adj2.length ];
				
				for( int n1 = 0; n1 < adj1.length; n1++ ) {
					degreesN1[n1] = adj1[n1].length ;
				}
				
				for( int n2 = 0; n2 < adj2.length; n2++ ) {
					degreesN2[n2] = adj2[n2].length;
				}
		
				return globalSimilarity( g1AdjMat, g2AdjMat, degreesN1, degreesN2 );
	}
	
	
	
	// input adjacency matrices
	public static DoubleMatrix2D globalSimilarity( DenseDoubleMatrix2D g1AdjMat, DenseDoubleMatrix2D g2AdjMat, double[] degreesN1, double[] degreesN2 ) {
		
		DoubleMatrix2D diagDegrees1 = DoubleFactory2D.dense.diagonal( new DenseDoubleMatrix1D(degreesN1) );
		DoubleMatrix2D diagDegrees2 = DoubleFactory2D.dense.diagonal( new DenseDoubleMatrix1D(degreesN2) );
		
		//System.out.println(diagDegrees1); )
		//System.out.println(g1AdjMat);
		//System.out.println( diagDegrees1.assign( g1AdjMat, Functions.minus ) );
		
		// eigenvector obtaining
		EigenvalueDecomposition evd1 = new EigenvalueDecomposition( diagDegrees1.assign( g1AdjMat, Functions.minus ) );
		EigenvalueDecomposition evd2 = new EigenvalueDecomposition( diagDegrees2.assign( g2AdjMat, Functions.minus ) );
		
		
		/* As the eigenvalues are already sorted for us (as are eigenvectors by eigenvalue), we simply apply this "shortcut" to ensure we're
		 * using the same number of columns per matrix, where the excluded columns (from the larger matrix) are those with the smallest
		 * eigenvalues 
		 */
		int atomCount = Math.min( degreesN1.length, degreesN2.length );
		int start1 = degreesN1.length - degreesN2.length;
		int start2 = degreesN2.length - degreesN1.length;
		
		// Ensure only the bigger matrix is adjusted (column-wise)
		if( start1 < 0 ) { start1 = 0; }
		if( start2 < 0 ) { start2 = 0; }
		
		
		//System.out.println( evd2.getV().viewPart(0, start2, evd2.getV().rows(), atomCount) );
		
		DoubleMatrix2D results = algebra.mult( 
				Transform.abs( evd1.getV().viewPart(0, start1, evd1.getV().rows(), atomCount) ), 
				Transform.abs( algebra.transpose( evd2.getV().viewPart(0, start2, evd2.getV().rows(), atomCount) ) )
		);
		
		//results.assign( Functions.mult(1.101332) );
		
		return results;
	}
	
	
	
	
	/**
	 * Labelled local similarity matrix
	 * 
	 * XXX changes:
	 * 
	 * Zhu et al's work described the creation of n-order subgraphs and matching nodes via a bipartite graph in the subgraphs.  
	 * The bipartite weight would be equivalent to the intersection of the labels in each subgraph, of the nodes being matched (thus
	 * nodes with different labels would have no weight).
	 * 
	 * I'm not doing this due to SMARTS (or any other degenerate) matching.  Instead I'll have to use 2 for loops to "count" the number
	 * of compatible atom matches per subgraph, and construct the bipartite weightings that way.
	 * 
	 * This method uses the neighbourhood graphs (that is, the vertex at the centre from which the neighbourhood was 
	 * contained, is excluded for computing the weights.
	 * 
	 * Note also that due to suggestions from the authors of the paper, degrees here (or equivalent) are obtained from the vertex in the 
	 * source graph, as opposed to the subgraph (i.e the degree values may be larger).  This was done by Zhu et al presumably to save 
	 * computation time, contrary to what the methodology in the original paper presented. 
	 * 
	 * @param g1
	 * @param g2
	 * @param nOrder
	 * @param labels
	 * @return
	 */
	protected DoubleMatrix2D localSimilarityAtoms( IAtomContainer g1, IAtomContainer g2, int nOrder ) {
		
		double[][] simMat = new double[ adj1.length ][ adj2.length ];  // default values of 0.0
		
		List<Set<Integer>> g1SgList = new ArrayList<Set<Integer>>( adj1.length );
		List<Set<Integer>> g2SgList = new ArrayList<Set<Integer>>( adj2.length );
		
		

		
		/* pre-calculation of neighbourhood information to save time
		 * 
		 *  XXX  Note that the central vertex is included in this "neighbourhood" information
		 *  
		*/
		for( int u = 0; u < adj1.length ; u++ ) {
			Set<Integer> g1SgElemsWithU = getNeighbourhood(adj1, u, nOrder);
			g1SgElemsWithU.add( u );
			g1SgList.add( g1SgElemsWithU );
			
			/*double uMatches = countNeighbourhoodAtomMatches(g1, adj1, u);
			uMatchesList[u] = uMatches;*/
		}
		
		for( int v = 0; v < adj2.length ; v++ ) {
			Set<Integer> g2SgElemsWithV = getNeighbourhood(adj2, v, nOrder);
			g2SgElemsWithV.add( v );
			g2SgList.add( g2SgElemsWithV );
			
			/*double vMatches = countNeighbourhoodAtomMatches(g2, adj2, v);
			vMatchesList[v] = vMatches;*/
		}
		
		
		
		for( int u = 0; u < adj1.length ; u++ ) {
			
			// "subgraph" of graph 1
			Set<Integer> g1SgElemsWithU = g1SgList.get(u);
			IAtom uElem = g1.getAtom( u );
			
			
			// first, count the number of matches between central vertex and the neighbourhood vertices
			// double uMatches = uMatchesList[u];
			
			for( int v = 0; v < adj2.length; v++ ) {
				
				// label checking - automatic value of 0 if labels don't match
				IAtom vElem = g2.getAtom( v );
				if( ! atomsMatch(g1, uElem, g2, vElem) ) {
					simMat[u][v] = avoidSimilarityValue;
					continue;
				}
				
				// "subgraph" for other graph
				Set<Integer> g2SgElemsWithV = g2SgList.get(v);
				
				//double vMatches = vMatchesList[v];

				//double mainMatches = Math.min(uMatches, vMatches);
				double mainMatches = countNeighbourhoodAtomMatches( g1, g2, adj1, adj2, u, v );
				

				double[][] localWeights = new double[ g1SgElemsWithU.size() ][ g2SgElemsWithV.size() ];

				int i = 0;
				for( Integer iIndex : g1SgElemsWithU ) {
					int j = 0;
					//double li = countNeighbourhoodAtomMatches(g1, adj1, iIndex);
					
					for( Integer jIndex : g2SgElemsWithV ) {

						
						//double lj = countNeighbourhoodAtomMatches(g2, adj2, jIndex);

						// negative sign due to the way the Hungarian algorithm is implemented (minimum cost)
						//localWeights[i][j] = -1 * Math.min( li, lj );
						localWeights[i][j] = -1 * countNeighbourhoodAtomMatches( g1, g2, adj1, adj2, iIndex, jIndex );
						
							

						++j;
					}
					++i;
				} 


				HungarianAlgorithm ha = new HungarianAlgorithm( localWeights );
				int[] optimalMatches = ha.execute();

				//if( verbose )	System.out.println( "optimal labels : " + Arrays.toString( optimalMatches ) );

				double weightSum = 0;
				for( int n = 0; n < optimalMatches.length; n++ ) {
					if( optimalMatches[n] >= 0 )
						weightSum += localWeights[n][ optimalMatches[n] ];
				}
				weightSum *= -1;  // restore positive value

				// XXX	The original authors divided this by 2.  This was removed (and still works)
				double neighbourhoodDistance = (mainMatches + weightSum) ;
				
				// XXX	degree values from graph's adjacency list, as opposed to the subgraph
				int sumU = adj1[u].length;
				int sumV = adj2[v].length;
				for( Integer vertex : g1SgElemsWithU ) { sumU += adj1[vertex].length; }
				for( Integer vertex : g2SgElemsWithV ) { sumV += adj2[vertex].length; }
				
				int minNeighbours = Math.min( g1SgElemsWithU.size() , g2SgElemsWithV.size() );
				
				// XXX	the "+1" in the original equation was removed
				double numerator = Math.pow( ( minNeighbours + neighbourhoodDistance ), 2 );
				double normaliser = ( g1SgElemsWithU.size() + sumU ) * ( g2SgElemsWithV.size() + sumV );
				
				simMat[u][v] = numerator / normaliser;
			}
		}
		
		return new DenseDoubleMatrix2D(simMat);
		
	}
	
	
	protected Set<Integer> getNeighbourhood( int[][] adjList, int vertex, int nOrder ) {
		
		Set<Integer> sgBonds = new HashSet<Integer>( adjList.length );
		sgBonds.add(vertex);  // add initial vertex to allow for finding neighbours
		for( int n = 0; n < nOrder; n++ ) {
			Set<Integer> nRadius = new HashSet<Integer>( adjList.length );
			
			for( Integer b : sgBonds ) {
				for( int b2 : adjList[b] ) {
					nRadius.add( b2 );
				}
			}
			
			sgBonds.addAll(nRadius);
		}
		
		
		sgBonds.remove(vertex);  // the neighbourhood should not include the central vertex
		
		return sgBonds;
	}
	
	
	/**
	 * Find adjacent nodes of the specificed node within the subgraph - count matches of the adjacents  
	 * 
	 * takes 2 indices (and their adjacency matrices) as the input
	 * for each L() of the 1st vertex being tested
	 *  - if a matching vertex exists in the other L(), increment matches
	 *  - remove vertex from consideration to prevent duplicates
	 * 
	 *  This' supposed to be the same as taking the intersection of the 2 multisets L(u) & L(v) in the paper
	 *  
	 *  
	 * @param g1
	 * @param g2
	 * @param neighbourhoodAtoms
	 * @param mainIndex
	 * @return
	 */
	
	protected double countNeighbourhoodAtomMatches( IAtomContainer g1, IAtomContainer g2, int[][] adjList1, int[][] adjList2, int ind1, int ind2 ) {
		 
				double matches = 0;
				
				Set<Integer> vSet2 = new HashSet<Integer>();
				for( int v2 : adjList2[ ind2 ] ) {
					vSet2.add(v2);
				}
				
				for( int v1 : adjList1[ ind1 ] ) {
					IAtom at1 = g1.getAtom(v1);
					
					for( int v2 : vSet2 ) {
						IAtom at2 = g2.getAtom(v2);
						
						// *one* match for atom in question is good enough
						if( atomsMatch(g1, at1, g2, at2 ) ) {
							matches++;
							vSet2.remove(v2);
							break;
						}
					}
				}
				
				return matches;
	}
	
	
	/**
	 * Finds "L(u)" - the number (originally the set) of compatible labels adjacent to the specified node
	 * 
	 * 
	 * 
	 * @param g
	 * @param adjList
	 * @param neighbourhoodAtoms
	 * @param mainIndex
	 * @return
	 */
	protected int countNeighbourhoodAtomMatches( IAtomContainer g, int[][] adjList, int mainIndex ) {
		// first, count the number of matches between central vertex and the neighbourhood vertices
				int matches = 0;
				IAtom mainElem = g.getAtom(mainIndex);
				for( int b : adjList[ mainIndex ] ) {
					IAtom oElem = g.getAtom(b);
					
					if( atomsMatch(g, mainElem, g, oElem ) )
						matches++;
				}
				
				return matches;
	}
	
	protected int countNeighbourhoodBondMatches( IAtomContainer g1, IAtomContainer g2, Set<Integer> neighbourhoodBonds, int mainIndex ) {
		// first, count the number of matches between central vertex and the neighbourhood vertices
				int matches = 0;
				IBond mainBond = g1.getBond(mainIndex);
				for( Integer b : neighbourhoodBonds ) {
					IBond bond = g2.getBond(b);
					
					if( GenerateCompatibilityGraphEdges.nodesMatch(g1, mainBond, g2, bond, matchBonds, null) )
						matches++;
				}
				
				return matches;
	}
	
	
	
	public static double largerAverageDegree( IAtomContainer g1, IAtomContainer g2 ) {
		
		return Math.max(
				(2 * g1.getBondCount()) / (double) g1.getAtomCount(),
				(2 * g2.getBondCount()) / (double) g2.getAtomCount()
		);
		
	}
	
	
	public static int maxDegree( int[][] adjList ) {
		int maxDegree = 0;
		for( int[] row : adjList ) {
			if( row.length > maxDegree )
				maxDegree = row.length;
		}
		
		return maxDegree;
	}
	
	
	/**
	 * Obtain anchors for expansion
	 * 
	 * XXX Changes from author's method:
	 * - I don't like the idea of specifying an arbitrary threshold as that's graph-dependent.  Here an arbitrary number of anchors
	 * is specified instead
	 * - the larger average degree is set to the smallest maximum degree (of the two graphs) if it exceeds said maximum degree
	 * 
	 * 
	 * 
	 * @param g1
	 * @param g2
	 * @param g1Used	Specify vertex indices in here (from g1) so that they are not considered as anchors
	 * @param g2Used
	 * @param simMat	
	 * @param numAnchors
	 * @param dConstraint	Whether to use the "larger average degree" constraint
	 * @return
	 */
	protected List< SimPair > findAnchors( IAtomContainer g1, IAtomContainer g2, Collection<Integer> g1Used, Collection<Integer> g2Used, DoubleMatrix2D simMat, int numAnchors, boolean dConstraint ) {
		
		
		//ArrayList<int[]> initialPairs = new ArrayList<int[]>();
		ArrayList< SimPair > anchors = new ArrayList< SimPair >( numAnchors );
		double lad = largerAverageDegree(g1, g2);
		
		if( g1Used == null )
			g1Used = new HashSet<Integer>();
		
		if( g2Used == null )
			g2Used = new HashSet<Integer>();
		
		
		/* XXX a CHANGE from original method:
		 * the larger average degree is unrealistic if one graph has far more bonds per atom than the other.
		 * So if this' the case it'll be changed to the maximum degree of that of the graph with the lower average degree
		 */
		int smallestMaxDegree = Math.max(maxDegree(adj1), maxDegree(adj2));
		
		if( smallestMaxDegree < lad )
			lad = smallestMaxDegree;

		
		//lad = lad - 0.0001;   
		
		List<SimPair> sortedPairs = new ArrayList<SimPair>( g1.getAtomCount() * g2.getAtomCount() );
		for( int a = 0; a < simMat.rows(); a++ ) {
			for( int b = 0; b < simMat.columns(); b++ ) {
				
				if( g1Used.contains(a) || g2Used.contains(b) )
					continue;
				
				double simValue = simMat.get(a,b);
				
				if( simValue < 0 )
					continue;
				
				if( dConstraint ) {
					if( Math.min( 
							adj1[ a ].length , 
							adj2[ b ].length
						) < lad ) {
						continue;
					}
				}
				
				
				sortedPairs.add( new SimPair( new int[]{a,b}, simValue ) );
			}
		}
		
		//System.out.println( sortedPairs );
		Collections.sort( sortedPairs );
		
		int n = 0;
		for( int i = 0; n < numAnchors && i < sortedPairs.size(); i++ ) {
			
			
			SimPair anchor = sortedPairs.get(i);
			
			
			if( g1Used.contains( anchor.getPair()[0] ) ||
				g2Used.contains( anchor.getPair()[1] )	
				)
				continue;

			anchors.add( anchor );
			g1Used.add( anchor.getPair()[0] );
			g2Used.add( anchor.getPair()[1] );
			
			++n;
		}
		
		if( verbose ) 	System.out.println("anchors: " + anchors) ;
		
		return anchors;
	}

	
	
	
	/**
	 * obtain optimal matching neighbour pairs, of the supplied vertex pair
	 * 
	 * @param g1
	 * @param g2
	 * @param simMat
	 * @param initialPair
	 * @return
	 */
	private Collection<int[]> obtainNeighbourPairs( IAtomContainer g1, IAtomContainer g2, final DoubleMatrix2D simMat, int[] initialPair ) {
		
		int[] g1Neighbours = adj1[ initialPair[0] ];
		int[] g2Neighbours = adj2[ initialPair[1] ];
		
		Set<Integer> nlUniqueRows = new HashSet<Integer>();
		Set<Integer> nlUniqueCols = new HashSet<Integer>();
		
		//Collection<int[]> nList = new LinkedList<int[]>();
		
		for( int i = 0; i < g1Neighbours.length; i++ ) {
			for( int j = 0; j < g2Neighbours.length; j++ ) {
				
				IAtom g1Atom = g1.getAtom( g1Neighbours[i] );
				IAtom g2Atom = g2.getAtom( g2Neighbours[j] );
				
				if( ! atomsMatch(g1, g1Atom, g2, g2Atom) )
					continue;
				
				/*double[] entry = new double[]{ 
						g1Neighbours[i], 
						g2Neighbours[j],
						simMat.get( g1Neighbours[i] , g2Neighbours[j] )
				};*/
				
				
				nlUniqueRows.add( g1Neighbours[i] );
				nlUniqueCols.add( g2Neighbours[j] );
			}
		}
		
		
		return similarityMatchList(simMat, nlUniqueRows, nlUniqueCols);
	}

	/**
	 * expand matching from anchors, neighbourhood-wise
	 * 
	 * Note that the return type is now int[] - this is for quick look-up in the refinement stages
	 * 
	 * @param g1
	 * @param g2
	 * @param simMat
	 * @param anchors
	 * @return
	 */
	private List<Integer> expandAnchors( IAtomContainer g1, IAtomContainer g2, final DoubleMatrix2D simMat, List< SimPair > anchors ) {
		
		TreeSet<SimPair> nList = new TreeSet<SimPair>(   );  // Q in the paper
		HashSet<Integer> g1Used = new HashSet<Integer>();
		HashSet<Integer> g2Used = new HashSet<Integer>();
		int upperBound = Math.min( g1.getAtomCount() , g2.getAtomCount() );
		
		// obtain Cartesian products of anchor neighbours
		for( int a = 0; a < anchors.size(); a++ ) {
			
			int[] anchorPair = anchors.get( a ).getPair();
			g1Used.add( anchorPair[0] );
			g2Used.add( anchorPair[1] );
			
			Collection<int[]> nPairs = obtainNeighbourPairs(g1, g2, simMat, anchorPair );
			
			for( int[] pair : nPairs ) {
				SimPair entry = new SimPair( pair, simMat.get( pair[0] , pair[1] ) ) ;
				
				nList.add( entry );
			}
		}
		
		
		
	/*	if( verbose )	System.out.println( "nList1: " + nList );
		
		if( verbose )	{
			System.out.println( "initial anchors: " + anchors );
			System.out.println( "nList1: " + nList );
		}
			

		if( verbose )	{
			System.out.println( );
			//Collections.sort( nList, rowComparator );
			for( ArrayList<Integer> ar : nList ) {
				System.out.print( ar + " " + simMat.get(ar.get(0), ar.get(1)) + " | " );
			}
			
			
			System.out.println();
		}*/
		
		//boolean test = true;
		
		
		while( ! nList.isEmpty() ) {
			
			SimPair spair = nList.pollFirst();  // pair with highest similarity is removed from queue
			int[] iPair = spair.getPair();
			
			if( ! g1Used.contains(iPair[0]) && ! g2Used.contains(iPair[1]) ) {

				anchors.add( spair );
				g1Used.add( iPair[0] );
				g2Used.add( iPair[1] );
				
				Collection<int[]> nPairs = obtainNeighbourPairs(g1, g2, simMat, iPair );
				
				for( int[] nPair : nPairs ) {
					SimPair entry = new SimPair( nPair, simMat.get( nPair[0] , nPair[1] ) ) ;
					
					nList.add( entry );
				}
			}
			
			/* XXX	Extra stage NOT in the paper:
			* Any nodes in either graph that've not been match, go through the anchor selection process
			* keep on finding anchors until all nodes in one of the graphs have been matched
			*/
			if( nList.isEmpty() && 
				anchors.size() < upperBound	
			) {
				nList.addAll( findAnchors(g1, g2, new HashSet<Integer>(g1Used), new HashSet<Integer>(g2Used), simMat, 1, false) );
			}
		}
		
		
		
		List<Integer> matches = new ArrayList<Integer>( g1.getAtomCount() );
		for( int n = 0; n < g1.getAtomCount(); n++ ) {
			matches.add(-1);  // autoboxing
		}
		
		for( SimPair sp : anchors ) {
			int[] pair = sp.getPair();
			matches.set( pair[0], pair[1] );
			//matches[ pair[0] ] = pair[1];
		}
		
		return matches;
	}
	
	
	
	
	
	/**
	 * XXX Approximate minimum vertex cover.  Isolation Algorithm by Ugurlu 2012 - "New heuristic algorithm for unweighted minimum vertex cover"
	 * 
	 * @param workBase  A starting adjacency list in the form of a Map (in which the algorithm performs removals on)
	 * @return
	 */
		public static List<Integer> findCoverIsolationAlgorithm( Map<Integer, Collection<Integer>> workBase, int[][] adjList ) {
	        // C <-- {}
			
			
			Set<Integer> cover = new HashSet<Integer>( workBase.size() );
			
			
			// removal of vertices with 0 degree.  Need a translation matrix to re-translate the cover back to the original ids
			List<Integer> removes = new ArrayList<>();
			for( Integer k : workBase.keySet() ) {
				if( workBase.get(k).size() == 0 )
					removes.add(k);
			}
			
			for( Integer r : removes ) {
				workBase.remove(r);
			}
			
			int origSize = workBase.size();
			
			int minDegree = origSize;
			int minDegreeAt = 0;
			
			/*List<Integer> minDegreeVertices = new ArrayList<Integer>( origSize );
			
			for( Integer at : workBase.keySet() ) {
				int degree = workBase.get(at).size();
        		if( degree < minDegree ) {
        			minDegreeVertices.clear();
        			minDegree = degree;
        			minDegreeVertices.add( at );
        		} else if( degree == minDegree ) {
        			minDegreeVertices.add( at );
        		}
        	}

			
			
			//int minIndex = (int) Math.floor( Math.random() * minDegreeVertices.size() );
			minDegreeAt = minDegreeVertices.get( 0 );*/
			
			//if( verbose )	System.out.println( "min vertices " + minDegreeVertices + " min vertex - " + minDegreeAt );
			
	        // while G' != {}
	        while ( ! workBase.isEmpty()  ) { 
	            // v <-- vertex with maximum degree in G'
	        	//int minDegreeAt = -1;
	        	minDegree = origSize;
	        	
	        	
	        	for( Integer at : workBase.keySet() ) {
	        		if( workBase.get(at).size() < minDegree && workBase.get(at).size() > 0 ) {
	        			minDegree = workBase.get(at).size();
	        			minDegreeAt = at;
	        		}
	        	}
	        	
	            
	            
	            // C <-- C U {v}
	            cover.add( minDegreeAt );
	            
	            // add all neighbours of the new vertex to the solution
	            Collection<Integer> neighbours = workBase.get( minDegreeAt );
	            cover.addAll( neighbours );
	            
	            
	            
	            // get neighbours of neighbours, store them
	            Map<Integer, Collection<Integer>> nNeighbours = new HashMap<Integer, Collection<Integer>>( );
	            for( Integer v : neighbours ) {
	            	if( workBase.containsKey( v ) ) {
	            		Collection<Integer> row = new ArrayList<Integer>( workBase.get( v ) );
	            		nNeighbours.put(v, row);
	            	}
            	}
	            
	            if(verbose)	System.out.println( minDegreeAt + " " + minDegree + " " + workBase.size() + " " + neighbours + " " + nNeighbours );

	            // remove from G' every edge incident on the neighbours (as well as original v)
	            Set<Integer> keys = new HashSet<Integer>( workBase.keySet() );
	            for( Integer rowKey : keys ) {
	            	Collection<Integer> row = workBase.get(rowKey);
	            	
	            	for( Integer v : nNeighbours.keySet() ) {
	            		Collection<Integer> neighbourRow = nNeighbours.get( v );
	            		
	            		row.removeAll( neighbourRow );  // removal of neighbours
	            	}
	            		
	            	row.removeAll( neighbours );
	            	
	            	if( row.isEmpty() ) {
	            		workBase.remove( rowKey );
	            	}
	            }
	            
	            for( Integer v : nNeighbours.keySet() ) {
            		Collection<Integer> neighbourRow = nNeighbours.get( v );
            		
            		 // removal of neighbours (as keys)
            		for( Integer nV : neighbourRow ) {
            			workBase.remove(nV);
            		}
            	}
	           
	            workBase.remove( minDegreeAt );
	            
	            
	        }
	        
	        

	        // removal of redundant vertices
	        Set<Integer> keys = new HashSet<Integer>( cover );
	        for( Integer v : keys ) {
	        	boolean redundant = true;
	        	
	        	for( int vNeighbour : adjList[v] ) {
	        		if( ! cover.contains( vNeighbour ) ) {
	        			redundant = false;
	        			break;
	        		}
	        	}
	        	
	        	if( redundant )
	        		cover.remove( v );
	        }

	        return new ArrayList<Integer>( cover );
	    }
		
	
	// XXX	Using greedy minimal cover rather than random minimal cover
	public static List<Integer> findGreedyCover( Map<Integer, Collection<Integer>> workBase ) {
        // C <-- {}
		
		
		ArrayList<Integer> cover = new ArrayList<Integer>( workBase.size() );
		

		
        // while G' != {}
        while ( ! workBase.isEmpty()  ) { 
            // v <-- vertex with maximum degree in G'
        	int maxDegreeAt = -1;
        	int maxDegree = -1;
        	
        	for( Integer at : workBase.keySet() ) {
        		if( workBase.get(at).size() > maxDegree ) {
        			maxDegree = workBase.get(at).size();
        			maxDegreeAt = at;
        		}
        	}
        	
            //System.out.println( maxDegreeAt.getProperty(idProperty) + " " + punchingBag.getBondCount() );
            
            // C <-- C U {v}
            cover.add( maxDegreeAt );

            // remove from G' every edge incident on v, and v itself
            //
            Collection<Integer> neighbours = workBase.get( maxDegreeAt );
            Set<Integer> keys = new HashSet<Integer>( workBase.keySet() );
            for( Integer rowKey : keys ) {
            	Collection<Integer> row = workBase.get(rowKey);
            	row.removeAll( neighbours );
            	row.remove( maxDegreeAt );
            	
            	if( row.isEmpty() ) {
            		workBase.remove( rowKey );
            	}
            }
            
            
           
            workBase.remove( maxDegreeAt );
            
            //punchingBag.removeAtom(maxDegreeAt);
            
        }
        
        //if( verbose )	System.out.println( "cover verification = " + isVertexCover(g, cover) );

        return cover;
    }
	
	
	
	/** XXX	Using greedy minimal cover rather than random minimal cover
	 * 
	 * Uses SRA algorithm by Balaji et al (2010) "An Effective Algorithm for Minimum Weighted Vertex Cover Problem"
	 * 
	 * @param workBase  Set of vertices to work with
	 * @return
	 */
		public static List<Integer> findCoverSRA( Map<Integer, Collection<Integer>> workBase ) {
	        // C <-- {}
			
			
			ArrayList<Integer> cover = new ArrayList<Integer>( workBase.size() );
			Map<Integer, Integer> supports = new HashMap<Integer, Integer>( workBase.size() );
			
			// calculate support values (sum of adjacent degrees to given vertex)
			for( Integer at : workBase.keySet() ) {
        		int support = 0;
        		for( int adj : workBase.get(at) ) {
        			support += workBase.get(adj).size();
        		}
        		supports.put( at, support );
        	}

			
	        // while G' != {}
	        while ( ! workBase.isEmpty()  ) { 
	            // v <-- vertex with maximum degree in G'
	        	int maxDegreeAt = -1;
	        	int maxDegree = -1;  // max support
	        	
	        	for( Integer at : workBase.keySet() ) {
	        		if( supports.get(at) > maxDegree ) {
	        			maxDegree = supports.get(at);
	        			maxDegreeAt = at;
	        		}
	        	}
	        	
	            //System.out.println( maxDegreeAt.getProperty(idProperty) + " " + punchingBag.getBondCount() );
	            
	            // C <-- C U {v}
	            cover.add( maxDegreeAt );

	            // remove from G' every edge incident on v, and v itself
	            //
	            Collection<Integer> neighbours = workBase.get( maxDegreeAt );
	            Set<Integer> keys = new HashSet<Integer>( workBase.keySet() );
	            for( Integer rowKey : keys ) {
	            	Collection<Integer> row = workBase.get(rowKey);
	            	row.removeAll( neighbours );
	            	row.remove( maxDegreeAt );
	            	
	            	if( row.isEmpty() ) {
	            		workBase.remove( rowKey );
	            	}
	            }
	            
	            
	           
	            workBase.remove( maxDegreeAt );
	            
	            //punchingBag.removeAtom(maxDegreeAt);
	            
	        }
	        
	        //if( verbose )	System.out.println( "cover verification = " + isVertexCover(g, cover) );

	        return cover;
	    }
		
		
		

		public static List<Integer> findRandomCover( Map<Integer, Collection<Integer>> workBase, int[][] adjList ) {
			
			Set<Integer> cover = new HashSet<Integer>();
			
			List<Integer> atomSet = new ArrayList<Integer>( workBase.keySet() );
	       Collections.shuffle( atomSet );  // the random part

	       for( Integer u : atomSet ) {
	    	   for( Integer v : workBase.get(u) ) {
	        	   if( ! cover.contains(v) ) {
	        		   cover.add( u );
	        	   }
	           }
	       }
	       
	       // refine cover
	       List<Integer> coverList = new ArrayList<Integer>( cover );
	       for( int i = 0; i < coverList.size(); i++ ) {
	    	   Integer vertex = coverList.remove(i);
	    	   
	    	   if( ! isVertexCover(adjList, coverList) )
	    		   coverList.add(vertex);
	       }
	       
	       
	       if( verbose )	System.out.println( "cover verification = " + isVertexCover(adjList, cover) + " | " + cover + " | " + adjList.length );
	       
	       return coverList;
			
		}
		
		
		/**
		 * Do the given list of atoms account for all the edges in the graph?
		 * 
		 * @param g
		 * @param cover
		 * @return
		 */
		protected static boolean isVertexCover( int[][] adjList, Collection<Integer> set ) {
			
			Set<String> totalBonds = new HashSet<String>();
			Set<String> coverBonds = new HashSet<String>();
			
			for( int n = 0; n < adjList.length; n++ ) {
				for( int j : adjList[n] ) {
					String edge = null;
					if( n < j ) {
						edge = n + "," + j;
					} else {
						edge = j + "," + n;
					}
					totalBonds.add(edge);
				}
			}
			
			
			for( Integer n : set ) {
				for( int j : adjList[n] ) {
					String edge = null;
					if( n < j ) {
						edge = n + "," + j;
					} else {
						edge = j + "," + n;
					}
					coverBonds.add(edge);
				}
			}
			
			
			if( coverBonds.size() == totalBonds.size() ) {
				return true;
			}
			
			return false;
		}
	
	
	/**
	 * initial refinement stage
	 * 
	 * @param g1
	 * @param g2
	 * @param matchList
	 * @param c1
	 * @param useF1NotC1
	 * @param labels
	 * @return
	 */
	private List< Integer > refineMatchCover( IAtomContainer g1, IAtomContainer g2, List<Integer> matchList, List<Integer> c1, boolean useF1NotC1 ) {
		
		List<Integer>  g1Matches = new ArrayList<Integer>( matchList.size() );
		//List<Integer>  g2Matches = new ArrayList<Integer>( matchList.length );
		
		for( int n = 0; n < matchList.size(); n++ ) {
			if( matchList.get(n) > 0 ) {
				g1Matches.add( n );
				//g2Matches.add( matchList[n] );
			}
		}
		
		//ArrayList<Integer> c1 = findRandomCover( g1 );
		//ArrayList<Integer> c2 = findRandomCover( g2 );
		
		Collection<Integer> f1 = new HashSet<Integer>( c1 );
		Collection<Integer> f2 = new HashSet<Integer>( c1.size() );
		f1.retainAll( g1Matches );
		
		//  construct F2 (corresponding matches from F1 in graph 2)
		for( int f : f1 ) {
			int correspondence = matchList.get( f );
			if( correspondence >= 0 )		f2.add( correspondence );
		}
		
		List<Integer> v1NotInC1 = new ArrayList<Integer>( g1.getAtomCount() );
		for( int n = 0; n < g1.getAtomCount(); n++ )	{ v1NotInC1.add( n ); }
		
		if( useF1NotC1 ) {
			v1NotInC1.removeAll( f1 );
		} else {
			v1NotInC1.removeAll( c1 );
		}
		
		List<Integer> v2NotInF2 = new ArrayList<Integer>( g2.getAtomCount() );
		for( int n = 0; n < g2.getAtomCount(); n++ )	{ v2NotInF2.add( n ); }
		v2NotInF2.removeAll( f2 );
		
		if( verbose )	System.out.println( c1 + " " + f1 + " " /*+ c2*/ + " " + f2 + " " + v1NotInC1 + " " + v2NotInF2 );
		
		
		// remove atoms from these lists which have no labels in common - avoids bogus matches
		/*if( labels ) {
			for( int v = 0; v < v2NotInF2.size(); v++ ) {
				IAtom vAt = g2.getAtom( v2NotInF2.get(v) );
				boolean commonLabel = false;
				
				for( int u = 0; u < v1NotInC1.size(); u++ ) {
					IAtom uAt = g1.getAtom( v1NotInC1.get(u) );
					if( uAt.getProperty( labelProperty ).equals( vAt.getProperty( labelProperty ) ) )
						commonLabel = true;
				}
				
				if( ! commonLabel ) {
					if( verbose )	System.out.println("removing " + vAt);
					v2NotInF2.remove(v);
					v--;
				}
			}
		}*/
		
		
		double[][] bigraphMatrix = new double[v1NotInC1.size()][v2NotInF2.size()];
		
		// construct weighted bipartite graph
		int uIndex = 0;
		for( int u : v1NotInC1 ) {
			
			IAtom atomU = g1.getAtom(u);
			
			
			Collection<Integer> nUF1 = new HashSet<Integer>( adj1[u].length );
			for( int neighbour : adj1[u] ) { nUF1.add(neighbour); }
			nUF1.retainAll( f1 );
			//System.out.println(nUF1 + " " + f1);
			
			Collection<Integer> mNuF1 = new ArrayList<Integer>( nUF1.size() );
			
			for( int f : nUF1 ) {
				int correspondence = matchList.get( f );
				if( correspondence >= 0 )		mNuF1.add( correspondence );
			}
			
			int vIndex = 0;
			for( int v : v2NotInF2 ) {
				
				// label-matching stage
				IAtom atomV = g2.getAtom(v);
				if( ! atomsMatch(g1, atomU, g2, atomV) ) {
					bigraphMatrix[uIndex][vIndex] = avoidSimilarityValue;
					continue;
				}
				
				Collection<Integer> nVF2 = new HashSet<Integer>( adj2[v].length );
				for( int neighbour : adj2[v] ) { nVF2.add(neighbour); }
				nVF2.retainAll( f2 );
				
				Collection<Integer> matchedEdgeContributions = new HashSet<Integer>( mNuF1 );
				matchedEdgeContributions.retainAll( nVF2 );
				
				bigraphMatrix[uIndex][vIndex++] = matchedEdgeContributions.size();
				
			}
			
			++uIndex;
		}
		
		DenseDoubleMatrix2D bigraphMat = new DenseDoubleMatrix2D( bigraphMatrix );
		
		
		HungarianAlgorithm ha = new HungarianAlgorithm( bigraphMat.assign( Functions.mult( -1 ) ).toArray() );
		int[] optimalMatches = ha.execute();
		
		// maximum weight matching - get indices of atoms
		ArrayList< SimPair > newMatched = new ArrayList< SimPair >( v1NotInC1.size() ); 
		for( int i = 0; i < optimalMatches.length; i++ ) {
			
			if( optimalMatches[i] >= 0 ) {
				int[] pair = new int[] {
						v1NotInC1.get(i),
						v2NotInF2.get( optimalMatches[i] )
				};
				
				//if( ! labels || g1.getAtom( matchPair.get(0) ).getProperty( labelProperty ).equals( g2.getAtom( matchPair.get(1) ).getProperty( labelProperty ) ) )
					newMatched.add( new SimPair( pair, 0.0 ) );
			}
			
		}
		
		
		Collection< SimPair > matchList2 = new HashSet< SimPair >( matchList.size() ); 
		for( int i = 0; i < matchList.size(); i++ ) {
			if( matchList.get( i ) >= 0 ) {
				matchList2.add( new SimPair( new int[]{ i, matchList.get( i ) }, 0.0 ) );
			}
		}
		
		Collection< SimPair > cartesianF1F2 = new HashSet< SimPair >( f1.size() * f2.size() );
		for( int i : f1 ) {
			for( int j : f2 ) {
				cartesianF1F2.add( new SimPair( new int[]{ i, j }, 0.0 ) );
			}
		}
		
		Collection< SimPair > refinedMatches = new HashSet< SimPair >( cartesianF1F2 );
		refinedMatches.retainAll( matchList2 );
		refinedMatches.addAll( newMatched );
		
		ArrayList<Integer> ml = new ArrayList<Integer>( refinedMatches.size() );
		for( int n = 0; n < g1.getAtomCount(); n++ ) { ml.add(-1); }
		for( SimPair p : refinedMatches ) { ml.set( p.getPair()[0], p.getPair()[1] ); }
		
		if( verbose )	{
			//System.out.println( "bigraphmat " + bigraphMat );
			System.out.println( "optimalMatches " + Arrays.toString(optimalMatches) );
			System.out.print( "newMatched (from Bipartite matching) " );
			for( int n = 0; n < newMatched.size(); n++ ) {
				//System.out.print( "(" + (newMatched.get(n).get(0) + 1) + g1.getAtom(newMatched.get(n).get(0)).getProperty( labelProperty ) + "," + (newMatched.get(n).get(1) + 1) + g2.getAtom(newMatched.get(n).get(1)).getProperty( labelProperty ) +"), " );
			}
			
			System.out.println("\n" + matchList);
			System.out.println("Original matches: " + matchList2);
			
			System.out.println( "F1 x F2: " + cartesianF1F2);
			System.out.println("refined: " + refinedMatches  );
			System.out.print("refined size: " + refinedMatches.size() + "  bonds - " + generateBondMap( ml ).size() );
			
			System.out.println();
			for( int n = 0; n < refinedMatches.size(); n++ ) {
				//System.out.print( "(" + (refinedMatches.get(n).get(0) + 1) + g1.getAtom(refinedMatches.get(n).get(0)).getProperty( labelProperty ) + "," + (refinedMatches.get(n).get(1) + 1) + g2.getAtom(refinedMatches.get(n).get(1)).getProperty( labelProperty ) +"), " );
			}
			System.out.println();
		}
		
		return ml;
	}

	
	/**
	 * removes unmapped atoms (and bonds) from consideration before generating the vertex cover
	 * 
	 * Relates to the "Randomly Refine Including C - F1" point in Zhu's work
	 * 
	 * @param g
	 * @param matchList
	 * @param coverSource2
	 * @return  adjacency list in the form of a map (suppliable to the cover generating thing)
	 */
	private static Map<Integer, Collection<Integer>> generateCoverSource( int[][] adjList, List<Integer> matchList, boolean coverSource2 ) {
		
		Map<Integer, Collection<Integer>> coverAdj = new HashMap<Integer, Collection<Integer>>();
		
		if( coverSource2 ) {
			for( int i = 0; i < adjList.length; i++ ) {
				if( matchList.contains(i) ) {
					Set<Integer> set = new HashSet<Integer>();
					for( int a : adjList[i] ) { 
						if( matchList.contains( a ) )
							set.add(a); 
					}
					
					coverAdj.put( i, set );
				}
			}
		} else {
			for( int i = 0; i < adjList.length; i++ ) {
				if( matchList.get(i) >= 0 ) {
					Set<Integer> set = new HashSet<Integer>();
					for( int a : adjList[i] ) { 
						if( matchList.get( a ) >= 0 )
							set.add(a); 
					}
					
					coverAdj.put( i, set );
				}
			}
		}
		
		return coverAdj;
		
		
	}
	
	
	
	/**
	 * Refinement of a given match list using randomly-generated minimal vertex covers
	 * 
	 * @param g1
	 * @param g2
	 * @param matchList
	 * @param includeC1MinusF1
	 * @return
	 */
	public List<Integer> coverRefine( IAtomContainer g1, IAtomContainer g2, List<Integer> matchList, int iterations, boolean includeC1MinusF1, boolean random  ) {
		
		
		List< Integer > tempML;
		
		int currentFitness = generateBondMap( matchList ).size();
		List<Integer> cover = null;
		
		Map<Integer, Collection<Integer>> cover1Source = null, cover2Source = null;
		
		
		
		
		
		for( int it = 0; it < iterations; it++ ) {
			
			// Exclude unmatched nodes from the cover (we want to match these later if possible hence why we're excluding them)
			if( includeC1MinusF1 ) {
				cover1Source = generateCoverSource(adj1, matchList, false);
				cover2Source = generateCoverSource(adj2, matchList, true);
			} else {
				cover1Source = new HashMap<Integer, Collection<Integer>>();
				for( int n = 0; n < adj1.length; n++ ) {
					Collection<Integer> set = new HashSet<Integer>();
					for( int s : adj1[n] ) { set.add(s); }
					cover1Source.put(n, set);
				}
				
				cover2Source = new HashMap<Integer, Collection<Integer>>();
				for( int n = 0; n < adj2.length; n++ ) {
					Collection<Integer> set = new HashSet<Integer>();
					for( int s : adj2[n] ) { set.add(s); }
					cover2Source.put(n, set);
				}
			}
			
			if( random ) {
				if( it % 2 == 0 ) {
					//cover = findGreedyCover(cover1Source);
					cover = findRandomCover(cover1Source, adj1);
				} else {
					//cover = findGreedyCover(cover2Source);
					cover = findRandomCover(cover2Source, adj2);
				}
			} else {
				if( it % 2 == 0 ) {
					if( cover1Source.isEmpty() )
						continue;
						
					cover = findCoverIsolationAlgorithm(cover1Source, adj1);
				} else {
					if( cover2Source.isEmpty() )
						continue;
					
					cover = findCoverIsolationAlgorithm(cover2Source, adj2);
				}
			}
			 
			if( ! cover.isEmpty() ) {
				tempML = refineMatchCover(g1, g2, matchList, cover, includeC1MinusF1);
			} else {
				tempML = matchList;
			}
			
			// score solution
			int newScore = generateBondMap( tempML ).size();
			
			if( newScore > currentFitness ) {
				matchList = tempML;
				currentFitness = newScore;
				it = -1;
				if( verbose )	System.out.println("new score " + newScore);
			}
		}
		
		if( verbose )	System.out.println( "BEST: " + matchList + " " + currentFitness );
		
		return matchList;
		
	}
	
	
	
	/**
	 * Maximum bipartite matching of a matrix, returns the highest-scoring pairs
	 * 
	 * @param simMat
	 * @return
	 */
	public List<int[]> similarityMatchList( DoubleMatrix2D simMat ) {
		
		List<int[]> matchList = new ArrayList<int[]>( simMat.rows() );
		DoubleMatrix2D bigraphMatrix = simMat.copy();
		
		if( simMat == null || simMat.size() == 0 )
			return matchList;
		
		
		
		HungarianAlgorithm ha = new HungarianAlgorithm(bigraphMatrix.assign( Functions.mult(-1) ).toArray() );
		int[] matches = ha.execute();
		
		// obtain the best pairs - convert array to ArrayList of pairs
		for( int n = 0; n < matches.length; n++ ) {
			
			if( matches[n] < 0 )
				continue;
			
			//int[] matchPair = new int[]{ n, matches[n] };
			int[] matchPair = new int[]{ n, matches[n] };
			matchList.add( matchPair );
		}
		
		//System.out.println( bigraphMatrix );
		
		
		return matchList;
	}
	
	
	/**
	 * Maximum bipartite matching of a matrix, returns the highest-scoring pairs
	 * 
	 * @param simMat
	 * @return
	 */
	public List<int[]> similarityMatchList( DoubleMatrix2D simMat, Collection<Integer> rows, Collection<Integer> cols ) {
		
		int[] nlRows = new int[ rows.size() ];
		int[] nlCols = new int[ cols.size() ];
		
		int r = 0;
		for( Integer n : rows ) {
			nlRows[r++] = n;
		}
		
		int c = 0;
		for( Integer n : cols ) {
			nlCols[c++] = n;
		}
		
		
		List<int[]> optimalPairs = similarityMatchList( simMat.viewSelection( nlRows, nlCols ) );
		
		List<int[]> optimalPairsTranslated = new ArrayList<int[]>( optimalPairs.size() );
		for( int[] p : optimalPairs ) {
			optimalPairsTranslated.add( 
					new int[]{ 
							nlRows[ p[0] ], 
							nlCols[ p[1] ] 
			} );
		}
		
		//printMatches( optimalPairsTranslated );
		return optimalPairsTranslated;
	}
	
	
	
	private Map<IBond, IBond> generateBondMap( List<Integer> matchList ) {
		// output common substructure
					Map<IAtom, IAtom> commonAtomMap = new HashMap<IAtom, IAtom> ();
					for(int n = 0; n < matchList.size(); n++ ) {
						
						if( matchList.get(n) < 0 )
							continue;
						
						IAtom hsElem = hsMol.getAtom( n );
						IAtom qElem = qMol.getAtom( matchList.get(n) );
						
						commonAtomMap.put( hsElem, qElem );
					}
					
					Map<IBond, IBond> commonBondMap = ConvenienceTools.makeBondMapOfAtomMap(hsMol, qMol, commonAtomMap, matchBonds);
					
					return commonBondMap;
	}
	
	
	@Override
	public void search(IAtomContainer graph1, IAtomContainer graph2) {
		
		
		
		
		// adjacency lists
		/*adj1 = ConvenienceTools.createBondAdjacencyList(graph1);
		adj2 = ConvenienceTools.createBondAdjacencyList(graph2);*/
		bMap1 = GraphUtil.EdgeToBondMap.withSpaceFor(graph1);
		bMap2 = GraphUtil.EdgeToBondMap.withSpaceFor(graph2);
		adj1 = GraphUtil.toAdjList(graph1, bMap1);  // atoms
		adj2 = GraphUtil.toAdjList(graph2, bMap2);
		
		if( verbose )	System.out.println( "Mapping graph 2 to graph 1..." );

		
		long gSimTime = System.currentTimeMillis();
		
		/*
		 * XXX  Java Colt API used for eigendecomposition
		 */
		DoubleMatrix2D globalSim = globalSimilarity( graph1, graph2 );
		if( verbose ) {	
			System.out.println( "globalSim: " + globalSim );
			System.out.println( System.currentTimeMillis() - gSimTime );
		}

		
		long lSimTime = System.currentTimeMillis();
		DoubleMatrix2D localSim = localSimilarityAtoms(graph1, graph2, neighbourhoodOrder);
		
		if( verbose )
			System.out.println( System.currentTimeMillis() - lSimTime );
		
		
		//DoubleMatrix2D localSim = globalSim;
		if( verbose )	System.out.println( "LocalSim: " + localSim );
		/*localSim = new DenseDoubleMatrix2D( new double[][]{{0.94, 0.6248077, 0.7120253, 0.7111688, 0.8234722, 0.4385965, 0.4308511, 0.8910843, 0.1511429, 0.8692063, 0.4907547, 0.5207143, 0.6469231, 0.4709302, 0.197027, 0.2288095, 0.1931429},
				{0.6107018, 0.9219048, 0.7555814, 0.8359184, 0.6912, 0.7020755, 0.7601282, 0.5485714, 0.2484507, 0.6315789, 0.7912329, 0.7823256, 0.7475556, 0.8012162, 0.3232692, 0.3875439, 0.3281818},
				{0.847, 0.6864286, 0.7912329, 0.7612857, 0.8121918, 0.4947059, 0.4708696, 0.8066667, 0.1666667, 0.8125, 0.4780328, 0.5217391, 0.676875, 0.5501887, 0.225, 0.2644898, 0.2278125},
				{0.7412281, 0.9155769, 0.8748, 0.9443902, 0.7953731, 0.6321951, 0.6384091, 0.680625, 0.2118919, 0.7125397, 0.648, 0.710597, 0.8032787, 0.6666667, 0.2849153, 0.3259524, 0.2831373},
				{0.8786441, 0.726, 0.7301408, 0.7623529, 0.961, 0.4998462, 0.5001923, 0.7876923, 0.1704348, 0.8767123, 0.550678, 0.5911364, 0.7383051, 0.5501887, 0.21875, 0.2656098, 0.2181132},
				{0.492093, 0.7020755, 0.6706349, 0.6801515, 0.5525397, 0.8691429, 0.8691379, 0.4609524, 0.3133898, 0.5142857, 0.7140984, 0.85, 0.6701538, 0.8692063, 0.4090741, 0.4843103, 0.414902},
				{0.6942857, 0.907561, 0.81, 0.9230769, 0.7593878, 0.6505797, 0.66, 0.6410256, 0.2181132, 0.7324615, 0.6933333, 0.735, 0.8301235, 0.7202381, 0.2898276, 0.3438095, 0.2934921},
				{0.8704819, 0.5976596, 0.6567925, 0.66, 0.7837288, 0.4108889, 0.4119149, 0.9216, 0.135641, 0.8448276, 0.441, 0.4767568, 0.6008889, 0.44, 0.1858065, 0.2180645, 0.1824138},
				{0.2016667, 0.30625, 0.28, 0.2803846, 0.2336111, 0.4530189, 0.4408333, 0.1959259, 0.735, 0.215641, 0.3944643, 0.3834783, 0.2989831, 0.4232, 0.9941379, 0.8421053, 1.001739},
				{0.9383117, 0.6506557, 0.6975294, 0.7092045, 0.8210227, 0.4389063, 0.448, 0.855625, 0.1451852, 0.9210465, 0.48, 0.5236364, 0.6506557, 0.4726563, 0.1955814, 0.2336111, 0.1912195},
				{0.70875, 0.8751515, 0.7505063, 0.8511392, 0.8010465, 0.6367442, 0.634717, 0.6272, 0.2101818, 0.7305195, 0.7118519, 0.7510345, 0.9094737, 0.6901408, 0.2891803, 0.3333333, 0.2817391},
				{0.5718919, 0.7585185, 0.7025806, 0.7442, 0.6515254, 0.7967368, 0.7837288, 0.5008696, 0.271129, 0.5851429, 0.8780488, 0.9411111, 0.86, 0.7935, 0.3562903, 0.4108889, 0.3491379},
				{0.6272, 0.8915714, 0.7735714, 0.8178723, 0.7161017, 0.6727119, 0.7185965, 0.5377778, 0.2414286, 0.6433784, 0.8069492, 0.8066667, 0.8816667, 0.7213636, 0.324386, 0.3723077, 0.316875},
				{0.29, 0.4302222, 0.4026471, 0.3911429, 0.3281818, 0.55, 0.6223529, 0.2694737, 0.5259524, 0.2943243, 0.567, 0.5308163, 0.421875, 0.5916981, 0.7074419, 0.8205128, 0.6993103},
				{0.197027, 0.3022857, 0.2852083, 0.2793878, 0.2266667, 0.4517647, 0.4545455, 0.1889286, 0.7401316, 0.2048, 0.405, 0.3716129, 0.2946939, 0.4153086, 0.9906349, 0.8501149, 0.9906349},
				{0.2336111, 0.3637879, 0.3266667, 0.3265306, 0.2625455, 0.5140984, 0.5215517, 0.2232143, 0.648, 0.24, 0.4761017, 0.4418, 0.3552632, 0.4928205, 0.8592754, 0.9959211, 0.8569014}
				} );*/
		
		DoubleMatrix2D finalSim = globalSim.copy().assign( localSim, Functions.mult );
		
		
		/* XXX  An adjustment that Zhu et al used (not mentioned in their publications)
		 * 
		 * EDIT:  Removed.  Seems this worsens MCS size, on empirical observation, at least with
		 * the local similarity in its current state.
		 */
		//finalSim.assign( Functions.sqrt );  
		
		if( verbose )	System.out.println( "FinalSim: " + finalSim );
		
		if( verbose ) {
			List<int[]> simMatchList = similarityMatchList( finalSim );
			printMatches( simMatchList );
		}
		
		List<Integer> coverRefined = null;
		
		for( int an = 1; an <= anchorLimit; an++ ) {
			
			Map<IBond, IBond> commonBondMap = new HashMap<IBond, IBond>();
			Map<IAtom, IAtom> commonAtomMap = new HashMap<IAtom, IAtom>();
			
			List<SimPair> anchors = findAnchors(graph1, graph2, null, null, finalSim, an, true);
			//int[] expanded = expandAnchors(graph1, graph2, localSim, anchors);
			List<Integer> expanded = expandAnchors(graph1, graph2, localSim, anchors);
			//int[] exp2 = new int[]{0, 6, 2, 3, 4, 5, 1, 7, 8, 9, 12, 11, 10, 13, 14, 15};
			//expanded.clear();
			//for( int n : exp2 ) { expanded.add(n); } 
			int expandedBonds = generateBondMap(expanded).size();
			if( verbose )	System.out.println( "Expanded " + expanded + " " + expandedBonds  );
			
			
			if( verbose )	System.out.println( "Cover refinement..." );
			coverRefined = expanded;
			
			
			
			
			/* XXX  Initial non-random cover refinement performed here seeing as this Isolation algorithm 
			 * is pretty good at approximating a minimum vertex cover
			 */
			
			
			coverRefined = coverRefine( hsMol, qMol, coverRefined, 2, false, false );
			coverRefined = coverRefine( hsMol, qMol, coverRefined, 2, true, false );
			//coverRefined = refineMatchCover( graph1, graph2, expanded, findCoverIsolationAlgorithm( cover2Source, adj2 ), true );
			int cRefinedBonds = generateBondMap(coverRefined).size();
			
			if( verbose )	System.out.println( "First cover size - " + cRefinedBonds + " vs expanded size - " + expandedBonds );
			if( expandedBonds > cRefinedBonds )
				coverRefined = expanded;
			
			/*
			 * XXX 
			 * Decided to avoid random refinement due to it reducing the reproducibility of an MCS
			 */
			/*if( verbose )	System.out.println( "Random Cover refinement..." );
			coverRefined = coverRefine( hsMol, qMol, coverRefined, refinementIterations, false, true );  // refine excluding V - C
			coverRefined = coverRefine( hsMol, qMol, coverRefined, refinementIterations, true, true );  // refine including V - C
			*/
			
			if( verbose )	System.out.println( coverRefined );
			
			
			
			// output common substructure
			for (int n = 0; n < coverRefined.size(); n++) {

				if (coverRefined.get(n) < 0)
					continue;

				IAtom hsElem = hsMol.getAtom(n);
				IAtom qElem = qMol.getAtom(coverRefined.get(n));

				// unfortunately some of the heuristics for unknown rings allow false matchings - don't allow them
				if( ConvenienceTools.atomsMatch( hsElem, qElem ) )
					commonAtomMap.put(hsElem, qElem);
			}

			commonBondMap = ConvenienceTools.makeBondMapOfAtomMap(graph1, graph2, commonAtomMap, matchBonds);
			
			if( mcsBondIsomorphisms.size() > 0 && 
				commonBondMap.size() > mcsBondIsomorphisms.get( 0 ).size() ) {
				mcsBondIsomorphisms.clear();
				mcsAtomIsomorphisms.clear();
			}

			if( mcsBondIsomorphisms.isEmpty() || 
					commonBondMap.size() >= mcsBondIsomorphisms.get( 0 ).size()) {
				mcsBondIsomorphisms.add(commonBondMap);
				mcsAtomIsomorphisms.add(atomMapToChromosome(hsMol, qMol, commonAtomMap));
				
				// create indices
				List<int[]> pairList = new ArrayList<int[]>( commonBondMap.size() );
				for( Entry<IBond, IBond> e : commonBondMap.entrySet() ) {
					int[] pair = new int[]{ hsMol.getBondNumber(e.getKey()), qMol.getBondNumber(e.getValue()) };
					pairList.add(pair);
				}
				mcsBondIndexIsomorphisms.add( pairList );
			}
		}

		
						
		

	}
	
	
	private void printMatches( List<int[]> matches ) {
		
		for( int[] m : matches ) {
			System.out.println( m[0] + "," + m[1] + "  " + hsMol.getAtom(m[0]).getSymbol() + " " + qMol.getAtom(m[1]).getSymbol() );
		}
		
	}
	
	
	public static void main( String[] args ) {
		
		String inputFileName = compoundPath + "/mos_mapping_test4.smi";
		//String inputFileName = compoundPath + "/zhu_graphs/zhu.smiles";
		//String inputFileName = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/hyperstructure_query_test1.smi";
		//String inputFileName = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/mddr/mddr_cox_random_10.sdf";
		//String inputFileName = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/aid466_actives_2.sdf";
		
		
		
		
		ArrayList<IAtomContainer> graphs = null;
		try {
			graphs = ConvenienceTools.getQueryMolecules( new File( inputFileName ), null );
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//for( IAtomContainer gr : graphs )	 assignLabels( gr );

		
		/*System.out.println( "graph 1: " + sGenerator.create(graphs.get(0)) );
		System.out.println( "graph 2: " + sGenerator.create(graphs.get(1)) );*/
		
		
	}
	
	/*protected List<List<Integer>> adj1;
	protected List<List<Integer>> adj2;*/
	int[][] adj1, adj2;
	GraphUtil.EdgeToBondMap bMap1, bMap2;
	int neighbourhoodOrder = 2;
	int anchorLimit = 3;
	//int refinementIterations = 6;
	
	protected int avoidSimilarityValue = -9999;
	
	protected static Algebra algebra = new Algebra();
	protected static String compoundPath = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/";
	protected static boolean verbose = false;
}
