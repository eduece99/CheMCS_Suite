package org.cisrg.mapping;


import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import org.cisrg.BitSetExtended;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;


/*
import com.chemaxon.search.mcs.McsSearchResult;
import com.chemaxon.search.mcs.maxclique.MaxCliqueMcs;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.doublealgo.Transform;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.jet.math.Functions;
import chemaxon.core.calculations.FindAllRings;
import chemaxon.struc.Molecule;
*/

public class DepolliCliqueDetection extends CliqueDetection {
	
	
    // stable sort - should retain original order on ties (hence no 0 handling)
	private Comparator<int[]> descendingDegreeComparator = new Comparator<int[]>() {
		
		public int compare(int[] o1, int[] o2) {
			if( o1[vDegree] == o2[vDegree] )
				return 0;
			
            return o1[vDegree] < o2[vDegree] ? 1 : -1;
		}
	 };
	 
	 
	 protected Comparator<int[]> descendingExDegreeComparator = new Comparator<int[]>() {
		 public int compare(int[] o1, int[] o2) {

	            /*int dComp;
	            
	            // sort in descending order of degree first
	            if( o1[vExDegree] == o2[vExDegree] ) {
	            	dComp = 0;
	            } else if( o1[vExDegree] > o2[vExDegree] ) {
	            	dComp = -1;
	            } else {
	            	dComp = 1;
	            }*/
	            
	            if( o1[vExDegree] == o2[vExDegree] ) {
	            	return 0;
	            }
	            
	            return o1[vExDegree] < o2[vExDegree] ? 1 : -1;
	            //return dComp;
		 }
	 };
	 
	

	public DepolliCliqueDetection(ModularProductOptions opts) {
		super( opts );
	}
	
	

	
	 
	protected LinkedList<int[]> assignColours( Collection<Integer> vertices, int maxSize, int cliqueSize ) {  // maxSize defaults to 0
		//List<int[]> colouredVertices = new ArrayList<int[]>( orderedVertices );
		int colour = 1;
		int kMin = maxSize - cliqueSize;
		
		BitSetExtended<Integer> vertices1 = new BitSetExtended<>( vertices );
		BitSetExtended<Integer> verticesCopy = new BitSetExtended<>( vertices );
		LinkedList<int[]> colourSorted = new LinkedList<>(   );
		//List<int[]> colourSorted = new ArrayList< >( vertices.size() );

		//System.out.println( "Colouring, LinkedHashSet creation TIME - " + (System.currentTimeMillis() - time1) );
		
		
		
		while( ! vertices1.isEmpty() ) {
			
			while( ! verticesCopy.isEmpty() ) {
				
				// get first element that's not null in the copy
				//int v = verticesCopy.iterator().next();
				int v = verticesCopy.getFirstBit();
				
				// retainAll (and) operation on the prepared inverse matrix is quicker than removeAll (andNot) on original adjMatrix
				verticesCopy.retainAll( adjListInverse.get(v) );
				
				verticesCopy.remove(v);
				
				//System.out.println( vertices.get(0) + " " + index  + " " + vertices1.contains(v) );
				vertices1.remove( v );
				

				//System.out.println( "size = " + vertices1.size() + " " + verticesCopy.size() + " | colour = " + colour  );
				
				
				//System.out.println( verticesCopy.size() );
				
				if( colour > kMin ) {
					/*int[] vertex = orderedVertexList.get(v);
					//vertex[vColour] = Math.min(colour, vertex[vColour]);
					vertex[vColour] = colour + 3;
					colourSorted.add(vertex);*/
					colourSorted.add( new int[]{ v, 0, 0, colour } );
				}
				
			}
			
			
			// vcopy = v
			verticesCopy = new BitSetExtended<>( vertices1 );
						
			++colour;
			
			/*
			for( int e = 0 ; e < orderedVertices.size(); e++ ) {
				int[] ea = orderedVertices.get(e);
				System.out.println( e + " [ " + ea[0] + "," + ea[1] + "," + ea[2] + "," + ea[3] + "," + ea[4] + "]"  );
			}
			*/
			
			
			
			
			
			
			/*Iterator<int[]> it = vertices1.iterator() ;
			while( it.hasNext() ) {
				verticesCopy.add( it.next() );
			}*/
			
		}
		
		/*System.out.println( "kMin " + maxSize );
		for( int[] ea : colourSorted ) {
			System.out.println( "b [ " + ea[vIndex] + "," + ea[vDegree] + "," + ea[vExDegree] + "," + ea[vColour] + " ]"  );
		}
		System.out.println();*/
		
		return colourSorted;
	}
	
	
	 
	protected void expand( Set<Integer> vertices, LinkedList<int[]> colourRestricted, Collection<int[]> clique  ) throws CDKException {
		
		++numberOfSteps;
		
		if( numberOfSteps % 10000 == 0 ) {
			long currentTime = System.currentTimeMillis();
			if( currentTime - mcsStartTime > expansionTimeLimitMs ) {
				throw new CDKException("Error - time limit of " + expansionTimeLimitMs + " ms reached!");
			}
		}
		
		//while( ! colourRestricted.isEmpty() && numberOfSteps < maxSteps ) {
		while( ! colourRestricted.isEmpty() ) {
			
			if( verbose )
				System.out.println( "counter: " + numberOfSteps + " vertices - " + vertices.size() + " " + colourRestricted.size() + " " + clique.size() + " maxCliqueSize - " + cliqueMax.size()  );
			 
			 
			
			 //if ( clique.size() + v[4] <= cliqueMax.size() ) {return;}
			 
			 // we don't want cliques of the same size or smaller to the current maximum
			 // base this on the number of remaining colours to use
			// if ( clique.size() + colourRestricted.get( colourRestricted.size() - 1 )[vColour] <= cliqueMax.size() ) {return;}
			 if ( clique.size() + colourRestricted.getLast()[vColour] <= cliqueMax.size() ) {return;}
			 
			 //int vcrIndex = colourRestricted.size() - 1;  // last vertex
			 //int[] v = colourRestricted.get( vcrIndex );  
			 int[] v = colourRestricted.getLast();
			 
			 
			 /*
			  * direct edge constraint:
			  * 
			  * - new edge must be connected to at least one edge in the clique
			  * 
			  * in other words:
			  * - if a vertex is an edge pair, both edges in the edge pair-to-be-added must be both edges of another edge pair in the clique
			  */
			/*boolean directEdgeConstraint = false;
			 for( int[] ce : clique ) {
				 
				 // See if considered vertex is a neighbour to both molecule's edges, for at least one clique node
				 System.out.println( "\n" + v[vIndex] + " " + ce[vIndex] + qMol.getBond( v[1] ) );
				 
				 int[] trv1 = modProd.getNodes().get( vertexTranslations.get( ce[vIndex] ) );
				 int[] trv2 = modProd.getNodes().get( vertexTranslations.get( v[vIndex] ) );
				 
				 if( hsMol.getBond( trv1[0] ).isConnectedTo( hsMol.getBond( trv2[0] ) ) && 
						 qMol.getBond( trv1[1] ).isConnectedTo( qMol.getBond( trv2[1] ) ) ) {
					 directEdgeConstraint = true;
					 System.out.println("connected " + clique.size() );
					 break;
				 }
			 }
			 
			 if( clique.size() > 0 && ! directEdgeConstraint )
				 return;*/
			 
			 clique.add(v);
			 
			 Set<Integer> neighbourIntersected = new BitSetExtended<Integer>( vertices );
			 //Set<Integer> neighbours = Arrays.asList( orderedAdjListMap.get(v) );
			 neighbourIntersected.retainAll( adjList.get( v[vIndex] ) );  // intersection
			 //neighbourIntersected.remove(v[vIndex]);
			 
			 if( verbose ) {
				 System.out.print( "clique: " + numberOfSteps + " "  );
				 for( int[] e : clique ) {
					 System.out.print( "[" + e[0] + "," + e[1] + "," + e[3] + "], " );
				 }
				 System.out.print( "\n" + "remaining: " + colourRestricted.size() + " " + neighbourIntersected.size() + "  " );
				 
				 
				 if( vertices.size() < 99 ) {
					 for( int e : vertices ) {
						 int[] ea = orderedVertexList.get(e);
						 System.out.print( "[" + ea[0] + "," + ea[1] + "], " );
					 }
					 System.out.println("");
				 }
			 }
			 
			 
			 if( neighbourIntersected.isEmpty() ) {
				 if ( clique.size() > cliqueMax.size() ) { 
					 
					 if( deltaYPossible ) {
						 List<Integer> cliqueTr = new ArrayList<Integer>( clique.size() );
						 
						 for( int[] cV : clique ) {
							 cliqueTr.add( vertexTranslations.get( cV[vIndex] ) );
						 }
				 
						 /* if( ! ConvenienceTools.deltaYExchangeOccured(modProd, hsMol, qMol, cliqueTr) )
						 	cliqueMax = new ArrayList<int[]>( clique );*/
						 
					 
						 if( ! ConvenienceTools.deltaYExchangeOccured(modProd, hsMol, qMol, cliqueTr) )
							 cliqueMax = new ArrayList<int[]>( clique );
					 } else {
						 cliqueMax = new ArrayList<int[]>( clique );
					 }	
				 }
			 } else {
				 //Set<int[]> cRestricted2 = getRelevantColouredVertices(neighbourIntersected, cliqueMax.size() - clique.size() );
				 LinkedList<int[]> cRestricted2 = assignColours( neighbourIntersected, cliqueMax.size(), clique.size() );
				// Set<int[]> cRestricted2 = assignColoursNeighbourhood( new ArrayList<int[]>( vertices ), cliqueMax.size() - clique.size() );
				 expand( neighbourIntersected, cRestricted2,  clique );
			 }
			 
			 int vIndexToRemove = colourRestricted.getLast()[vIndex];
			 clique.remove(v);  
			 colourRestricted.removeLast();
			 vertices.remove( vIndexToRemove );
		}
		
	}
	
	
	/**
	 * Give a set of vertices, return those (and their neighbours) of suitable minimum colour (max clique - current clique size)
	 * 
	 * @return
	 */
	protected Set<int[]> getRelevantColouredVertices( Set<int[]> inputVertices, int minColour ) {
		
		Set<int[]> vertexSet = new LinkedHashSet<int[]>( inputVertices.size() );
		
		for( int[] v : inputVertices ) {
			if( v[4] > minColour ) {
				vertexSet.add(v);
			}
		}
		
		return vertexSet;
	}
	
	
	/**
	 * Unsure if this' actually correct, but it seems to work
	 * 
	 * @param node
	 * @param adjList
	 * @return
	 */
	protected int exDegree( int node, List<Collection<Integer>> adjList ) {
		
		int exDegree = 0;
		
		
		//exDegree = adjList.get(node).size() * 2;  // hack to make it on-par with original source code results
			
		for( Integer neighbour : adjList.get(node) ) {
				exDegree += adjList.get(neighbour).size();
		}
		
		return exDegree;
	}

	
	

	
	protected void findCliquesDP( List<int[]> nodes, List<Collection<Integer>> edgeList, IAtomContainer graph1, IAtomContainer graph2 ) {
		
		
		
		hsMol = graph1;
		qMol = graph2;
		
		cliqueMax = new LinkedHashSet<int[]>();
		
		int nNodes = edgeList.size();
		
		
		// adjacency List initialisation
		this.adjList = edgeList;
		// initialise inverse adj list
		adjListInverse = new ArrayList<Collection<Integer>>( nNodes );
		for (int i = 0; i < nNodes; ++i) {
			BitSetExtended<Integer> origBitSet = (BitSetExtended<Integer>) adjList.get(i);
            BitSet inversed = (BitSet) origBitSet.getBitSet().clone();
            inversed.flip( 0, nNodes );
            adjListInverse.add( new BitSetExtended<Integer>( inversed ) );
        }
		
		
		// vertex list initialisation
		orderedVertexList = new ArrayList<int[]>( nNodes );
		
		/*
		double[][] simMatrix1 = null;
		if( graph1 != null && graph2 != null) {
			//DoubleMatrix2D simMatrixPart1 = globalBondSimilarity( graph1, graph2 )  ;
			//DoubleMatrix2D simMatrixPart1 = partialChargeSimilarity( graph1, graph2 )  ;
			DoubleMatrix2D simMatrixPart2 = DepolliCliqueDetectionHashes2.localPathSimilarity2( graph1, graph2, 3, false );  // could be sped up
			//DoubleMatrix2D simMatrixPart2 = partialChargeSimilarity( graph1, graph2 );
			//DoubleMatrix2D simMatrixPart4 = localPathSimilarity2( graph1, graph2, 3, true );  // could be sped up
			//simMatrixPart2 = simMatrixPart2.assign( simMatrixPart4, Functions.mult );
			
			DoubleMatrix2D simMatrixPart3 = DepolliCliqueDetectionHashes2.globalBondSimilarity( graph1, graph2 )  ;
			
			System.out.println( simMatrixPart2 );
			
			//DoubleMatrix2D temp = simMatrixPart1.copy().assign( simMatrixPart2, Functions.mult );
			//simMatrix1 = temp.copy().assign( simMatrixPart3, Functions.mult ).toArray();
			simMatrix1 = simMatrixPart3.copy().assign( simMatrixPart2, Functions.mult ).toArray();
			//simMatrix1 = simMatrixPart2.toArray();
			//simMatrix1 = partialChargeSimilarity( graph1, graph2 ).toArray()  ;
		}
		*/
		
		
		
		// Attempt to use minimum vertex cover to bias search space.  Didn't work.
		/*int[][] adjListTemp = ConvenienceTools.listOfListsToMatrix(adjList);
		Map<Integer, Collection<Integer>> adjListMapTemp = ConvenienceTools.adjListToMap( adjListTemp );
		List<Integer> cover = ZhuSpectralMCES.findCoverIsolationAlgorithm(adjListMapTemp, adjListTemp);
		*/
		
		

		for( int i = 0; i < nNodes; i++ ) {
			int degree = adjList.get(i).size();  // degree
			
			/*if( cover.contains(i) ) {
				degree += adjList.size();
			}*/
			
			int[] vertex = new int[4];
			vertex[vIndex] = i;  // index
			vertex[vDegree] = degree;  // degree
			//vertex[vDegree] = (int) (10000 * simMatrix1[ modProd.getNodes().get(i)[0] ][ modProd.getNodes().get(i)[1] ] );  // degree
			vertex[vExDegree] = exDegree( i, adjList );  // used for ex-degree
			orderedVertexList.add(vertex);
		}
		
		
		// create a vertex mapping table that will be used to renumber vertices back to original
        vertexTranslations = new ArrayList<Integer>( nNodes );
        
        // create default mapping that maps i → i
        for (int i = 0; i < nNodes; ++i)
        	vertexTranslations.add( i );
        
        
		
		
        
        
		/*double[][] simMatrix1 = null;
		int sleepTime = 1;*/
		
		//long startTime = System.currentTimeMillis();
		orderingTime = System.currentTimeMillis();
		
		
		//orderVerticesByTwoDegreeTypes( nodes, adjList, null );
		
		if( performOrdering )
			initialSort();
		else {
			try {
				//Collections.sort( orderedVertexList, descendingDegreeComparator );
				
				List<Integer> vertexTranslations2 = new ArrayList<>( nNodes );
			        
			        // create default mapping that maps i → i
			        for (int i = 0; i < nNodes; ++i)
			        	vertexTranslations2.add( orderedVertexList.get(i)[0] );
				
				orderVertices(vertexTranslations);
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		orderingTime = System.currentTimeMillis() - orderingTime;
		
		expansionTime = System.currentTimeMillis();
		
		// dud vertex list (all vertices present)
		ArrayList<Integer> initialSet = new ArrayList<Integer>( orderedVertexList.size() );

		for( int n=0; n < orderedVertexList.size(); n++ ) {
			initialSet.add(n);
		}
		
		
		try {
			if( performOrdering ) {
				expand( 
						new BitSetExtended<Integer>( initialSet ), 
						//assignColours( new BitSetCollectionAL<Integer>( initialSet ), 0 ), 
						new LinkedList<int[]>( orderedVertexList ),
						//new ArrayList<int[]>( assignColours( initialSet, cliqueMax.size(), 0 ) ),
						new LinkedList<int[]>()
						//cliqueMax
				);
			} else {
				expand( 
						new BitSetExtended<Integer>( initialSet ), 
						new LinkedList<int[]>( assignColours( initialSet, cliqueMax.size(), 0 ) ),
						new LinkedList<int[]>()
				);
			}
		} catch (CDKException e) {
			e.printStackTrace();
		}
		
		expansionTime = System.currentTimeMillis() - expansionTime;
		
		List<Integer> clique = new ArrayList<Integer>( cliqueMax.size() );
		
		
		for( int[] cV : cliqueMax ) {
			clique.add( vertexTranslations.get( cV[vIndex] ) );
		}
		Collections.sort( clique );
		
		if( verbose ) {
			System.out.println( "\n clique size - " + cliqueMax );
			System.out.print( "clique members - " );
			for( int[] v : cliqueMax ) { System.out.print( v[vIndex] + ", " ); }
			System.out.println( "\nordering time " + orderingTime );
			System.out.println( "expansion time " + expansionTime );
			System.out.println( "expansions - " + numberOfSteps );
			System.out.println( "translated clique - " + clique );
		}
		
		/*int n = 0;
		for( int[] cv : orderedAdjListMap.keySet() ) {
			if( cliqueMax.contains( cv ) )
				clique.add(n);
			
			n++;
		}
		*/
		
		bestCliques.add( clique );
	
		
	}
	
	
	
	@Override
	protected void findCliques() {
		
		List<int[]> nodes;
		List<Collection<Integer>> edges;
		
		/*nodes = new ArrayList<int[]>(6);
		for( int n = 0; n < 6; n++ ) {
		  nodes.add( new int[]{ n,n } );
		}

		edges = new ArrayList<Collection<Integer>>(6);
		edges.add( Arrays.asList( new Integer[]{ 1,2,3,4,5 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0, 2,3,4 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0,1, 3 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0,1,2 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0,1 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0 } ) );*/
		
		/*
		nodes = new ArrayList<int[]>(5);
		for( int n = 0; n < 4; n++ ) {
		  nodes.add( new int[]{ n,n } );
		}

		edges = new ArrayList<Collection<Integer>>(6);
		edges.add( Arrays.asList( new Integer[]{ 1,2,4 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0,2,4 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0,1,3,4 } ) );
		edges.add( Arrays.asList( new Integer[]{ 2 } ) );
		edges.add( Arrays.asList( new Integer[]{ 0,1,2 } ) );*/
		
		nodes = modProd.getNodes();
		edges = modProd.getAdjacencyList();
		
		//ConvenienceTools.writeDIMACSGraph( modProd );
		
		findCliquesDP( nodes, edges, hsMol, qMol );
		//findCliquesDP( modProd.getNodes(), modProd.getAdjacencyList(), hsMol, qMol );
	}
	
	
 
	
	
	protected void initialSort() {
		
		int vSize = adjList.size();
		
		//adjList = modProd.getAdjacencyList();
		int maxDegree = adjList.get(0).size();
		//Collection<Integer> r = new ArrayList<Integer>( vSize );
		List<Integer> vOrder = new ArrayList<Integer>( vSize );  // re-ordering vector
		List<Integer> storedColour = new ArrayList<Integer>( vSize );  // stores colours
		List<int[]> r = new ArrayList<int[]>( vSize );
		
		
		
		for( int i = 0; i < vSize; i++ ) {
			int[] vertex = orderedVertexList.get(i);
			//vertex[vDegree] *= 2; // the authors' original method had "doubled degrees"
			maxDegree = Math.max(maxDegree, vertex[vDegree]);
			r.add( vertex );
			vOrder.add(i);
		}
		
		
		// sort by degree in descending order
		Collections.sort( r, descendingDegreeComparator );
		
		// index in vertices
		int vi = vSize - 1;
		
		int rMinIndex = r.size() - 1;
		
		while( rMinIndex > 0 ) {
			// locate the vertices with minimum degree - yields a subset "Rmin" which is a subarray from rMinIndex to the end of r
			rMinIndex = r.size() - 1;
			int minDegree = r.get( rMinIndex )[vDegree];  // last one
			
			
			// calculate index of last (degree-sorted) vertex of the same degree
			while( (rMinIndex > 0) && (r.get(rMinIndex-1)[vDegree] == minDegree ) ) {
				--rMinIndex;
			}
			
			if( rMinIndex == 0 )
				break;
			
			// if Rmin has more than one element, then we sort by ex-degree
			if( rMinIndex < r.size()-1 ) {
				List<int[]> subList = r.subList(rMinIndex, r.size() );
				Collections.sort( subList, descendingExDegreeComparator );
				//boolean lol = true;
			}
			
			// vertex with minimum ex-degree goes into the ordered set of vertices (filled from back to front)
			int[] p = r.get( vi );  // back of r
			vOrder.set(vi, p[vIndex]);
			r.remove(vi);
			--vi;
			
			// decrease the degree of remaining vertices that are adjacent to p
			Collection<Integer> adjacent = adjList.get(p[vIndex]);

			int rs = r.size();
			
		
			
			// might optimise this later
			for( int i = rs; i > 0; --i ) {
				if( adjacent.contains( r.get(i-1)[vIndex] ) ) {  // is in adjacency list/matrix
					int rId = --r.get(i-1)[vDegree];  // modified degree value
					
					// sort the modified vertex immediately, so that it's at the end
					for( int j = i; (j < rs) && (rId < r.get(j)[vDegree]); ++j ) {
						Collections.swap(r, j-1, j);
					}
				}
			}
			
		}
		
		List<Integer> dummySet = new ArrayList<Integer>( r.size() );
		for( int i = 0; i < r.size(); i++ ) { dummySet.add( r.get(i)[vIndex] ); }
		//maxDegree = 24;
		
		List<int[]> colourSorted = assignColours( dummySet, 0, 0 );
		int m = colourSorted.get( colourSorted.size() - 1 )[ vColour ];  // largest identified colour
		int mMax = r.size() + maxDegree - m;
		for( int i = r.size(); i > 0; --i ) {
			vOrder.set(i-1, colourSorted.get(i-1)[vIndex]);
			storedColour.add(0, colourSorted.get(i-1)[vColour] );
			colourSorted.remove( i-1 );
		}
		
		/* XXX  This lower bound stage I translated from their C++ source code.  I don't recall seeing this lower bound thing in the paper though...
		 * 
		* if the degree of vertices in r (which all have the same degree) equals r.size() - 1, then r is a clique
		* this is important for determining the lower bound for the clique detection algorithm
		*/
        if (r.size() == r.get(0)[vDegree]+1) {
        	// fill the clique with the new vertex numbers (these are actually numbers [0..r.size()-1]), 
        	// the ones that will be set by the orderVertices function also for the rest of the graph
        	
        	Set<int[]> initClique = new LinkedHashSet<int[]>( r.size() );
        	for (int i = 0; i < r.size(); ++i) {
        		initClique.add( new int[]{ i, 0, 0, r.get(i)[vColour] } );
        		//initClique.add( new int[]{ i, 0, 0, i } );
        	}
        	
        	// filler, hack for lower bound
        	/*for (int i = r.size(); i < 37; ++i) {
        		initClique.add( new int[]{ i, 0, 0, i } );
        		//mMax++;
        		++m;
        	}*/
        	
        	if( verbose )
        		System.out.println( "initial clique size - " + initClique.size()  );
        	
        	if( deltaYPossible ) {
        		if( ! ConvenienceTools.deltaYExchangeOccured(hsMol, qMol, initClique) )
        			cliqueMax = initClique;
        	} else {
        		cliqueMax = initClique;
        	}
        }
		
		/*vOrder.clear();
		vOrder.add(4);
		vOrder.add(2);
		vOrder.add(3);
		vOrder.add(1);
		vOrder.add(0);*/
		
		// re-order vertices
		try {
			orderVertices( vOrder );
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		orderedVertexList.clear();
		// numbering - first few colors remain as they are, the rest are filled in with up to mmax
		for (int i = 0; i < vSize; ++i) {
            if (i < r.size()) {
                //assignVertexNumber(vertices, color, i, i, storedColor[i]);
                orderedVertexList.add( new int[]{ i, adjList.get(i).size(), 0, storedColour.get(i) }  );
            } else if (i < mMax) {
                ++m;
                //assignVertexNumber(vertices, color, i, i, m);
                orderedVertexList.add( new int[]{ i, adjList.get(i).size(), 0, m }  );
            } else {
                //assignVertexNumber(vertices, color, i, i, maxDegree + 1);
            	orderedVertexList.add( new int[]{ i, adjList.get(i).size(), 0, maxDegree + 1 }  );
            }
            
            if( verbose ) {
            	System.out.println( orderedVertexList.get(i)[0] + " " + orderedVertexList.get(i)[1] );
            }
            
            //std::cout << i << ": " << this->graph->mapping[color[i].first] << "," << color[i].second << "  ";
        }
		
		// debug 
        if( verbose ) {
        	for( Collection<Integer> row : adjList ) {
        		System.out.println( row.size() + " elements, " + row.toString() );
        	}
        	
        	try {
				Thread.sleep(1200);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }
		
	}
	
	
	

	
	
	
	
	private void orderVertices( List<Integer> order ) throws Exception {
		
		int n = adjList.size();
        if (order.size() != n) 
            throw new Exception( "Invalid size vector in orderVertices" );
       
        // reverse look-up of order
        List<Integer> orderInverse = new ArrayList<Integer>( order );  // dummy copy
        for( int i = 0; i < order.size(); i++ ) {
        	orderInverse.set( order.get(i), i );
        }

        adjListInverse = new ArrayList<Collection<Integer>>( n );
        List<Integer> vTranslations2 = new ArrayList<Integer>( n );
        
        // remap to temporary adjacencyMatrix
        List<Collection<Integer>> adjacencyMatrix2 = new ArrayList<Collection<Integer>>( n );
       
        for (int i = 0; i < n; ++i) {
            adjacencyMatrix2.add( new BitSetExtended<Integer>( n ) );
            vTranslations2.add(i, vertexTranslations.get( order.get(i) ) );
           
            Collection<Integer> adjRowI = adjList.get( order.get(i) );
            //for (int j = 0; j < adjRowI.size(); ++j) {
            for( Integer j : adjRowI ) {
                //adjacencyMatrix2[i][j] = adjRowI[order[j]] == true;
            	
            	adjacencyMatrix2.get(i).add( orderInverse.get( j ) );
            	
                //invAdjacencyMatrix[i][j] = (i != j) & !adjacencyMatrix2[i][j];
            }
            
            BitSetExtended<Integer> origBitSet = (BitSetExtended<Integer>) adjacencyMatrix2.get(i);
            BitSet inversed = (BitSet) origBitSet.getBitSet().clone();
            inversed.flip( 0, n );
            adjListInverse.add( new BitSetExtended<Integer>( inversed ) );
        }
        
        adjList = adjacencyMatrix2;
        vertexTranslations = vTranslations2;

        
        
		
	}
	
	
	
	
	

	
	// input settings
	//protected final int maxSteps = 2002000;
	
	// internal data for processing
	protected List<int[]> orderedVertexList;
	protected List<Collection<Integer>> adjList;
	protected List<Collection<Integer>> adjListInverse;
	protected List<Integer> vertexTranslations;
	//protected List<int[]> vertexInfo;
	
	// vertex index info
	private static final int vIndex = 0;
	private static final int vDegree = 1;
	private static final int vExDegree = 2;
	private static final int vColour = 3;
	
	// outputs
	protected Collection<int[]> cliqueMax;
	
	protected long orderingTime = 0;
	protected long colouringTime = 0;
	protected long expansionTime = 0;
	

	protected boolean performOrdering = true;
	

	

	
}
