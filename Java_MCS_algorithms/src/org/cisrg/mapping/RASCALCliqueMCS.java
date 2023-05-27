package org.cisrg.mapping;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;


import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;




/**
 * 
 * @author edmund duesbury
 * @date 4th December 2014
 */
public class RASCALCliqueMCS extends CliqueDetection {

	/*
	 * Differences between this implementation and John Raymond's paper:
	 * 
	 * - MSI isn't used - we don't care if two molecules are very different, we still want the MCES.  Lower bound is therefore maximum clique size
	 * - We don't use the colouring algorithm, don't see the need if the partitioning is better
	 * 
	 * 
	 */
	
	public RASCALCliqueMCS(ModularProductOptions opts) {
		super( opts );

	}
	
	
	public RASCALCliqueMCS(ModularProductOptions opts, boolean swapTest) {
		super( opts, swapTest );

	}
	

	
	/**
	 * 
	 * XXX  I don't like this lower bound estimate, as specifying a "Minimum Similarity Index" close to zero will yield a negative lower bound.  
	 * In addition this restricts the algorithm from comparing structurally different molecules 
	 * 
	 * @param graph1
	 * @param graph2
	 */
	private void calculateLowerBound( IAtomContainer graph1, IAtomContainer graph2 ) {
		
		int labelDifferences = 0;
		
		Map<String, Integer> aTypes = new HashMap<String,Integer>();
		Map<String, Integer> bTypes = new HashMap<String,Integer>();


		// classify hsMol atoms
		for( IAtom aAt : graph1.atoms() ) {
			
			String atProp = (String) aAt.getProperty( ConvenienceTools.atomTypeProperty );
			if( atProp == null )
				atProp = aAt.getSymbol();
			
			if( ! aTypes.containsKey( atProp ) ) {
				aTypes.put( atProp, 1 );
			} else {
				aTypes.put( atProp, aTypes.get(atProp) + 1 );
			}
		}
		
		// classify qMol atoms
		for( IAtom bAt : graph2.atoms() ) {
			
			String atProp = (String) bAt.getProperty( ConvenienceTools.atomTypeProperty );
			if( atProp == null )
				atProp = bAt.getSymbol();
			
			if( ! bTypes.containsKey( atProp ) ) {
				bTypes.put( atProp, 1 );
			} else {
				bTypes.put( atProp, bTypes.get(atProp) + 1 );
			}
		}
		
		
		// get all possible keys
		Set<String> keys = new HashSet<String>( bTypes.keySet() );
		keys.addAll( aTypes.keySet() );
		
		for( String k : keys ) {
			if( ! aTypes.containsKey(k) ) 
				aTypes.put(k, 0);
			
			if( ! bTypes.containsKey(k) ) 
				bTypes.put(k, 0);
			
			
			if( aTypes.containsKey(k) ) {
				labelDifferences += Math.abs( aTypes.get(k) - bTypes.get(k) );
			}
		}
		
		
		lowerBound = (int) Math.sqrt( 
				MinSimilarityIndex * 
				( hsMol.getAtomCount() + hsMol.getBondCount() ) * 
				( qMol.getAtomCount() + qMol.getBondCount() ) 
		) - hsMol.getAtomCount() + labelDifferences;
		
		if( lowerBound < 0 )
			lowerBound = 0;
		
	}
	
	
	
	private int labeledProjectionUpperBound( List< ArrayList<Integer> > parts, IAtomContainer graph1, IAtomContainer graph2 ) {
		
		int bound = 0;
		
		// we're projecting the modular product graph (via partition sets) back onto the original graphs 
		// Start by obtaining unique edge IDs (with Sets)
		Set<Integer> edges1 = new HashSet<Integer>( graph1.getBondCount() );
		Set<Integer> edges2 = new HashSet<Integer>( graph2.getBondCount() );
		
		
		// go through all the partition vertices
		for( List<Integer> pList : parts ) {
			for( Integer p : pList ) {
				int[] pair = modProd.getNodes().get(p);
				edges1.add( pair[0] );
				edges2.add( pair[1] );
			}
		}
		
		/*for( int[] pair : modProd.getNodes() ) {
			edges1.add( pair[0] );
			edges2.add( pair[1] );
		}*/
		
		/*for( int i = 0; i < graph1.getBondCount(); i++ )
			edges1.add(i);
		
		for( int i = 0; i < graph2.getBondCount(); i++ )
			edges2.add(i);*/
		
		
		
		// count normal bond types.  Store any abnormal bond types
		// in hs, any SMARTS bonds remaining - test with any "excesses" from db mol and count.  If there's a match, remove from stored list
		// ignore special atom types for now
		
		// so yeah - loop thru DB molecule (as it's smaller) and compare with HS bonds.  If match, delete both bonds from lists and add 1 to thing
		List<Integer> edges2List = new ArrayList<Integer>(Arrays.asList( ( edges2.toArray( new Integer[0] ) ) ));

		dbl: for( int d = 0; d < edges2List.size(); d++ ) {
			Integer db = edges2List.get(d);
			IBond dbBond = graph2.getBond(db);
			
			for( Integer hs : edges1 ) {
				IBond hsBond = graph1.getBond(hs);
				
				if( GenerateCompatibilityGraphEdges.nodesMatch(graph1, hsBond, graph2, dbBond, true, null) ) {
					bound++;
					edges2List.remove( d );
					edges1.remove( hs );
					d = -1;  // gets reset back to 0 on the continue statement
					continue dbl;
				}
				
			}
		}
		
		return bound;
	}
	
	
	private int colourUpperBound( List< ArrayList<Integer> > parts ) {
		
		int colour = 1;
		
		Map<Integer, ArrayList<Integer>> adjSet = new HashMap<Integer, ArrayList<Integer>>();
		
		// Ordered sets instead of Lists to store actual vertices with constant time access 
		List<int[]> vertices = new ArrayList<int[]>();
		LinkedHashSet<Integer> vertices1 = new LinkedHashSet<Integer>(  );
		LinkedHashSet<Integer> verticesCopy = new LinkedHashSet<Integer>(  );
		
		// conversion of partitions to sets
		for( ArrayList<Integer> part : parts ) {
			for( Integer n : part ) {
				vertices1.add(n);
			}
		}
		
		// make adjacency list (or whatever) relevant to the partition vertices only
		for( int n : vertices1 ) {
			// copy to avoid modifying the original modular product
			
			ArrayList<Integer> adj = new ArrayList<Integer>( modProd.getAdjacencyList().get(n) )  ;
			
			adj.retainAll( vertices1 );  // perform intersection
			
			adjSet.put(n, adj); 
			
			vertices.add( new int[]{ n, adj.size() } );
			//vertexDegrees.put( n, adj.size() );
		}

		Collections.sort(vertices , descendingDegreeComparator);
		vertices1.clear();
		
		for( int[] p : vertices ) {
			vertices1.add( p[0] );
			verticesCopy.add( p[0] );
		}
		
		//System.out.println( "Colouring, LinkedHashSet creation TIME - " + (System.currentTimeMillis() - time1) );
		
		
		
		
		while( ! vertices1.isEmpty() ) {
			
			while( ! verticesCopy.isEmpty() ) {
				
				// get first element that's not null in the copy
				Integer v = null;
				v = verticesCopy.iterator().next();
				
				//int index = vertices.indexOf(v);
				
				//System.out.println( vertices.get(0) + " " + index  + " " + vertices1.contains(v) );
				vertices1.remove( v );
				
				//System.out.println( "v " + v + " " + v[0] + "," + v[1] );
				//System.out.println( "neighbours " + Arrays.asList(adjSet.get(index)) + " " + adjSet.get(index).length );
				
				//System.out.println( "neighbours2 " + orderedAdjListMap.get(v) + " " + orderedAdjListMap.get(v).length );
				

				
				// delete neighbours of v from copy, as well as v itself
				for( Integer neighbour : adjSet.get( v ) ) {
					//System.out.print( neighbour[0] + " " + neighbour[1] + ", " );
					verticesCopy.remove( neighbour );  // remove the object at the index, hence the "cast"
					//System.out.println( "size2 = " + vertices1.size() + " " + verticesCopy.size() + " | colour = " + colour  );
				}
				//System.out.println("");
				verticesCopy.remove(v);
				

				//System.out.println( "size = " + vertices1.size() + " " + verticesCopy.size() + " | colour = " + colour  );
				//verticesCopy.remove( v );  // delete the first element as well as its neighbours
				
				//System.out.println( verticesCopy.size() );
				
				//if( colour > maxSize ) {
				//	v[5] = colour;
				//}
				
			}
			
			colour++;
			
			/*
			for( int e = 0 ; e < orderedVertices.size(); e++ ) {
				int[] ea = orderedVertices.get(e);
				System.out.println( e + " [ " + ea[0] + "," + ea[1] + "," + ea[2] + "," + ea[3] + "," + ea[4] + "]"  );
			}
			*/
			
			
			
			// vcopy = v
			verticesCopy.clear();
			Iterator<Integer> it = vertices1.iterator() ;
			while( it.hasNext() ) {
				verticesCopy.add( it.next() );
			}
			
		}
		
		
		
		return colour;
		
	}
	
	
	private int determineUpperBound( GenerateCompatibilityGraphEdges modProd, List< ArrayList<Integer> > parts ) {
		
		
		
		//List<ArrayList<Integer>> initialPartitions = initialPartitioning( modProd );
		Collections.sort( parts, sizeDescendingComparator );
		
		// now to re-partition nodes (thus shrinking them and obtaining an upper bound)
		int upperBound1 = refinePartitions(parts);

		Collections.sort( parts, sizeDescendingComparator );
		
		int upperBound2 = labeledProjectionUpperBound(parts, hsMol, qMol);
		
		int upperBound3 = upperBound2;  // represents colouring upper bound, if the number of nodes in the partition set is small enough
		
		int psNodes = 0;
		for( ArrayList<Integer> part : parts ) {
			psNodes += part.size();
		}
		
		/*
		 * XXX  I'm keeping the colouring algorithm, though Dave Cosgrove and I both find that this generally is still outperformed
		 * by the current partition-refine method for "colouring."
		 */
		if( psNodes >= 150 ) {
			upperBound3 = colourUpperBound(parts) ;
		}
		
		// debug section
		/*int ps = parts.size();
		for( int p = 0; p < ps; p++ ) {
			
			for( Integer n : parts.get(p) ) {
				int[] node = modProd.getNodes().get( n );
				System.out.print( "[" + node[0] + "," + node[1] + "]"  );
			}
			
			System.out.println();
		}*/
		
		//System.out.println( "upper bounds: " + upperBound1 + ", " + upperBound2 + ", " + upperBound3 );
				
		int upperBound12 = Math.min( upperBound1, upperBound2 );
		
		return Math.min( upperBound12, upperBound3 );
		
	}
	
	
	/**
	 * XXX As we don't know how the labelled projection "partitioning" actually works, I've implemented Dave Cosgrove's method.
	 * 
	 * In this case, partitions are created based on the fact that a node pair (or whatever) in the modular product cannot be connected to another node pair
	 * where one of the elements is contained in the given node pair.  Thus, they are guaranteed partitions.  
	 * 
	 * Should yield a number of partitions equal to min(edges1, edges2)
	 * 
	 * NOTE - assumes modular product nodes are sorted in ascending order of first node
	 * 
	 * @param modProd
	 * @return
	 */
	private List<ArrayList<Integer>> initialPartitioning() {
		
		ArrayList< ArrayList<Integer> > parts = new ArrayList< ArrayList<Integer> >( 50 );
		
		int prevIndex = -1;
		ArrayList<Integer> partition = null;
		List<int[]> nodes = modProd.getNodes();
		
		// assumes modular product nodes are sorted in ascending order of first node (atom index of first compound)
		
		for( int n = 0, ns = nodes.size(); n < ns; n++ ) {
			int[] node = nodes.get(n);
			
			// new atom index on first molecule - new partition (sorting assumption)
			if( node[0] > prevIndex ) {
				prevIndex = node[0];
				
				partition = new ArrayList<Integer>(50);
				parts.add(partition);
				//System.out.println();
			}
			
			partition.add( n );
			//System.out.print( "[" + node[0] + "," + node[1] + "]"  );
		}
		
		return parts;
	}
	
	
	
	/**
	 * 
	 * @param modProd
	 * @param parts
	 * @param mBound  lower bound of clique size (I imagine this'll be the size of the current clique)
	 */
	@Deprecated
	public int refinePartitions( GenerateCompatibilityGraphEdges modProd, List< ArrayList<Integer> > parts, int mBound ) {
		
		// no point in re-shuffling if there're no pEx partitions (ones greater than the lower bound)
		//if( mBound >= parts.size() )
		//	return parts.size();
		
		int ps = parts.size();
		System.out.println("pre-refinement, level " + ps );
		//for( int pex = ps -1; pex > 0; pex-- ) {
		for( int pex = ps - 1; pex >= mBound; pex++ ) {
			
			pexv: for( int vInd = 0; vInd < parts.get(pex).size(); vInd++ ) {
				Integer v = parts.get(pex).get(vInd);
				
				for( int pm = 0; pm < pex; pm++ ) {
				
					//List<Integer> vNeighbours = modProd.getAdjacencyList().get(v);
					int size = modProd.intersectNodeNeighbours(v, parts.get(pm) ).size();
					
					// if no neighbours intersect with pM partition, reshuffle the vertex into the pM partition
					if( size == 0 ) {
						parts.get(pm).add(v);
						parts.get(pex).remove(v);
						
						if( parts.get(pex).size() == 0 ) {
							parts.remove(pex);
							pex--;
							ps--;
							//break outer;
						}
						
						continue pexv;
					}
					
				}
			}
			
		}
		
		System.out.println("post-refinement, level " + parts.size() );
		
		
		
	
		return parts.size();
	}
	
	
	
	/**
	 * Usually called when partitions have been modified in some way (e.g. a vertex has been removed from the partition set)
	 * 
	 * Re-assigns any vertices in each partition to existing partitions, if they are found to have no neighbours in the given partition to assign to.
	 * This assumes that the partitions have been sorted in descending order of size
	 * 
	 * 
	 * 
	 * @param parts
	 * @return	size (number of partitions) of new partitioning
	 */
	public int refinePartitions( List< ArrayList<Integer> > parts ) {
		
		// no point in re-shuffling if there're no pEx partitions (ones greater than the lower bound)
		//if( mBound >= parts.size() )
		//	return parts.size();
		
		//System.out.println("pre-refinement, level " + parts.size() );
		
		if( parts.isEmpty() )
			return 0;
		
		int ps = parts.size();
		for( int pex = ps -1; pex >= 0; pex-- ) {  // start from smallest partition
			
			pexv: for( int vInd = 0; vInd < parts.get(pex).size(); vInd++ ) {
				Integer v = parts.get(pex).get(vInd);
				
				for( int pm = 0; pm < pex; pm++ ) {  // start from largest partition - want to assign to biggest first
				
					//List<Integer> vNeighbours = modProd.getAdjacencyList().get(v);
					int size = modProd.intersectNodeNeighbours(v, parts.get(pm) ).size();
					
					// if no neighbours intersect with pM partition, reshuffle the vertex into the pM partition
					if( size == 0 ) {
						parts.get(pm).add(v);
						parts.get(pex).remove(v);
						
						if( parts.get(pex).size() == 0 ) {
							parts.remove(pex);
							pex--;
							ps--;
							//break outer;
						}
						continue pexv;
					}
					
				}
			}
			
		}
		
		//System.out.println("post-refinement, level " + parts.size() );
		
		
		
	
	
		return parts.size();
	}



	private void sortPartitionsOnDegree( List< ArrayList<Integer> > parts ) {
		
		for( int pIndex = 0, ps = parts.size(); pIndex < ps; pIndex++ ) {
			
			ArrayList<Integer> part = parts.get(pIndex);
			ArrayList<int[]> newOrder = new ArrayList<int[]>( part.size() );
			ArrayList<Integer> newPart = new ArrayList<Integer>( part.size() );
			
			// pairs - to contain indices
			for( Integer p : part ) {
				int[] pair = new int[]{ p, modProd.getAdjacencyList().get(p).size() };
				newOrder.add( pair );
			}
			
			// perform degree sorting
			Collections.sort( newOrder, descendingDegreeComparator );
			//Collections.reverse(newOrder);
			
			for( int[] pair : newOrder ) {
				newPart.add( pair[0] );
			}
			
			parts.set(pIndex, newPart);
		}
		
	}
	
	
	
	/**
	 * keep partitions & vertices in partitions which are neighbours of vertex "v"
	 * 
	 * @param parts  current partitioning of vertices
	 * @param size  number of partitions in parts
	 * @param v  the (index of the) vertex to find neighbours of in the modular product graph
	 * @return 
	 */
	private List< ArrayList<Integer> > vertexNeighbourhoodPartitions( List< ArrayList<Integer> > parts, int size, int v ) {
		
		List< ArrayList<Integer> > newParts = new ArrayList< ArrayList<Integer> >( size );
		
		// explore neighbours of the vertex (in all partitions) to continue the search 
		for( int pIndex = 0; pIndex < size; pIndex++ ) {
			//ArrayList<Integer> nPart = new ArrayList<Integer>( parts.get(pIndex) );
			
			ArrayList<Integer> nPart = new ArrayList<Integer>( modProd.intersectNodeNeighbours( v , parts.get(pIndex) ) );
			
			if( nPart.size() > 0 /*&& ! newParts.contains( nPart )*/ ) {
				//boolean unique = true;
				
				
				// test whether the "new partition set" already contains this vertex
				// this part might not be needed
				/*for( ArrayList<Integer> enPart : newParts ) {
					for( Integer nPartInd : nPart ) {
						if( enPart.contains( nPartInd ) ) {
							unique = false;
							System.out.println( "unique is false at " + parts.size() );
							break;
						}
					}
				}*/
				
				// test whether the "new partition set" already contains this vertex
				// this part might not be needed
				/*for( Integer nPartInd : nPart ) {
					if( partitionContainsVertex( newParts, nPartInd ) ) {
						unique = false;
						//System.out.println( "unique is false at " + parts.size() );
						break;
					}
				}*/
					
					
				//if( unique ) {
					newParts.add(nPart);
				//}
			}
		}
		
		return newParts;
	}
	
	
	/*
	 * Kikusts pruning:
	 * 
	 * On part where we start removing nodes from the partitions, 
	 * test each vertex in pEx
	 * 
	 * condition 1:
	 * find a vertex v_j that's a neighbour of the neighbours of v_i (the clique vertex) but NOT in the neighbourhood
	 * if it's next to another vertex not in the neighbours of v_i, keep it
	 * 
	 * condition 2:
	 * v_j is a neighbour of v_i.  
	 * Find two neighbours to v_j not in the neighbourhood of v_i, which are adjacent to each other.
	 * Keep v_j if this' the case
	 * 
	 * 
	 * 
	 * So, get all vertices that're NOT (1) and AND (2) of parts and neighbours(v_i) (newParts), removing clique vertex from both sets
	 * 
	 * (1) - for each vertex v_j in (1) with a neighbour in (1, being v_m) & a neighbour in (2), intersect neighbours of v_j or v_m with newParts.  If size of result > 0 then keep v_j in parts
	 * (2) - for each vertex v_j, intersect neighbours of v_j with (1).  From resulting set, generate all vertex pairs and look for adjacent pairs.  
	 * if there is one, keep v_j 
	 * 
	 */
	private void kikustsPruning( List< ArrayList<Integer> > parts, List< ArrayList<Integer> > newParts, int lb ) {
		
		int ps = parts.size();

		// Kikusts pruning
		HashSet<Integer> partsEx = new HashSet<Integer>( modProd.getNodes().size() );
		for( int pIndex = lb; pIndex < ps; pIndex++ ) {
			partsEx.addAll( parts.get(pIndex) );
		}

		HashSet<Integer> nvPartsEx = new HashSet<Integer>( ps );
		for( int pIndex = 0; pIndex < newParts.size(); pIndex++ ) {
			nvPartsEx.addAll( newParts.get(pIndex) );
		}

		HashSet<Integer> nonvNeighboursEx = new HashSet<Integer>( partsEx );
		nonvNeighboursEx.removeAll(nvPartsEx);  // get non-neighbours of clique vertex in the partitions




		// condition 1 - a vertex must be adjacent to 
		for( Integer vj : nonvNeighboursEx ) {
			
			HashSet<Integer> neighboursVj = new HashSet<Integer>( modProd.getAdjacencyList().get( vj ) );  // neighbours of vj
			neighboursVj.retainAll(partsEx);  // keep in domain of the partitions being explored
			neighboursVj.removeAll(nvPartsEx);  // keep those that aren't neighbours of the clique

			boolean adjPairNonCliqueNeighbour = false;
			
			for( Integer vm : nonvNeighboursEx ) {
				if( neighboursVj.contains(vm) ) {
					adjPairNonCliqueNeighbour = true;
					break;
				}
			}
			
			/*boolean cliqueNeighboursSatisfied = false;
			boolean nonCliqueNeighboursSatisfied = false;

			for( Integer vTest : neighboursVj ) {
				if( nvPartsEx.contains(vTest) )
					cliqueNeighboursSatisfied = true;  

				if( nonvNeighboursEx.contains(vTest) )
					nonCliqueNeighboursSatisfied = true;
			}

			// v_j must be connected to the neighbours of clique vertex, and connected to a vertex that isn't in the neighbourhood
			if( ! cliqueNeighboursSatisfied || ! nonCliqueNeighboursSatisfied ) {
				partsEx.remove(vj);
				System.out.println("Kikusts condition 1 deletion");
			}*/
			
			if( ! adjPairNonCliqueNeighbour ) {
				partsEx.remove(vj);
				//System.out.println("Kikusts condition 1 deletion");
			}
		}

		// condition 2
		for( Integer vj : nvPartsEx ) {
			HashSet<Integer> nonCliqueNeighboursVj = new HashSet<Integer>( modProd.getAdjacencyList().get( vj ) );
			nonCliqueNeighboursVj.retainAll(nonvNeighboursEx);  // only keep neighbours that aren't in the clique-vertex neighbourhood

			boolean adjacentPair = false;

			// now we have neighbours of v_j that're not in the clique-vertex neighbourhood, we keep the vertex if there're any adjacent pairs
			outer:  for( Integer vn : nonCliqueNeighboursVj ) {
				for( Integer vm : nonCliqueNeighboursVj ) {
					if( modProd.getAdjacencyList().get(vn).contains(vm) ) {  // test for a pair of vertices
						adjacentPair = true;
						break outer;
					}
				}
			}

			if( ! adjacentPair ) {
				partsEx.remove(vj);
				//System.out.println("Kikusts condition 2 deletion");
			}
		}



		// translating removals to actual partitions
		for( int pIndex = lb; pIndex < ps; pIndex++ ) {
			List<Integer> partition = parts.get(pIndex);
			for( int p = partition.size() - 1; p >= 0; p-- ) {
				if( ! partsEx.contains( partition.get(p) ) ) {
					partition.remove(p);
					//System.out.println("Kikusts - partition vertex removed");
					if( partition.isEmpty() ) {
						parts.remove(pIndex);
						ps--;

					}
				}
			}
		}
		
	}
	
	
	
	
	
	private boolean partitionContainsVertex( List< ArrayList<Integer> > parts, Integer vertex ) {
		
		for( ArrayList<Integer> part : parts ) {
			if( part.contains(vertex) )
				return true;
		}
		
		return false;
	}
	
	
	
	
	public void branch( List< ArrayList<Integer> > parts, List<Integer> currentClique ) throws CDKException {
		
		++numberOfSteps;
		
		if( numberOfSteps % 10000 == 0 && System.currentTimeMillis() > ( mcsStartTime + expansionTimeLimitMs ) )
			throw new CDKException("RASCAL time limit exceeded for clique search");
	
		int ps = parts.size();
		
		//lowerBound = Math.max( 53, lowerBound - 1 );
		
		if( currentClique.size() + ps <= lowerBound )
			return;  // can't be better than what we have currently
		
		
		
		
		// deal with upper bound first (step 2 of RASCAL in paper)
		// also sorts partitions and performs repartitioning
		int upperBound = determineUpperBound( modProd, parts );
		ps = parts.size();
		
		//int upperBound2 = refinePartitions(modProd, parts, lowerBound);
		
		// account for lower bound
		/*if( mBound < currentClique.size() )
			mBound = currentClique.size();*/
		
		
		
		// start of "Step 3" from here - branching and clique expansion (looping)
		
		while( ! parts.isEmpty() && (upperBound + currentClique.size() >= lowerBound) ) {
			
			//ps = parts.size();
			
			
			
			// choose a vertex from the final partition (should be in pEx)
			// the 'pEx' condition is the second part of the while loop condition
			List<Integer> lastPartition = parts.get( ps - 1 );  // choose last vertex in last partition
			Integer v = lastPartition.get( lastPartition.size() - 1 ) ;
			
			currentClique.add( v ); // add this vertex to clique
			
			
			//System.out.println("pre-neighbours, level " + parts.size() );
			
			// keep partitions & vertices in partitions which are neighbours of vertex "v"
			List< ArrayList<Integer> > newParts = vertexNeighbourhoodPartitions( parts, ps, v );
			
			//System.out.println("post-neighbours, level " + newParts.size() + ", before was " + parts.size() );
			
			// sorting and re-shuffling (is this needed?)
			//Collections.sort( newParts, sizeDescendingComparator );
			//refinePartitions(newParts);
			
			//System.out.println("post-refinement, level " + newParts.size() + " at lower bound " + lowerBound + ", upper bound is " + upperBound + ", clique size is " + currentClique.size() );
			
			
			// go to next depth of search tree (using modified partitions and boundaries)
			if( newParts.isEmpty() ) {
				saveClique( currentClique );
				//System.out.println("clique size " + currentClique.size() );
			} else {
				//System.out.println("pre-recursion, level " + parts.size() );
				branch( newParts, currentClique );
				//System.out.println("post-recursion, level " + parts.size() );
			}
			
			
	
			
			// removal of said vertex and updating lower bound
			lastPartition.remove( lastPartition.size() - 1 );
			if( lastPartition.isEmpty() ) {
					parts.remove( ps - 1 );
					ps--;
					//mBound--;
			}
			
			
			
			currentClique.remove( currentClique.size() - 1 );
			
			
			/*
			 * XXX I didn't implement the equivalence class pruning due to time constraints
			 * 
			 * In addition, probably due to me not thoroughly checking the Kikusts pruning part - turning it off actually increases its speed.
			 */
			
			/*// perform Kikusts pruning at backtrack
			//if( currentClique.size() <= 1 ) {
				kikustsPruning(parts, newParts, currentClique.size() );
				ps = parts.size();
			//}
			 */			
			
			
					
			
			//lowerBound = currentClique.size();  // important to decrease lower bound
			
			// do other backtracking stuff (pEq and KP)
		
		}
		
	}
	
	
	
	
	
	@Override
	protected void findCliques() {
		//calculateLowerBound(hsMol, qMol) ;
		//lowerBound = 34;
		//System.out.println("lowerBound " + lowerBound);
		List<ArrayList<Integer>> initialPartitions = initialPartitioning( );
		Collections.sort( initialPartitions, sizeDescendingComparator );
		
		if( modProd.getNodes().size() >= 400 ) {
			refinePartitions(initialPartitions);
			sortPartitionsOnDegree( initialPartitions );
		}
		
		try {
			//branchStartTime = System.currentTimeMillis();
			branch(initialPartitions, new LinkedList<Integer>() );
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	// descending order of cardinality
	protected Comparator<List<?>> sizeDescendingComparator = new Comparator<List<?>>() {
			 public int compare(List<?> o1, List<?> o2) {
				 
		         return -1 * Integer.compare( o1.size(), o2.size() );   
				
			 } 
		};
		
		// descending order of degree
		protected Comparator<int[]> descendingDegreeComparator = new Comparator<int[]>() {
				 public int compare(int[] o1, int[] o2) {

			            int dComp;
			            // sort in descending order of degree first
			            if( o1[1] == o2[1] ) {
			            	dComp = 0;
			            } else if( o1[1] > o2[1] ) {
			            	dComp = -1;
			            	//Collections.swap( sortedVertices2, o1, o2 );
			            } else {
			            	dComp = 1;
			            	//Collections.swap( sortedVertices2, o2, o1 );
			            }
			            
			            return dComp;
				 } 
			};
			
			
			
	
	/*
	 *  XXX  This isn't used to allow comparisons between structurally different molecules
	 *  Really, a better lower bound is needed (much like the one used by Depolli et al)
	 */
	protected double MinSimilarityIndex = 0.5;






	
	
	
}
