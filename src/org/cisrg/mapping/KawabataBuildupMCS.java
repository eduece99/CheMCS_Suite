package org.cisrg.mapping;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import org.openscience.cdk.isomorphism.matchers.smarts.AnyAtom;




/**
 * 
 * 
 * (Known) Differences between this implementation and the paper for the Kawabata algorithm:
 * 
 * - MCES, not MCIS is being found.  As such, a list of bond pairs, NOT atom pairs is being analysed.  In fact, no atom pairs really exist in this implementation
 * - BondEquivalence and BondConnection functions obviously look for whether bonds are adjacent or not - using a common vertex
 * - no method of calculating the shortest path is mentioned in the paper
 * - Heuristic Score dNeighbour is different
 * 
 * @author edmund
 *
 */




public class KawabataBuildupMCS extends MCSMethods {
	
	
	
	public KawabataBuildupMCS(int topoDistanceLimit, boolean isConnectedMCS, boolean useRaymondHeuristics) {
		
		if( topoDistanceLimit >= 0 ) {
		
			maxTopologicalDifference = topoDistanceLimit;
			useTopoDistance = true;
			
		}
		
		connectedMCS = isConnectedMCS;
		this.useRaymondHeuristics = useRaymondHeuristics;
	}
	
	
	class MatchList {
		
		public MatchList() {
			pairs = new ArrayList<int[]>(1);
		}
		
		public MatchList( int[] initial ) {
			pairs = new ArrayList<int[]>(initSize);
			//elems1 = new HashSet<Integer>(initSize);
			//elems2 = new HashSet<Integer>(initSize);
			atomTypeCounts = new HashMap<Integer,Integer>();
			addPair( initial );
		}
		
		public MatchList( MatchList other ) {
			pairs = new ArrayList<int[]>( other.getPairs() );
			score = other.getScore();
			atomTypeCounts = new HashMap<Integer,Integer>(other.getAtomTypeCounts());
			//elems1 = new HashSet<Integer>( other.elems1 );
			//elems2 = new HashSet<Integer>( other.elems2 );
		}
		
		
		/**
		 * total up the individual heuristic scores of the pairs
		 */
		public void calculateScore() {
			score = 0;
			scoreNeighbour = 0;
			scoreEC = 0;
			scoreTopology = 0;
			
			if( useTopoDistance ) {
				for( int i = 0; i < pairs.size(); i++ ) {
					int[] iPair = pairs.get(i);
					
					for( int j = i + 1; j < pairs.size(); j++ ) {
						int[] jPair = pairs.get(j);
						
						int iDistance =  pathDistancesHsMol[ iPair[0] ][ jPair[0] ];
						int jDistance =  pathDistancesQMol[ iPair[1] ][ jPair[1] ];
						//int topoDistance = Math.abs( pathDistancesHsMol[ iPair[0] ][ jPair[0] ] - pathDistancesQMol[ iPair[1] ][ jPair[1] ] );
						
						if( Math.min( iDistance, jDistance ) <= topoHeuristicConstant ) {
							scoreTopology += Math.abs( iDistance - jDistance );
						}
					}
				}
			}
			
			for( int[] p : pairs ) {
				scoreNeighbour += p[2];
				scoreEC += p[3];
				//score += p[5];
			}
			
			score = scoreNeighbour + scoreEC + scoreTopology;
		}
		
		
		public void addPair( int[] np ) {
			
			// hash set additions
			//elems1.add(np[0]);
			//elems2.add(np[1]);
			
			// update scores
			scoreNeighbour += np[2];
			scoreEC += np[3];
			
			if( useTopoDistance ) {
				for( int n = 0; n < pairs.size(); n++ ) {
					int[] otherPair = pairs.get(n);
					
						
						int iDistance =  pathDistancesHsMol[ otherPair[0] ][ np[0] ];
						int jDistance =  pathDistancesQMol[ otherPair[1] ][ np[1] ];
						//int topoDistance = Math.abs( pathDistancesHsMol[ iPair[0] ][ jPair[0] ] - pathDistancesQMol[ iPair[1] ][ jPair[1] ] );
						
						if( Math.min( iDistance, jDistance ) <= topoHeuristicConstant ) {
							scoreTopology += Math.abs( iDistance - jDistance );
						}
				}
			}
			
			score = scoreNeighbour + scoreEC + scoreTopology;
			
			
			// register atom types & frequency in the MatchList to test
			//for( int[] pair : pairs ) {
				
				for( IAtom aAt : hsMol.getBond( np[0] ).atoms() ) {
					if( ! atomTypeCounts.containsKey( aAt.getAtomicNumber() ) ) {
						atomTypeCounts.put( aAt.getAtomicNumber(), 1 );
					} else {
						atomTypeCounts.put( aAt.getAtomicNumber(), atomTypeCounts.get(aAt.getAtomicNumber()) + 1 );
					}
				}
			//}
			
			pairs.add(np);
		}
		
		
		
		public int getScore() {
			return score;
		}
		
		public int getScoreNeighbour() {
			return scoreNeighbour;
		}
		
		public int getScoreEC() {
			return scoreEC;
		}
		
		public int getScoreTopology() {
			return scoreTopology;
		}
		
		public List<int[]> getPairs() {
			return pairs;
		}
		
		public Map<Integer, Integer> getAtomTypeCounts() {
			return atomTypeCounts;
		}
		
		@Override
		public int hashCode() {
			return super.hashCode();
		}
		
		
		
		
		/**
		 * Compares objects.  If the other object is a MatchList then it compares them based on what it contains - two MatchLists
		 * are equal if they contain the same pairs (regardless of order)
		 */
		public boolean equals( Object other ) {
			if( other instanceof MatchList ) {
				MatchList o = (MatchList) other;
				
				if( getPairs().size() != o.getPairs().size() ) 
					return false;
				
				for( int[] p : getPairs() ) {
					if( ! o.getPairs().contains( p ) )
						return false;
				}
				
				return true;
			}
			
			return super.equals( other );
		}
		
		
		/*private boolean contains( int which, int n ) {
			
			if( which == 0 ) {
				return elems1.contains(n);
			}
			
			return elems2.contains(n);
			
		}*/
		
		private boolean contains( int which, int n ) {
			
			for( int[] pair : pairs ) {
				if( pair[which] == n ) {
					return true;
				}
			}
			
			return false;
		}
		
		protected ArrayList<int[]> pairs;
		//protected Set<Integer> elems1, elems2;
		Map<Integer, Integer> atomTypeCounts;
		
		
		
		protected int score;
		protected int scoreNeighbour = 0;
		protected int scoreEC = 0;
		protected int scoreTopology = 0;
		
		protected static final int topoHeuristicConstant = 4;
		protected static final int initSize = 50;
	}
	
	
	// ascending order of score ideally (lower is better)
	protected Comparator<MatchList> matchListComparator = new Comparator<MatchList>() {
		 public int compare(MatchList o1, MatchList o2) {

	            /*int dComp;
	            // sort in descending order of degree first
	            if( o1.getScore() == o2.getScore() ) {
	            	dComp = 0;
	            } else if( o1.getScore() > o2.getScore() ) {
	            	dComp = 1;
	            	//Collections.swap( sortedVertices2, o1, o2 );
	            } else {
	            	dComp = -1;
	            	//Collections.swap( sortedVertices2, o2, o1 );
	            }
	            
	            return dComp;*/
			 
			 if( o1.getScore() == o2.getScore() )
				 return 0;
			 
			 return o1.getScore() > o2.getScore() ? 1 : -1;
		 } 
	};
	
	
	protected Comparator<int[]> pairComparator = new Comparator<int[]>() {
		 public int compare(int[] o1, int[] o2) {

	            int dComp;
	            // sort in descending order of degree first
	            if( o1[5] == o2[5] ) {
	            	dComp = 0;
	            } else if( o1[5] > o2[5] ) {
	            	dComp = 1;
	            	//Collections.swap( sortedVertices2, o1, o2 );
	            } else {
	            	dComp = -1;
	            	//Collections.swap( sortedVertices2, o2, o1 );
	            }
	            
	            return dComp;
		 } 
	};
	
	
	/**
	 * XXX change:
	 * "bond" equivalence in this case actually refers to the presence of a common vertex (of identical label) between two bonds
	 * 
	 * Original paper tested for a common bond between 2 atoms.
	 * 
	 * @param matchList
	 * @param x
	 * @param y
	 * @return
	 */
	private boolean bondEquivalence( List<int[]> matchList, int x, int y ) {
		
		IBond bondA1 = hsMol.getBond(x);
		IBond bondA2 = qMol.getBond(y);
		
		for( int[] pair : matchList ) {
			
			IBond bondB1 = hsMol.getBond( pair[0] );
			IBond bondB2 = qMol.getBond( pair[1] );
			
			IAtom commonVertex1 = null;  // null = bonds aren't connected
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

	    	//if( commonVertex1 != null && commonVertex2 != null )
	    		//return true;
	    	
			if( (commonVertex1 == null) ^ (commonVertex2 == null) )
				return false;
		}
		
		return true;
	}
	
	
	/**
	 * XXX
	 * Unsurprisingly this had to be modified as we're testing for bond adjacency.  As with edge creation in the modular
	 * product, we seek common vertices between two bonds.
	 * 
	 * @param matchList
	 * @param x
	 * @param y
	 * @return
	 */
	private boolean bondConnection( List<int[]> matchList, int x, int y ) {
		
		IBond bondA1 = hsMol.getBond(x);
		IBond bondA2 = qMol.getBond(y);
		
		for( int[] pair : matchList ) {
			
			IBond bondB1 = hsMol.getBond( pair[0] );
			IBond bondB2 = qMol.getBond( pair[1] );
			
			IAtom commonVertex1 = null;  // null = bonds aren't connected
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

	    	if( commonVertex1 != null && commonVertex2 != null )
	    		return true;
	    	
		}
		
		return false;
	}
	
	
	private boolean differenceOfTopologicalDistance( List<int[]> matchList, int x, int y ) {
		
		for( int[] pair : matchList ) {

	    	int topoDifference = Math.abs( pathDistancesHsMol[ x ][ pair[0] ] - pathDistancesQMol[ y ][ pair[1] ] );
	    	
	    	if( topoDifference > maxTopologicalDifference )
	    		return false;
	    	
		}
		
		return true;
	}
	
	
	private void defineAtomTypes( IAtomContainer mol ) {
		
		for( IAtom at : mol.atoms() ) {
			
			String atomType = "" + at.getAtomicNumber();
			
			if( at.getFlag( CDKConstants.ISINRING ) ) {
				atomType += ".r";
			} else if( at.getAtomicNumber() != 6 ) {
				atomType += "." + mol.getConnectedBondsCount(at);
			}
			
			at.setProperty(atomTypeProperty, atomType);
		}
		
	}
	
	/**
	 * Difference in the numbers of each atom type neighbouring a bond, for a bond-pair
	 * Any-atoms are subtracted from this difference (will implement other SMARTS atom types later)
	 * 
	 * 
	 * @param a  hyperstructure bond
	 * @param b	 query molecule bond
	 * @return
	 */
	@Deprecated
	private int dNeighbour( IBond a, IBond b ) {
		
		int score = 0;
		
		Map<String, Integer> aTypes = new HashMap<String,Integer>();
		Map<String, Integer> bTypes = new HashMap<String,Integer>();
		
		// unique set of neighbouring atoms of a bond
		HashSet<IAtom> aNeighbours = new HashSet<IAtom>( hsMol.getConnectedAtomsList( a.getAtom(0) ) );
		aNeighbours.addAll( hsMol.getConnectedAtomsList( a.getAtom(1) ) );
		
		HashSet<IAtom> bNeighbours = new HashSet<IAtom>( qMol.getConnectedAtomsList( b.getAtom(0) ) );
		bNeighbours.addAll( qMol.getConnectedAtomsList( b.getAtom(1) ) );
		
		
		// classify hsMol atoms
		for( IAtom aAt : aNeighbours ) {
			
			if( aAt instanceof AnyAtom ) {
				if( ! aTypes.containsKey(anyAtomID) ) {
					aTypes.put( anyAtomID, 1 );
				} else {
					aTypes.put( anyAtomID, aTypes.get( anyAtomID ) + 1 );
				}
			}
			
			String atProp = (String) aAt.getProperty(atomTypeProperty);
			if( ! aTypes.containsKey( atProp ) ) {
				aTypes.put( atProp, 1 );
			} else {
				aTypes.put( atProp, aTypes.get(atProp) + 1 );
			}
		}
		
		// classify qMol atoms
		for( IAtom bAt : bNeighbours ) {
					
					if( bAt instanceof AnyAtom ) {
						if( ! bTypes.containsKey(anyAtomID) ) {
							bTypes.put( anyAtomID, 1 );
						} else {
							bTypes.put( anyAtomID, bTypes.get( anyAtomID ) + 1 );
						}
					}
					
					String atProp = (String) bAt.getProperty(atomTypeProperty);
					if( ! bTypes.containsKey( atProp ) ) {
						bTypes.put( atProp, 1 );
					} else {
						bTypes.put( atProp, bTypes.get( atProp ) + 1 );
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
				score += Math.abs( aTypes.get(k) - bTypes.get(k) );
			}
		}
		
		// account for any atoms - subtract number of any atoms from incompatibilities
		if( aTypes.containsKey( anyAtomID) )
			score -= aTypes.get( anyAtomID );
		
		
		// negative scores are not allowed
		return Math.max(score, 0);
		
	}
	
	
	/**
	 * XXX a change:
	 * the original "dNeighbour" function has been re-defined in the context of degenerate bond matching.
	 * 
	 * For each of the 2 bonds:
	 * - look at each neighbouring bond.  Count the number of times a bond fails to match at least 1 bond in the other bond's neighbourhood
	 * 
	 * @param aIndex
	 * @param bIndex
	 * @return
	 */
	private int dNeighbour2( int aIndex, int bIndex ) {
		
		int score = 0;
		
		
		// get neighbouring bonds
		List<Integer> aNeighbourIndices = hsMolAdjList.get(aIndex);
		List<Integer> bNeighbourIndices = qMolAdjList.get(bIndex);
		
		IBond[] aNeighbours = new IBond[ aNeighbourIndices.size() ];
		IBond[] bNeighbours = new IBond[ bNeighbourIndices.size() ];		
		
		for( int anI = 0; anI < aNeighbours.length; anI++ ) {
			aNeighbours[anI] = hsMol.getBond( aNeighbourIndices.get(anI) );
		}
		
		for( int bnI = 0; bnI < bNeighbours.length; bnI++ ) {
			bNeighbours[bnI] = qMol.getBond( bNeighbourIndices.get(bnI) );
		}
		
		// count matches for first neighbourhood to second
		for( IBond aB : aNeighbours ) {
			
			boolean matchFound = false;
			
			for( IBond bB : bNeighbours) {
				if( GenerateCompatibilityGraphEdges.nodesMatch(hsMol, aB, qMol, bB, true, atomTypeProperty) ) {
					matchFound = true;
					break;
				}
			}
			
			if( ! matchFound )
				score++;
			
		}
		
		// now to reverse match count
		for( IBond bB : bNeighbours ) {

			boolean matchFound = false;

			for( IBond aB : aNeighbours) {
				if( GenerateCompatibilityGraphEdges.nodesMatch(hsMol, aB, qMol, bB, true, atomTypeProperty) ) {
					matchFound = true;
					break;
				}
			}

			if( ! matchFound )
				score++;

		}
		
		return score;
		
	}
	
	
	
	private int[] calculateExtendedConnectivityValues( List<List<Integer>> adjList ) {
		
		int radius = 2;
		
		int[] ecValues = new int[ adjList.size() ];
		
		// initial values
		for( int b1 = 0; b1 < adjList.size(); b1++ ) {
			ecValues[b1] = adjList.get(b1).size();
		}
		
		// iterations
		for( int rad = 0; rad < radius; rad++ ) {
			int[] newEcValues = new int[ adjList.size() ];
			for( int b1 = 0; b1 < adjList.size(); b1++ ) {
				newEcValues[b1] = 0;
				for( int b2 = 0; b2 < adjList.get(b1).size(); b2++ ) {
					newEcValues[b1] += ecValues[ adjList.get(b1).get(b2) ];
				}
			}
			
			ecValues = newEcValues;
		}
		
		return ecValues;
	}
	
	

	private boolean nonredundancyBySelectionScore( MatchList nml, Collection<MatchList> cs ) {
		
		Map<Integer, Integer> nmlTypes = nml.getAtomTypeCounts();
		
		for( MatchList cml : cs ) {
			if( (cml.getScoreNeighbour() == nml.getScoreNeighbour()) && 
					(cml.getScoreEC() == nml.getScoreEC()) && 
					(cml.getScoreTopology() == nml.getScoreTopology()) 
				) {
				
				Map<Integer, Integer> cmlTypes = cml.getAtomTypeCounts();
				
				// get all possible keys
				Set<Integer> keys = new HashSet<Integer>( nmlTypes.keySet() );
				keys.addAll( cmlTypes.keySet() );
				
				for( Integer k : keys ) {
					
					if( cmlTypes.get(k) != nmlTypes.get(k) ) {
						//System.out.println("redundant match list identified");
						return false;
					}
				}
				
			}
		}
		
		return true;
	}
	
	
	
	
	
	
	
	/*public AllPairsShortestPaths(int[][] adjList, IAtomContainer container) {

        int n = container.getAtomCount();

        //this.container     = container;
        ShortestPaths[] shortestPaths = new ShortestPaths[n];

        // for each atom construct the ShortestPaths object
        for (int i = 0; i < n; i++) {
            shortestPaths[i] = new ShortestPaths(adjList, container, i);
        }
    }*/
	
	/*
	 * each "node" should have edge in mol A, edge in mol B and a heuristic score
	 * lower score is beter
	 */
	private void buildup( List<int[]> nodes ) {
		
		int worstScore = correspondences.get( correspondences.size() - 1 ).getScore();
		//int worstScore = correspondences.last().getScore();
		int n = 1;
		
		// upper bound to MCS 
		bu: while( n < Math.max(hsMol.getBondCount(), qMol.getBondCount()) ) {
			
			List<MatchList> nextCorrespondences = new ArrayList<MatchList>( correspondenceNumber ); 
			//TreeSet<MatchList> nextCorrespondences = new TreeSet<MatchList>( matchListComparator ); 
			
			// go through correspondence list
			for( MatchList matchList : correspondences ) {
				
				// look for compatible pairs
				for( int[] p : nodes ) {
					

					boolean prenewGenericConditions = ! matchList.contains( 0, p[0] ) &&
							! matchList.contains( 1, p[1] ) &&
							bondEquivalence(matchList.getPairs(), p[0], p[1]) &&
							nonredundancyBySelectionScore(matchList, nextCorrespondences); 
					
					if( ! prenewGenericConditions )
						continue;
					
					MatchList nml = new MatchList( matchList );
					/*nml.getPairs().add(p);
					nml.calculateScore();*/
					nml.addPair(p);
							
					boolean postnewGenericConditions = 
							(nextCorrespondences.size() < correspondenceNumber || nml.getScore() <= worstScore );
					
					if( ! postnewGenericConditions )
						continue;
					
					if( connectedMCS && ! bondConnection(matchList.getPairs(), p[0], p[1]) )
						continue;
					
					if( useTopoDistance && ! differenceOfTopologicalDistance( matchList.getPairs(), p[0], p[1] ) )
						continue;
					
					
					
					// test conditions passed - execute this block
					{
						if( nextCorrespondences.size() < correspondenceNumber ) {
							nextCorrespondences.add( nml );
						} else {
							nextCorrespondences.set( nextCorrespondences.size() - 1, nml );
						}
						
						// sort correspondences by score
						if( nml.getScore() != worstScore )
							Collections.sort( nextCorrespondences, matchListComparator );
						
						// update worst score to that of last match list
						worstScore = nextCorrespondences.get( nextCorrespondences.size() - 1 ).getScore();
						//worstScore = nextCorrespondences.last().getScore();
					}
					
					
				}
			}
			
			if( nextCorrespondences.isEmpty() ) {
				break bu;
			}
			
			++n;
			//System.out.println( n + " " + correspondences.size() );
			
			correspondences = nextCorrespondences;
		}
		
	}
	

	@Override
	public void search( IAtomContainer graph1, IAtomContainer graph2 ) {
		
		hsMol = graph1;
		qMol = graph2;
		
		
		
		
		
		//ConvenienceTools.countRings(hsMol);
		//ConvenienceTools.countRings(qMol);
		
		// Conflicts with SMARTS matching, so we ignore this
		/*defineAtomTypes(hsMol);
		defineAtomTypes(qMol);*/
		atomTypeProperty = null;
		
		if( ! useTopoDistance )
			maxTopologicalDifference = 1000;
		
		modProdTime = System.currentTimeMillis();
		
		modProd = null;
		try {
			modProd = new GenerateCompatibilityGraphEdges(
					hsMol, 
					qMol, 
					matchBonds, 
					useRaymondHeuristics,
					false,
					maxTopologicalDifference, 
					atomTypeProperty, 
					false 
			);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		modProdTime = System.currentTimeMillis() - modProdTime;
		
		int modProdNodeCount = modProd.getNodes().size();
		modProdEdgeDensity = (modProd.numberOfEdges ) / (double)( modProdNodeCount * (modProdNodeCount - 1) );
		
		
		long descriptorTime = System.currentTimeMillis();
		
		// adjacency lists for extended connectivities
		hsMolAdjList = ConvenienceTools.createBondAdjacencyList( hsMol );
		qMolAdjList = ConvenienceTools.createBondAdjacencyList( qMol );
		
		pathDistancesHsMol = PathTools.computeFloydAPSP( ConvenienceTools.bondAdjacencyMatrix(hsMol) );
		pathDistancesQMol = PathTools.computeFloydAPSP( ConvenienceTools.bondAdjacencyMatrix(qMol) );
		
		int[] hsMolEC = calculateExtendedConnectivityValues(hsMolAdjList);
		int[] qMolEC = calculateExtendedConnectivityValues(qMolAdjList);
		
		
		
		
		if( verbose )
			System.out.println( hsMol.getBondCount() + " " + qMol.getBondCount() + " mod prod nodes = " + modProd.getNodes().size() + "  edges = " + modProd.numberOfEdges );
		
		List<int[]> nodes = new ArrayList<int[]>( modProd.getNodes().size() );
		for( int[] pair : modProd.getNodes() ) {
			int[] newPair = new int[6];
			newPair[0] = pair[0];
			newPair[1] = pair[1];
			//newPair[2] = dNeighbour( hsMol.getBond(pair[0]), qMol.getBond(pair[1]) );  // first heuristic score
			newPair[2] = dNeighbour2( pair[0], pair[1] );  // first heuristic score
			newPair[3] = Math.abs( hsMolEC[ pair[0] ] - qMolEC[ pair[1] ] );  // second heuristic score
			newPair[4] = 0;  // third heuristic score (calculated on the fly)
			newPair[5] = newPair[2] + newPair[3] + newPair[4];  // initial heuristic score
			
			nodes.add( newPair );
		}
		
		// sort nodes in ascending order on heuristic scores
		Collections.sort( nodes, pairComparator );
		
		
		if( verbose ) {
			for( int[] ea : nodes ) {
				IBond bondA = hsMol.getBond( ea[0] );
				IBond bondB = qMol.getBond( ea[1] );
				System.out.println(  " [ " + ea[0] + "," + ea[1] + "," + ea[2] + "," + ea[3] + "," + ea[4] + "," + ea[5] +  "]  " + 
						bondA.getOrder() + " " + bondA.getAtom(0).getSymbol() + " " +  bondA.getAtom(1).getSymbol() + " " +
						bondB.getOrder() + " " + bondB.getAtom(0).getSymbol() + " " +  bondB.getAtom(1).getSymbol() + " " );
			}
		}
		
		/*try {
			Thread.sleep(1200);
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		*/
		
		descriptorTime = System.currentTimeMillis() - descriptorTime; 
		
		
		long startBuildupTime = System.currentTimeMillis();
		
		// initialise correspondences
		int corLim = Math.min( correspondenceNumber, nodes.size() );
		correspondences = new ArrayList<MatchList>( corLim );
		//correspondences = new TreeSet<MatchList>( matchListComparator );
		for( int p = 0; p < corLim; p++ ) {
			MatchList ml = new MatchList( nodes.get(p) );
			ml.calculateScore();
			correspondences.add(ml);
		}
		
		if( nodes.size() > 0 ) {
			buildup(nodes);
		} else {
			correspondences.add( new MatchList() );  // null placeholder
		}
		
		long buildupTime = System.currentTimeMillis() - startBuildupTime;


		Map<IBond, IBond> commonSubgraphMap = new HashMap<IBond, IBond>();

		for( int[] pair : correspondences.get(0).getPairs() ) {

			IBond hsBond = hsMol.getBond(pair[0]);
			IBond qBond = qMol.getBond(pair[1]);

			if( verbose ) {
				System.out.println( 
						hsMol.getAtomNumber( hsBond.getAtom(0) ) + "," +
								hsMol.getAtomNumber( hsBond.getAtom(1) )  + "=" +
								qMol.getAtomNumber( qBond.getAtom(0) ) + "," +
								qMol.getAtomNumber( qBond.getAtom(1) )
						);
			}
			
			commonSubgraphMap.put(hsBond, qBond);
		}


		Map<IAtom, IAtom> commonAtomMap = createAtomMapFromBondMap(graph1, graph2, commonSubgraphMap);

		// query to hyperstructure mapping
		List<Integer> mapping = atomMapToChromosome(hsMol, qMol, commonAtomMap);
		
		if( verbose ) {
			System.out.println( "mapping - " + commonAtomMap  );
			System.out.println( "chromosome - " + mapping  );
		}
		
		mcsAtomIsomorphisms.add( mapping );
		mcsBondIsomorphisms.add( commonSubgraphMap );


		//mcsBondIndexIsomorphisms.add( correspondences.first().getPairs() );
		mcsBondIndexIsomorphisms.add( correspondences.get(0).getPairs() );
		
		if( verbose ) {
			System.out.println( "Modular Product creation time = " + modProdTime );
			System.out.println( "Descriptor creation time = " + descriptorTime );
			System.out.println( "buildup time - " + buildupTime );
		}
	}
	
	
	
	
	
	protected boolean useTopoDistance = false;
	protected boolean connectedMCS = false;
	protected boolean useRaymondHeuristics = false;
	
	protected int maxTopologicalDifference = -1;
	protected int correspondenceNumber = 40;
	protected List<MatchList> correspondences;  // 
	
	/* 
	 * TODO  IMO the TreeSet method actually better than the original list solution, due to it being "constantly sorted".
	 * 
	 * Also due to sorting (and hashing) via heuristic score, much less solutions are stored (including less fit ones), providing a 
	 * better diversity of solutions & ultimately sometimes yielding better answers
	 */
	//protected TreeSet<MatchList> correspondences;   
	
	protected int[][] pathDistancesHsMol;
	protected int[][] pathDistancesQMol;
	protected List<List<Integer>> hsMolAdjList;
	protected List<List<Integer>> qMolAdjList;
	
	protected static String anyAtomID = "*";
	protected static String atomTypeProperty = "_at";
	
}
