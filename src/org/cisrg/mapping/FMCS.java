package org.cisrg.mapping;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;
import java.util.Map.Entry;

import org.cisrg.mapping.VentoFoggia.Mappings;
import org.cisrg.mapping.VentoFoggia.VentoFoggia2;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
 



/**
 * An adaptation of Andrew Dalke's fMCS Python algorithm into Java
 * 
 * @author Edmund Duesbury
 * @date January 2015
 */

public class FMCS extends MCSMethods {
	
	
	// descending order of bond score
	protected Comparator<IBond> bondScoreComparator = new Comparator<IBond>() {
		 public int compare(IBond o1, IBond o2) {

			 	int score1 = bondRingScore( o1 );
			 	int score2 = bondRingScore( o2 );
	            
	            return (score1 == score2) ? 0 : (score1 > score2) ? 1 : -1;
	            
	           
		 } 
	};
	
	
	// "descending" order of number of bonds - programmed as ascending, but PriorityQueue reverses the order (to make it descending)
		protected Comparator<ExtendableSubgraph> bondNumberComparator = new Comparator<ExtendableSubgraph>() {
			 public int compare(ExtendableSubgraph o1, ExtendableSubgraph o2) {

				 	int nb1 = o1.getInternalBonds().size();
				 	int nb2 = o2.getInternalBonds().size();
				 	
		           // int dComp;
		            if( nb1 == nb2 ) {
		            	//dComp = 0;
		            	
		            	int eb1 = o1.getVisitedBonds().size();
		            	int eb2 = o2.getVisitedBonds().size();
		            	
		            	return eb1 > eb2 ? 1 : -1;
		            	
		            	/*if( eb1 > eb2 ) {
		            		dComp = 1;
		            	} else if( eb1 < eb2 ) {
		            		dComp = -1;
		            	}*/
		            	
		            } /*else if( nb1 < nb2 ) {
		            	dComp = 1;
		            } else {
		            	dComp = -1;
		            }*/
		            
		            return nb1 < nb2 ? 1 : -1;
		            //return dComp;
			 } 
		};
		
	
	
	public FMCS() {
		
	}
	
	
	/**
	 * Is the query a subgraph of the target molecule?
	 * 
	 * XXX
	 * Uses a SMARTS (non-canonical) look-up table in an attempt to save time.  I'm unsure if this helps given how fast VentoFoggia subgraph isomorphism actually is.
	 * 
	 * @param query
	 * @param target
	 * @return
	 */
	private boolean isSubgraph( ExtendableSubgraph query, IAtomContainer target ) {
		
		boolean isSg = false;
		
		String smarts = query.getSMARTS();
		if( SMARTSMatch.containsKey(smarts) ) {
			isSg = SMARTSMatch.get( smarts );
			return isSg;
		}
		
		
		org.cisrg.mapping.VentoFoggia.Pattern pattern = VentoFoggia2.findSubstructure(query.getAtomContainer(), !matchBonds );
		
		isSg = pattern.matches(target);
		SMARTSMatch.put( smarts, isSg );
		
		
		return isSg;
		
	}
	
	
	/**
	 * Tests to see if a bond in one molecule has any matching bonds in another molecule
	 * 
	 * @param qMol
	 * @param oMol
	 * @param qBond
	 * @return
	 */
	private boolean molMatchesBond( IAtomContainer qMol, IAtomContainer oMol, IBond qBond ) {
		
		for( IBond bond : oMol.bonds() ) {
			
			boolean matches;
			if( oMol instanceof IQueryAtomContainer )
				matches = GenerateCompatibilityGraphEdges.nodesMatch(oMol, bond, qMol, qBond, matchBonds, null);
			else 
				matches = GenerateCompatibilityGraphEdges.nodesMatch(qMol, qBond, oMol, bond, matchBonds, null);
			
			if( matches )
				return true;
		}
		
		return false;
	}
	
	
	
	private int bondRingScore(IBond bond) {
		
		int rBond = ConvenienceTools.isRingBond(bond) ? 1 : 0;
		int rAtom1 = bond.getAtom(0).getFlag( CDKConstants.ISINRING ) ? 1 : 0;
		int rAtom2 = bond.getAtom(1).getFlag( CDKConstants.ISINRING ) ? 1 : 0;
		//int rBondArom = ConvenienceTools.isAromatic(bond) ? 1 : 0;
		
		//System.out.println(rBond + rAtom1 + rAtom2);
		return rBond + rAtom1 + rAtom2;
	}
	
	
	
	/**
	 * Create a power set via bit sets.
	 * 
	 * Credit to rolfl on StackOverflow for this code - http://stackoverflow.com/questions/20035547/generate-binary-representation-of-numbers-from-1-to-2k-1
	 * 
	 * @param bits
	 * @return
	 */
	private static final ArrayList<BitSet> powerset(final int bits) {
        int size = 1 << bits;
        ArrayList<BitSet> results = new ArrayList<BitSet>(size);
        for (int val = 1; val < size; val++) {
            BitSet bs = new BitSet(bits);
            results.add(bs);
            int v = val;
            int b = 0;
            while (v != 0) {
                if ( (v & 1) == 1) {
                   bs.set(b);
                }
                b++;
                v >>>= 1;
            }
        }
        return results;
    }
	
	
	/**
	 * Create a power set via hash sets.
	 * 
	 * Credit to Andrew Mao on StackOverflow for this code - http://stackoverflow.com/questions/1670862/obtaining-a-powerset-of-a-set-in-java
	 * 
	 * @param bits
	 * @return
	 */
	private static final <T> Set<Set<T>> powerset( Collection<T> originalSet ) {
		
		List<T> list = new ArrayList<T>(originalSet);
		int n = list.size();

		Set<Set<T>> powerSet = new HashSet<Set<T>>();

		for( long i = 0; i < (1 << n); i++) {
		    Set<T> element = new HashSet<T>();
		    for( int j = 0; j < n; j++ )
		        if( (i >> j) % 2 == 1 ) element.add(list.get(j));
		    powerSet.add(element); 
		}

		return powerSet;
		
	}
	
	/**
	 * A depth-first search to estimate the number of possible extendable bonds
	 * 
	 * Though - surely just identifying the fragment that the subgraph is in and taking the difference of
	 *  bonds between the fragment & subgraph will give the same answer?  Also subtract excluded bonds
	 * 
	 * @param esg
	 * @param adj
	 * @return
	 */
	private int findExtensionSize( ExtendableSubgraph esg, List<List<Integer>> adj ) {
		int numBonds = 0;
		
		Set<Integer> visitedBonds = new HashSet<Integer>( esg.getVisitedBonds() );
		//visitedBonds.addAll( esg.getVisitedBonds() );
		//stack.pop();
		
		for( Integer e : esg.getAdjacentBonds() ) {
			
			Stack<Integer> stack = new Stack<Integer>();
			stack.push( e );
			
			while( ! stack.isEmpty() ) {
				Integer edge = stack.pop();
				
				for( Integer nextEdge : adj.get(edge) ) {
					if( visitedBonds.contains( nextEdge ) )
						continue;
					
					numBonds++;
					visitedBonds.add( nextEdge );
					stack.push( nextEdge );
				}
			}
			
		}
		
		return numBonds;
	}
	
	/*
	private int findExtensionSize2( ExtendableSubgraph esg, IAtomContainer mol, IAtomContainerSet fragments, List<List<Integer>> adj ) {
		 
		
		// find which fragment it's in
		IBond fBond = mol.getBond( esg.getInternalBonds().get(0) );
		IAtomContainer currentFrag = null;
		for( IAtomContainer frag : fragments.atomContainers() ) {
			if( mol.contains(fBond) ) {
				currentFrag = frag;
				break;
			}
		}
		
		int extSize = currentFrag.getBondCount() - esg.getVisitedBonds().size();
		
		return extSize;
	}
	*/
	 
	
	
	private void saveSubgraph( ExtendableSubgraph sg ) {
		
			sg.createBondMaps();
			
			// FIXME  There's a better way of doing this.  Unfortunately I haven't found the fix yet for why some yield a size of 0
			if( sg.getCommonSubgraphMap() == null || sg.getCommonSubgraphMap().size() == 0 )
				return;
		
			if( bestSolutions.isEmpty() ) {
				bestSolutions.add( sg );
			} else if( sg.size() >= bestSolutions.get(0).size() ) {
				
					if( sg.getCommonSubgraphMap().size() > bestSolutions.get(0).getCommonSubgraphMap().size() ) {  // delete smaller solutions than newly-identified best
						bestSolutions.clear();
					}
					
					bestSolutions.add( sg );
				
			}
		
	}
	
	
	private int largestSubgraphSize() {
		
		if( ! bestSolutions.isEmpty() )
			return bestSolutions.get(0).size();
			
		return 0;
	}
	
	
	/**
	 * Create an adjacency list using the specified bond indices
	 * 
	 * Looks like I don't actually need this - quicker to just create the adjacency list from the modified IAtomContainer object.
	 * 
	 * @param mol
	 * @param fMol
	 * @param origAdj		adjacency list to use as template
	 * @param bondIndices	bond indices to keep
	 * @return
	 */
	@Deprecated
	private List<List<Integer>> getRefinedAdjacencyList( IAtomContainer mol, IAtomContainer fMol, List<List<Integer>> origAdj, List<Integer> bondIndices ) {
		List<List<Integer>> adj = new ArrayList<List<Integer>>( bondIndices.size() );
		
		for( Integer e : bondIndices ) {
			List<Integer> adjRow = new ArrayList<Integer>( origAdj.get(e).size() );
			//IBond bond1 = hsMol.getBond(e);
			for( Integer bi : origAdj.get(e) ) {
				IBond bond2 = mol.getBond(bi);
				
				
				int newIndex = fMol.getBondNumber( bond2 );
				
				if( newIndex >= 0 )
					adjRow.add( newIndex );
			}

			adj.add( adjRow );
		}
		
		return adj;
	}
	
	
	
	
	/**
	 * declare seeds array (set of bonds to start extending from)
	 * prune seeds
	 * find extensions:
		*	remove last element from seeds
		* 	create new (larger) subgraphs - add to list of seeds
		* 	track largest common fragments
	 * @throws CDKException 
	 */
	private void enumerateSubgraphs() throws CDKException {
		
		IAtomContainer smallestLargestFragment = null;
		bestSolutions = new ArrayList<ExtendableSubgraph>();
		SMARTSMatch = new HashMap<String, Boolean>();
		
		/*hsAdj = ConvenienceTools.createBondAdjacencyList(hsMol);
		qAdj = ConvenienceTools.createBondAdjacencyList(qMol);*/
		
		// List of subgraphs to extend
		PriorityQueue<ExtendableSubgraph> extendables = new PriorityQueue<ExtendableSubgraph>( hsMol.getBondCount() + 1, bondNumberComparator );
		
		hsBondIndices = new ArrayList<Integer>( hsMol.getBondCount() );
		qBondIndices = new ArrayList<Integer>( qMol.getBondCount() );
		
		for( int n = 0, bs = hsMol.getBondCount(); n < bs; n++ ) {
			if( molMatchesBond( hsMol, qMol, hsMol.getBond(n) ) )
				hsBondIndices.add(n);
		}
		
		for( int n = 0, bs = qMol.getBondCount(); n < bs; n++ ) {
			if( molMatchesBond( qMol, hsMol, qMol.getBond(n) ) )
				qBondIndices.add(n);
		}
		
		if( verbose ) {
			System.out.println( hsBondIndices + " " + hsBondIndices.size() );
			System.out.println( qBondIndices + " " + qBondIndices.size() );
		}
		
		// no MCS possible
		if( hsBondIndices.isEmpty() || qBondIndices.isEmpty() )
			return;
		
		
		swapMols = true;
		
		// obtain fragments
		mol1 = ConvenienceTools.createSubgraph( hsMol, hsBondIndices );
		mol2 = ConvenienceTools.createSubgraph( qMol, qBondIndices );
		
		/*// want mol1 (the thing we're creating subgraphs from) to be the smallest molecule
		if( mol1.getBondCount() > mol2.getBondCount() ) {
			System.out.println("old swap " + mol1.getBondCount() + " " + mol2.getBondCount() );
			IAtomContainer tempMol = mol1;
			mol1 = mol2;
			mol2 = tempMol;
		}*/
		
		
		IAtomContainerSet mol1Frags = ConvenienceTools.partitionIntoMolecules(mol1);
		IAtomContainerSet mol2Frags = ConvenienceTools.partitionIntoMolecules(mol2);
		
		
		
		// delete excluded edges from the adjacency representations
		//List<List<Integer>> adj1 = getRefinedAdjacencyList(hsMol, mol1, hsAdj, hsBondIndices);
		List<List<Integer>> adj1 = ConvenienceTools.createBondAdjacencyList(mol1);
		List<List<Integer>> adj2 = ConvenienceTools.createBondAdjacencyList(mol2);
		
		
		// which molecule has the smallest of largest fragments?
		
		
		int maxFSize1 = 0;
		int maxFSize2 = 0;
		IAtomContainer slf1 = null, slf2 = null;
		
		try {
		for( IAtomContainer mol : mol1Frags.atomContainers() ) {
			ConvenienceTools.countRings(mol);
			
			if( mol.getBondCount() > maxFSize1 ) {
				maxFSize1 = mol.getBondCount();
				slf1 = mol;
			}
		}
		
		for( IAtomContainer mol : mol2Frags.atomContainers() ) {
			ConvenienceTools.countRings(mol);
			
			int bCount = mol.getBondCount();
			if( bCount > maxFSize2 ) {
				maxFSize2 = bCount;
				slf2 = mol;
			}
		}
		
		} catch ( IllegalArgumentException e2 ) {
			e2.printStackTrace();
		}
		
		// if there's no smallest-largest fragment, there's no MCS
		if( slf1 == null || slf2 == null ) {
			saveSubgraph( new ExtendableSubgraph( new ArrayList<Integer>(), null, null) );
			return;
		}
		
		// choose the smallest of the largest fragments found
		if( slf2.getBondCount() > slf1.getBondCount() )
			smallestLargestFragment = slf1;
		else
			smallestLargestFragment = slf2;
		
		if( smallestLargestFragment instanceof IQueryAtomContainer )
			ConvenienceTools.correctAtomBondTypes(smallestLargestFragment);
		
		
		if( maxFSize2 > maxFSize1 || (hsMol instanceof IQueryAtomContainer) )
			swapMols = false;
		
		if( swapMols ) {
			IAtomContainer tempMol = mol1;
			mol1 = mol2;
			mol2 = tempMol;
			
			IAtomContainerSet tempSet = mol1Frags;
			mol1Frags = mol2Frags;
			mol2Frags = tempSet;
			
			List<List<Integer>> tempAdj = adj1;
			adj1 = adj2;
			adj2 = tempAdj;
		}
		
		
		
		if(verbose) {
			System.out.println("swap molecules - " + swapMols);
			System.out.println( sGenerator.create(mol1) );
			System.out.println( sGenerator.create(mol2) );
		}
		
		// smallest largest fragment check
		// if the smallest largest fragment is a valid subgraph, save the solution and exit
		if( ConvenienceTools.isSubgraph(smallestLargestFragment, mol1) && ConvenienceTools.isSubgraph(smallestLargestFragment, mol2) ) {
			
			Pattern pattern = VentoFoggia.findSubstructure( smallestLargestFragment );
			int[] mapping = pattern.match(mol1);
			
			
			HashMap<IAtom, IAtom> commonAtomMap = new HashMap<IAtom, IAtom>();
			for( int n = 0; n < mapping.length; n++ ) {
				
				IAtom hsAtom = smallestLargestFragment.getAtom( n );
				IAtom qAtom = mol1.getAtom( mapping[n] );
				commonAtomMap.put( hsAtom, qAtom );

			}
			
			ArrayList<Integer> internals = new ArrayList<Integer>( mapping.length );
			 
			Map<IBond, IBond> commonSubgraphMap = ConvenienceTools.makeBondMapOfAtomMap(smallestLargestFragment, mol1, commonAtomMap, true);
			for( IBond b : commonSubgraphMap.values() ) {
				internals.add( mol1.getBondNumber( b ) );
			}
			
			if(verbose)
				System.out.println("Smallest-largest fragment match");
			
			saveSubgraph( new ExtendableSubgraph(internals, null, null) );
			return;
		}
		
		
		
		Set<Integer> initialVisited = new HashSet<Integer>( adj1.size() );
		
		// initial seeds - set of single bond subgraphs which have been sorted in canonical order, from one of the molecules
		// use the molecule with the smallest of largest fragments as the point of enumeration
		for( IAtomContainer frag : mol1Frags.atomContainers() ) {
			
			List<IBond> bondList = new ArrayList<IBond>( frag.getBondCount() );
			for( IBond fBond : frag.bonds() ) { bondList.add( fBond );  /*System.out.println( frag.getBondNumber(fBond) + " " + bondRingScore(fBond) );*/ }
			Collections.sort( bondList, bondScoreComparator );
			//Collections.reverse(bondList);
			
			for( IBond fBond : bondList ) {
				int index = mol1.getBondNumber(fBond);
				
				List<Integer> internals = new ArrayList<Integer>( maxFSize1 );
				internals.add(index);
				initialVisited.add(index);
				List<Integer> externals = adj1.get(index);
				
				// things are added to the queue in equal priority (number of internal bonds)
				extendables.add( new ExtendableSubgraph(internals, externals, new HashSet<Integer>(initialVisited) ) );
			}
		}
		
		
		
		// now to track through the seeds, and extend
		while( ! extendables.isEmpty() ) {
			
			++numberOfSteps;
			
			if( numberOfSteps % 100 == 0 ) {
				long currentTime = System.currentTimeMillis();
				if( currentTime - mcsStartTime > expansionTimeLimitMs ) {
					throw new CDKException("Error - time limit of " + expansionTimeLimitMs + " ms reached!");
				}
			}
			
			ExtendableSubgraph es = extendables.poll( );  // removes from head
			
			Set<Integer> visitedBonds = new HashSet<Integer>( es.getVisitedBonds() );
			visitedBonds.addAll( es.getInternalBonds() );
			visitedBonds.addAll( es.getAdjacentBonds() );  // we'll be visiting each of these too, so mark them as visited
			
			
			
			// create power sets of all bonds to extend from
			// Bitset manipulation seems pretty fast (in Java) for the purposes of power sets
			Collection<BitSet> extCombos = powerset( es.getAdjacentBonds().size() );
			//Collection<Set<Integer>> extCombos = powerset( es.getAdjacentBonds() );
			
			for( BitSet extBs : extCombos ) {
				
				// continue if there're no new extensions
				 /*if( extBs.cardinality() == 0 )
					 continue;*/
				
				// BitSet to bond indices
				List<Integer> bondsToAdd = new ArrayList<Integer>( es.getAdjacentBonds().size() );
				 for (int i = extBs.nextSetBit(0); i >= 0; i = extBs.nextSetBit(i+1)) {
					 bondsToAdd.add( es.getAdjacentBonds().get(i) );
				 }
				 
				 
			 	 // new internals - basically the new subgraph
				 List<Integer> nInternals = new ArrayList<Integer>( es.getInternalBonds() );
				 nInternals.addAll( bondsToAdd );
				 
				 // new externals - neighbours of original externals, but removing the "visited" bonds.
				 // must be unique
				 Set<Integer> nExternals = new HashSet<Integer>( adj1.size() );
				 //for( Integer e : es.getAdjacentBonds() ) {  // interestingly this "bug" can yield useful disconnected answers
				 for( Integer e : nInternals ) {
					 nExternals.addAll( adj1.get(e) );
				 }
				 nExternals.removeAll( visitedBonds );
				 
				 
				 ExtendableSubgraph newEs = new ExtendableSubgraph( nInternals, new ArrayList<Integer>( nExternals ), visitedBonds );
				 //System.out.println( "subgraph - " + newEs + " " + es.getInternalBonds() + " " + newEs.getInternalBonds() );
				 
				 
				 
				 // pruning - don't continue if extending this subgraph won't yield anything of comparable size to the current best solution(s)
				 int extSize = findExtensionSize(newEs, adj1);
				 //int extSize2 = findExtensionSize2(newEs, mol1, mol1Frags, adj1);
				 int largestPossibleSize = extSize + nInternals.size();
				 
				 if( verbose )
					 System.out.println( "possible ext size & int size & largest = " + extSize +  " " + newEs.size() + " " + largestPossibleSize + " " + largestSubgraphSize() + " " + smallestLargestFragment.getBondCount() );
				 
				 if( largestPossibleSize <= largestSubgraphSize() && bestSolutions.size() > 1 )
					 continue;
				 
				 // does extension match target?
				 /* XXX NOTE - we don't bother with the canonical SMARTS in this algorithm. 
				  * The VF algorithm in CDK is so fast at checking for subgraph isomorphism that I felt it wasn't
				  * necessary to use anything else
				  */
				 boolean matches = isSubgraph(newEs, mol2);
				 //boolean matches = isSubgraph(newEs.getAtomContainer(), mol2);
				 
				 /*if( largestPossibleSize >= 12 ) {
					 System.out.println( newEs + " size = " + newEs.size() + " " + matches );
				 }*/
				
				 if( verbose )
					 System.out.println( "matches - " + matches + "  " + newEs );
				 
				 if( ! matches )
					 continue;
				 
				 
				 

				 
				 extendables.add( newEs );
				 
				 saveSubgraph( newEs );
				 
			}
		}
		
		
		if( verbose ) {
			System.out.println( "fMCS best solutions: " );
			System.out.println( bestSolutions );
		}
		
		
		
	}

	
	@Override
	public void search(IAtomContainer graph1, IAtomContainer graph2) {
		
		
		
		/*try {
			ConvenienceTools.initializeMolecule(hsMol);
			ConvenienceTools.initializeMolecule(qMol);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
		numberOfSteps = 0;
		//timeOut = false;
		
		try {
			enumerateSubgraphs();
		} catch (CDKException e) {
			e.printStackTrace();
		} catch (IllegalArgumentException e2) {
			System.err.println("Illegal Argument Exception in fMCS: ");
			e2.printStackTrace();
			System.err.println("Error end. ");
		}
		
		
		for( ExtendableSubgraph es : bestSolutions ) {
			
			es.createBondMaps();
			Map<IBond, IBond> commonSubgraphMap = es.getCommonSubgraphMap();
			Map<IAtom, IAtom> commonAtomMap = es.getCommonAtomMap();
			
			
			
			
			List<int[]> pairList = new ArrayList<int[]>( commonSubgraphMap.size() );
			for( Entry<IBond, IBond> e : commonSubgraphMap.entrySet() ) {
				int[] pair = new int[]{ hsMol.getBondNumber(e.getKey()), qMol.getBondNumber(e.getValue()) };
				pairList.add(pair);
			}
			mcsBondIndexIsomorphisms.add( pairList );
			//System.out.println( pairList );

			mcsAtomIsomorphisms.add( atomMapToChromosome( hsMol, qMol, commonAtomMap) );
			mcsBondIsomorphisms.add( commonSubgraphMap );
		}
		
		
		
	}
	
	/*
	public static void main( String[] args ) {
		
		FMCS fmcs = new FMCS();
		//String smiles1 = "CNc1nccc2n(C)c3c(ncnc3c12)N1CCN(CCc2ccc(F)c(F)c2)CC1";
		String smiles1 = "c1ccc2c(c1)ccc1cc3c(ccc4ccccc34)cc21";
		String smiles2 = "Cc1c2cccc3CC4c5cccc6c7CNCc8ccc9c%10ccc%11ccc1c1c(c4c(c%10c%111)c(c56)c9c78)c23";
		//String smiles2 = "CCOC(=O)Cn1c2ccc(cc2c2ncnc(N3CCN(CCc4ccc(F)c(F)c4)CC3)c12)[N+]([O-])=O";
		
		long mcsTime = System.currentTimeMillis();
		try {
			IAtomContainer m1 = fmcs.sParser.parseSmiles( smiles1 );
			IAtomContainer m2 = fmcs.sParser.parseSmiles( smiles2 );
			fmcs.setMainMol(m1);
			fmcs.setQueryMol(m2);
			System.out.println( ConvenienceTools.isSubgraph(m1, m2) );
			 fmcs.enumerateSubgraphs();
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println( "fMCS took " + ( System.currentTimeMillis() - mcsTime) + " milliseconds" );
		
	}
	*/
	
	
	class ExtendableSubgraph {
		
		public ExtendableSubgraph( List<Integer> internals, List<Integer> externals, Collection<Integer> visited ) {
			internalBonds = internals;
			adjacentBonds = externals;
			visitedBonds = visited;
			smarts = null;
			container = null;
		}
		
		public IAtomContainer getAtomContainer() {
			if( container == null ) {
				container = ConvenienceTools.createSubgraph(mol1, internalBonds);
				//ConvenienceTools.calculateImplicitHydrogens(container);
			}
			
			return container;
		}
		
		public String getSMARTS() {
			
			if( smarts == null ) {
				IAtomContainer sg = getAtomContainer();
				
				
				if( sg instanceof QueryAtomContainer ) {
					smarts = smaH.toSmarts( (QueryAtomContainer) sg );
				} else {
					try {
						smarts = sGenerator.create(sg);
					} catch (CDKException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			
			return smarts;
			
		}
		
		
		public void createBondMaps() {
			

			commonSubgraphMap = new HashMap<IBond, IBond>();
			commonAtomMap = new HashMap<IAtom, IAtom>();
			
			if(verbose)
				System.out.println( "final sg - " +  this.toString() + " " + size() );
			
			/*IAtomContainer esac = es.getAtomContainer();
			
			try {
				ConvenienceTools.initializeMolecule(esac);
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}*/
			
			org.cisrg.mapping.VentoFoggia.Pattern pattern = VentoFoggia2.findSubstructure( this.getAtomContainer(), !matchBonds );
			
			//pattern.matchAll(target)
			
			// atom mappings
			int[] mapping1 = pattern.match(hsMol);
			int[] mapping2 = pattern.match(qMol);
			 
			/*
			 * optimistically test first mapping to see if all the bonds are captured.
			 * 
			 * If not, test all mappings instead
			 */
			if( ! bondMapFromAtomMap( mapping1, mapping2 ) ) {
				
				Mappings mappings1 = pattern.matchAll(hsMol);
				Mappings mappings2 = pattern.matchAll(qMol);
				
				outer: for( int[] mappingHs : mappings1 ) {
					for( int[] mappingQ : mappings2 ) {
						if( bondMapFromAtomMap( mappingHs, mappingQ ) )
							break outer;
					}
				}
			}
			
			
			
			if( verbose )
				System.out.println( "common sg bond count - " + commonSubgraphMap.size() + " | " + mapping1.length + " " + mapping2.length );
			
		}
		
		
		private boolean bondMapFromAtomMap( int[] mapping1, int[] mapping2 ) {
			
			commonAtomMap.clear();
			commonSubgraphMap.clear();
			
			// avoid null maps
			if (mapping1.length == 0 || mapping2.length == 0)
				System.err
						.println("Warning - no mapping for subgraph identified!");

			for (int n = 0; n < mapping1.length; n++) {
				IAtom at1 = hsMol.getAtom(mapping1[n]);
				IAtom at2 = qMol.getAtom(mapping2[n]);
				commonAtomMap.put(at1, at2);
			}

			// restore original bond indices (from potentially fragmented
			// version of mol1 or mol2)
			List<Integer> esBonds = new ArrayList<Integer>(size());
			if (swapMols) {
				for (Integer bInd : getInternalBonds()) {
					esBonds.add(qBondIndices.get(bInd));
				}
			} else {
				for (Integer bInd : getInternalBonds()) {
					esBonds.add(hsBondIndices.get(bInd));
				}
			}

			commonSubgraphMap = ConvenienceTools.bondMapFromOtherGraph( commonAtomMap, hsMol, qMol, esBonds, swapMols );

			if( commonSubgraphMap.size() == esBonds.size() )
				return true;
			
			return false;
			
		}
		
		
		// return SMARTS string
		public String toString() {
			return getSMARTS();
		}
		
		public int size() {
			return internalBonds.size();
		}
		
		public List<Integer> getAdjacentBonds() {
			return adjacentBonds;
		}
		
		public List<Integer> getInternalBonds() {
			return internalBonds;
		}
		
		public Collection<Integer> getVisitedBonds() {
			return visitedBonds;
		}
		
		public Map<IBond, IBond> getCommonSubgraphMap() {
			return commonSubgraphMap;
		}

		public Map<IAtom, IAtom> getCommonAtomMap() {
			return commonAtomMap;
		}
		
		
		private List<Integer> adjacentBonds, internalBonds;
		private Collection<Integer> visitedBonds;  
		private String smarts;
		private IAtomContainer container;
		private Map<IBond, IBond> commonSubgraphMap;
		private Map<IAtom, IAtom> commonAtomMap;
	}

	
	//private List<List<Integer>>  hsAdj, qAdj;
	private boolean swapMols;
	private IAtomContainer mol1, mol2;
	private List<ExtendableSubgraph> bestSolutions;
	private Map<String, Boolean> SMARTSMatch;


	private List<Integer> hsBondIndices;


	private List<Integer> qBondIndices;
	
}
