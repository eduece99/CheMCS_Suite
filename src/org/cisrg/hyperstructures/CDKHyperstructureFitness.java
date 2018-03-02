package org.cisrg.hyperstructures;

//import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.cisrg.knime.ExtendedIsomorphism;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.matchers.smarts.AliphaticSymbolAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AnyOrderQueryBond;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticOrSingleQueryBond;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.silent.AtomContainer;
//import org.openscience.cdk.interfaces.IBond.*;
//import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;




public class CDKHyperstructureFitness extends GAPlugins {
	
	/**
	 * Don't invoke this constructure - only override!  (This constructur does absolutely nothing, only a place-holder for overriding constructors in subclasses)
	 */
	protected CDKHyperstructureFitness() {
		
	}
	
	public CDKHyperstructureFitness( IAtomContainer hs, IAtomContainer qmol ) {
		this.hs = hs;
		this.queryMol = qmol;
		
		random = new Random();
		
		int hsAtomCount = hs.getAtomCount();
		int qAtomCount = queryMol.getAtomCount();
		
		/*
		 *  To save the chemistry API from searching for the neighbours of an atom in a molecule each time
		 *  we need to find the neighbours, we'll find them first them store them for later as we know
		 *  the molecule's won't change until we modify the hyperstructure - necessary for fitness function
		 */
		hAtomBonds = new HashMap<Integer, List<IBond>>();
		qAtomBonds = new HashMap<Integer, List<IBond>>();
		
		// Create map of atomic numbers - lists of all atoms per atomic number
		hAtomsMap = new HashMap<Integer, ArrayList<Integer>>();  
		qAtomsMap = new HashMap<Integer, ArrayList<Integer>>();  
		
		for( int n = 0; n < hsAtomCount; n++ ) {
			int atNum = hs.getAtom(n).getAtomicNumber();
			
			if( ! hAtomsMap.containsKey( atNum ) ) {
				hAtomsMap.put(atNum, new ArrayList<Integer>( hsAtomCount ) );
			}
			hAtomsMap.get( atNum ).add( n );
			
			hAtomBonds.put( n , hs.getConnectedBondsList( hs.getAtom(n) ) );
		}
		
		for( int n = 0; n < qAtomCount; n++ ) {
			int atNum = queryMol.getAtom(n).getAtomicNumber();

			// dummy list for any atom types in the query that don't exist in the hyperstructure.  Allows dummy matches.
			if( ! hAtomsMap.containsKey( atNum ) ) {
				hAtomsMap.put(atNum, new ArrayList<Integer>( hsAtomCount ) );
			}
			
			if( ! qAtomsMap.containsKey( atNum ) ) {
				qAtomsMap.put(atNum, new ArrayList<Integer>( (int) qAtomCount ) );
			}
			qAtomsMap.get( atNum ).add( n );
			
			qAtomBonds.put( n , queryMol.getConnectedBondsList( queryMol.getAtom(n) ) );
		}
		
		// We now need to make the sizes of the hAtomsMap ArrayLists equal - fill the excess
		// with asterixes to allow for null mappings
		
		Integer[] qKeys = qAtomsMap.keySet().toArray( new Integer[0] );
		
		for( int k = 0; k < qKeys.length; k++ ) {
			
			int qAtSize = qAtomsMap.get( qKeys[k] ).size();
			int hAtSize = hAtomsMap.containsKey( qKeys[k] ) ? hAtomsMap.get( qKeys[k] ).size() : 0; 
			
			if( qAtSize > hAtSize ) {
				for( int s = 0; s < (qAtSize - hAtSize); s++ ) {
					hAtomsMap.get( qKeys[k] ).add( -1 - s );
				}
			}
		}
		
		Integer[] hKeys = hAtomsMap.keySet().toArray( new Integer[0] );
		
		// placeholder so that anything can potentially "not match"
		for( int k = 0; k < hKeys.length; k++ ) {
			hAtomsMap.get( hKeys[k] ).add( -99 );
		}
		
		// set atom ids in the query for the fitness function to use
		for( IAtom at : queryMol.atoms() ) {
			at.setID( queryMol.getAtomNumber(at) + "" );;
		}
		
		// Find number of rings in the molecules
		ringFinder = new AllRingsFinder();
		ringFinder.setTimeout(1000);
		sg = new SmilesGenerator();
		/*
		try {
			//hsRings = ringFinder.findAllRings(hs).getAtomContainerCount();
			
			String hsSMILES = sg.createSMILES(hs);
			Pattern hsP = Pattern.compile("%\\d+");
			Matcher hsM = hsP.matcher(hsSMILES);

			int maxSize = 9;
			while( hsM.find() ) {
				String group = hsM.group().substring(1);
				int tempSize = Integer.parseInt(group);
				
				if( tempSize > maxSize)
					maxSize = tempSize;
			}
			hsRings = maxSize;
			
			queryRings = ringFinder.findAllRings(queryMol).getAtomContainerCount();
			
			System.out.println(hsRings + " @ ");
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();

			
			
			try {
				queryRings = ringFinder.findAllRings(queryMol).getAtomContainerCount();
			} catch (CDKException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		} */
	}
	
	protected IAtom[] getAtomNeighbours( IAtomContainer mol, int atIndex ) {
		/*
		Atom hsAtom = mol.getAtomWithIdx( atIndex );
		Bond_Vect atomBonds = hsAtom.getBonds();
		Atom[] neighbours = new Atom[ (int) atomBonds.size() ];
		
		for( int b = 0; b < atomBonds.size(); b++ ) {
			int index1 = (int) atomBonds.get(b).getBeginAtomIdx();
			int index2 = (int) atomBonds.get(b).getEndAtomIdx();
			
			if( index1 != atIndex ) {
				neighbours[b] = mol.getAtomWithIdx( index1 );
			} else {
				neighbours[b] = mol.getAtomWithIdx( index2 );
			}
		}
		
		atomBonds.delete();
		*/
		List<IAtom> atomNeighbours = mol.getConnectedAtomsList( mol.getAtom(atIndex) );
		return atomNeighbours.toArray( new IAtom[0] );
	}
	
	protected int[] getAtomNeighboursIndices( IAtomContainer mol, int atIndex ) {
		IAtom[] neighbours = getAtomNeighbours( mol, atIndex );
		int[] n2 = new int[ neighbours.length ];
		
		for( int n = 0; n < neighbours.length; n++ ) {
			n2[n] = mol.getAtomNumber( neighbours[n] );
		}
		
		return n2;
	}
	
	private IAtom getAtomFromId( IAtomContainer mol, String id ) {
		for( IAtom at : mol.atoms() ) {
			if( at.getID().equals(id) )
				return at;
		}
		
		return null;
	}
	
	
	public IAtomContainer calculateMCS( ArrayList<Integer> chromosome ) {

		IAtomContainer mcs = new AtomContainer();
		/*
		for( IAtom at : hs.atoms() ) {
			at.setID( hs.getAtomNumber(at) + "" );
			//System.out.print( at.getAtomicNumber() + " " );
			
		}
		*/

		for( IAtom at : queryMol.atoms() ) {
			mcs.addAtom( at );
		}
		
		for( int i = 0; i < chromosome.size(); i++ ) {
			if( chromosome.get(i) >= 0 ) {
				
				List<IBond> qatomBonds = qAtomBonds.get( i );  // very slow
				
				// go through all the bonds and add bonds where both atoms exist in the chromosome
				for( IBond b : qatomBonds ) {
					int qat1 = Integer.parseInt( b.getAtom(0).getID() );
					int qat2 = Integer.parseInt( b.getAtom(1).getID() );
					//int qat1 = queryMol.getAtomNumber( b.getAtom(0) );
					//int qat2 = queryMol.getAtomNumber( b.getAtom(1) );
					
					if( chromosome.get(qat1) >= 0 && chromosome.get(qat2) >= 0 ) {
						//IAtom mAt1 = getAtomFromId( mcs, b.getAtom(0).getID() );
						//IAtom mAt2 = getAtomFromId( mcs, b.getAtom(1).getID() );
						
						IBond hsBond = hs.getBond( hs.getAtom(chromosome.get(qat1)), hs.getAtom(chromosome.get(qat2) ));
						
						if( hsBond != null ) {
							if( mcs.getBond( mcs.getAtom(qat1), mcs.getAtom(qat2)) == null )
								mcs.addBond( qat1, qat2, Order.SINGLE );
						}
					}
				}
				
			}
		}
		
		return mcs;
	}
	
	/**
	 * Add all query atoms that have been mapped
	 * 
	 * Find the bonds between the query atoms mapped to generate the MCS
	 * 
	 * @param chromosome
	 * @return
	 */
	public double fitnessAllBonds(ArrayList<Integer> chromosome) {
		
		IAtomContainer mcs = calculateMCS(chromosome);
		
		/*
		for( int i = 0; i < chromosome.size(); i++ ) {
			for( int j = 0; j < chromosome.size(); j++ ) {
				if( j != i && chromosome.get(i) >= 0 && chromosome.get(j) >= 0 ) {
					
					IBond hsBond = hs.getBond( hs.getAtom(chromosome.get(i)), hs.getAtom(chromosome.get(j) ));
					if( hsBond != null ) {
						//if( mcs.getBond( mcs.getAtom(i), mcs.getAtom(j)) == null )
						
						IBond qBond = queryMol.getBond( queryMol.getAtom(i), queryMol.getAtom(j ));
						
						if( qBond != null && qBond.getOrder().equals(hsBond.getOrder()) ) {
							
							if( mcs.getBond( mcs.getAtom(i), mcs.getAtom(j)) == null )
								mcs.addBond( i,  j, hsBond.getOrder() );
						}
					}
					
				}
			}
		}
		*/
		
		/*
		for( int i = 0; i < chromosome.size(); i++ ) {
			//System.out.println( hs.getAtomNeighbors( hs.getAtomWithIdx(1) ).getClass() );

			if( chromosome.get(i) >= 0 ) {
				//IAtom hsAtom = hs.getAtom( chromosome.get(i) );
				List<IBond> atomBonds = hAtomBonds.get( chromosome.get(i) );  // very slow
				List<IBond> qatomBonds = qAtomBonds.get( i );  // very slow
				
				if( atomBonds != null ) {
					for( int b = 0; b < atomBonds.size(); b++ ) {
						
						int index1 = hs.getAtomNumber( atomBonds.get(b).getAtom(0) );
						int index2 = hs.getAtomNumber( atomBonds.get(b).getAtom(1) );
						
						if( index1 != chromosome.get(i) && chromosome.contains(index1) ) {
						
							IAtom at1 = getAtomFromId( mcs, atomBonds.get(b).getAtom(0).getID() );
							IAtom at2 = getAtomFromId( mcs, atomBonds.get(b).getAtom(1).getID() );
							
							if( at1 == null ) {
								at1 = atomBonds.get(b).getAtom(0);
							} 
							
							if( at2 == null ) {
								at2 = atomBonds.get(b).getAtom(1);
							}
							
							int id1 = Integer.parseInt( at1.getID() );
							int id2 = Integer.parseInt( at2.getID() );
							
							int qat1 = chromosome.indexOf( index1 );
							int qat2 = chromosome.indexOf( index2 );
							
							if( qat1 < 0 || qat2 < 0 )
								continue;
							
							System.out.println( (1+qat1) + " -|- " + (1+qat2) );
							
							if( queryMol.getBond( queryMol.getAtom(qat1), queryMol.getAtom(qat2) ) != null ) {
								//System.out.print( " | " + b2.getBeginAtom().getSymbol() + b2.getBondTypeAsDouble() + b2.getEndAtom().getSymbol() + "(" + qat1 + " " + qat2 + ")" );
								System.out.print( " | " + "(" + (qat1 + 1) + " " + (qat2 + 1) + ")" );
								mcs.addBond( mcs.getAtomNumber(at1), mcs.getAtomNumber(at2), Order.SINGLE);
							}
							
							
						}
					}
				}
			}
			
		}*/
		
		int mcsRingCount = 0;
		try {
			mcsRingCount = ringFinder.findAllRings(mcs).getAtomContainerCount();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//double fitness = ( mcs.getBondCount() + mcsRingCount ) / (double) ( queryMol.getBondCount() + queryRings );
		double fitness = ( mcs.getBondCount() ) / (double) ( queryMol.getBondCount() );
		
		/*
		if( fitness > 0.099 ) {
			SmilesGenerator sg = new SmilesGenerator(true);
			System.out.println( sg.createSMILES(mcs) + " " + fitness + " " + mcs.getBondCount() + " " + mcsRingCount );
		}
		*/
		return fitness;
		
	}
	
	
	/**
	 * Returns the number of correctly mapped bonds, in which the atom correspondences matter
	 * (that is, atoms must match with respect to chromosome order, rather than just forming a connected graph).
	 */
	public double fitnessStrict(ArrayList<Integer> chromosome) {
		
		/*
		 * Look for mapped atoms to the hyperstructure from the query, and see how many
		 * "bonds" are present from the chromosome - note which pairs of atoms these are
		 */
		ArrayList<int[]> hsFoundBonds = new ArrayList<int[]>( hs.getBondCount() * 1 );
		//ArrayList<long[]> qFoundBonds = new ArrayList<long[]>( (int) hs.getNumBonds() * 5 );
		
		for( int i = 0; i < chromosome.size(); i++ ) {
			//System.out.println( hs.getAtomNeighbors( hs.getAtomWithIdx(1) ).getClass() );

			if( chromosome.get(i) >= 0 ) {
				//IAtom hsAtom = hs.getAtom( chromosome.get(i) );
				List<IBond> atomBonds = hAtomBonds.get( chromosome.get(i) );  // very slow
				
				if( atomBonds != null ) {
					for( int b = 0; b < atomBonds.size(); b++ ) {
						int index1 = hs.getAtomNumber( atomBonds.get(b).getAtom(0) );
						int index2 = hs.getAtomNumber( atomBonds.get(b).getAtom(1) );
						
						if( index1 != chromosome.get(i) && chromosome.contains(index1) ) {
							//hsFoundBonds.add( atomBonds.get(b).getIdx() );
							hsFoundBonds.add( new int[]{index1, index2} );
							//qFoundBonds.add( new long[]{i, chromosome.indexOf(index1)} );
						} //else if( index2 != chromosome.get(i) && chromosome.contains(index2) ) {
							//hsFoundBonds.add( new long[]{index1, index2} );
							//qFoundBonds.add( new long[]{i, chromosome.indexOf(index2)} );
						//}
						
					}
				}
				//hsAtom = null;
				atomBonds = null;
			}
			
		}
		
		/*
		 * This part verifies if the identified bonds have corresponding bonds in the query,
		 * instead of simply being atom pairs (tests if the bonds map back onto the query molecule basically)
		 */
		HashSet<int[]> hfb = new HashSet<int[]>( hsFoundBonds );
		//System.out.println( hfb.size() + " " + hsFoundBonds.size() );
		int count = 0;
		for( int[] i : hfb ) {
			
			int qat1 = chromosome.indexOf( i[0] );
			int qat2 = chromosome.indexOf( i[1] );
			
			if( qat1 < 0 || qat2 < 0 )
				continue;
			
			if( queryMol.getBond( queryMol.getAtom(qat1), queryMol.getAtom(qat2) ) != null ) {
				//System.out.print( " | " + b2.getBeginAtom().getSymbol() + b2.getBondTypeAsDouble() + b2.getEndAtom().getSymbol() + "(" + qat1 + " " + qat2 + ")" );
				count++;
			}
		}
		//System.out.println( "bonds: " + count );
		
		//hsFoundBonds.clear();
		//hfb.clear();
		//hsFoundBonds = null;
		//hfb = null;
		
		// Now see which bonds map correctly to the query
		// Thus find the index of each atom per bond, then see if any of these are actually bonds
		// Then check for bond types
		
		
		return (double) (count / (double) queryMol.getBondCount()) + 0.00001;
		//return count;
		//return Math.random() - 0.8;
	}
	
	
	@Override
	public double fitness(ArrayList<Integer> chromosome) {
		// TODO Auto-generated method stub
		return fitnessAllBonds(chromosome);
	}
	
	/**
	 * @param chromosomes
	 * @param cindex
	 * @param probability	probability that a mutation is performed in a chromosome
	 * @param upperBound	the maximum index
	 * 
	 * Instead of changing random numbers, we select only atoms of the same type in the hyperstructure
	 */
	public void mutate( ArrayList<ArrayList<Integer>> chromosomes, int cindex, double probability, double upperBound  ) {
		ArrayList<Integer> c2 = new ArrayList<Integer>(chromosomes.get(cindex));
		
		for( int n = 0; n < c2.size(); n++ ) {

			if( Math.random() < probability ) {
				c2.set(n, (int) Math.round( (Math.random() * upperBound)) );
			}
		}

		chromosomes.add(c2);
	}
	
	public int mutationType( int queryIndex ) {
		
		int queryAtNum = queryMol.getAtom( queryIndex ).getAtomicNumber();  // atomic number of query
		
		int hIndex = (random.nextInt( hAtomsMap.get(queryAtNum).size() ) );
		
		
		//c2.set(queryIndex, atomsMap.get( queryAtNum ).get(hIndex) );
		
		int newValue = hAtomsMap.get( queryAtNum ).get(hIndex);
		
		if( newValue < 0 )
			newValue = ( -1 * random.nextInt(99) ) - 1;  
		
		//System.out.println( queryAtNum + " | " + hAtomsMap.get( queryAtNum ).size() + " | " + hIndex + " | " + newValue );
		
		return newValue;
	}
	
	/**
	 * Using the already-defined 2 molecules as hyperstructure and query, use the supplied dMCES
	 * to append query atoms to hyperstructure.
	 * 
	 * Look at all unmapped atoms in the dMCES chromosome and identify the corresponding neighbours 
	 * in the hyperstructure (translating of course, from the unmapped atom's neighbours in the query).
	 * Add a bond to each of these
	 * 
	 * To add any other bonds inferred from the dMCES, we go through all the atoms in the dMCES.
	 * For each atom, obtain the neighbouring atoms indices from the query structure, and also translate
	 * these into the corresponding hyperstructure indices using the dMCES to yield 2 lists.
	 * The union of these 2 lists tells us which bonds should now be present (and not all will be, so 
	 * ensure these "new" bonds are added!)
	 * 
	 * @param dmces
	 * @return
	 */
	public IAtomContainer createHyperstructure( ArrayList<Integer> dmces ) {

		
		/* 
		 * Construct a hash which translates the "old" index from the chromosome to the newly assigned
		 * index in the hyperstructure, for an unmapped atom that's been added to the hyperstructure
		 */
		HashMap<Integer, Integer> usedAtoms = new HashMap<Integer, Integer>(); 
		IAtomContainer novelHs = null;
		novelHs = hs;
		/*try {
			novelHs = hs.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
		// go through all (unmapped) atoms in chromosome
		for( int n = 0; n < dmces.size(); n++ ) {
			if( dmces.get( n ) < 0 ) {
				novelHs.addAtom( queryMol.getAtom( n ) );
				int index = novelHs.getAtomNumber( queryMol.getAtom( n ) );
				usedAtoms.put( dmces.get( n ), index );
				
				int[] qNeighbours = getAtomNeighboursIndices(queryMol, n); 
				//long[] hsNeighbours = new long[ qNeighbours.length ];
				
				for( int h = 0; h < qNeighbours.length; h++ ) {
					
					// identify corresponding neighbour(s) in hyperstructure
					int hsNeighbour = dmces.get( qNeighbours[h] );
					
					// get bond order of the original bond in the query molecule
					Order neighbourOrder = queryMol.getBond( 
						queryMol.getAtom(n), 
						queryMol.getAtom(qNeighbours[h]) 
					).getOrder();
					
					if( hsNeighbour >= 0 ) {
						hs.addBond( index, hsNeighbour, neighbourOrder );
					} else if( hsNeighbour < 0 && usedAtoms.containsKey(hsNeighbour) ) {
						// this' for sets of previously unmapped atoms linked to each other
						int missingAtom = usedAtoms.get(hsNeighbour);  
						novelHs.addBond( index, missingAtom, neighbourOrder );
					}
				}
			}
			
			/*
			long[] neighbours = getAtomNeighboursIndices(queryMol, n);
			System.out.println( "neighbours of " + n + " is " + Arrays.toString(neighbours) );
			
			// search for atoms in the mapping with an unmapped neighbour (in the query)
			for( int i = 0; i < neighbours.length; i++ ) {
				if( dmces.get( (int) neighbours[i] ) < 0 ) {
					long index = 0;
					
					// A measure to avoid adding a duplicate atom
					if( ! usedAtoms.containsKey(neighbours[i]) ) {
						index = hs.addAtom( queryMol.getAtomWithIdx( neighbours[i] ) );
						usedAtoms.put( neighbours[i] , index );
					} else {
						index = usedAtoms.get( neighbours[i] );
					}
					
					// identify the atom in the hyperstructure to attach the new one to
					long hsAtomToConnect = dmces.get(n);
					if( hsAtomToConnect < 0 ) {
						hsAtomToConnect = index;
					}
					
					System.out.println( index + "-" + n + " of " + hsAtomToConnect );  // n & dmces.get(n) are what we must bond the new atom to
					
					if( hs.getBondBetweenAtoms( hsAtomToConnect, n ) == null )
						hs.addBond(hsAtomToConnect, n, BondType.SINGLE);
					
				}
			}
			*/
			
		}
		
		/*
		 * Here we add any other "new" bonds inferred from the dMCES which do not necessarily involve
		 * the unmapped atoms
		 */
		for( int n = 0; n < dmces.size(); n++ ) {
			int hsIndex = dmces.get(n);  // 7 in query, 4 in hs
			if( hsIndex < 0 ) {
				hsIndex = usedAtoms.get(hsIndex);
			}
			
			int[] hsNeighbours = getAtomNeighboursIndices(novelHs, hsIndex);
			int[] qsNeighbours = getAtomNeighboursIndices(queryMol, n);
			//long[] otherHsNeighbours = new long[ qsNeighbours.length ];
			Arrays.sort( hsNeighbours );
			
			for( int a = 0; a < qsNeighbours.length; a++ ) {
				int translatedNeighbour = dmces.get( qsNeighbours[a] );
				
				//System.out.println( n + " " + translatedNeighbour + ": " + Arrays.toString(hsNeighbours) + Arrays.toString(qsNeighbours) );
				
				// find neighbours which exist in the query atom's neighbourhood but not in the hyperstructure's neighbourhood
				if( translatedNeighbour > -1 && Arrays.binarySearch( hsNeighbours , translatedNeighbour ) < 0 ) {
					
					Order neighbourOrder = queryMol.getBond( 
						queryMol.getAtom(n), 
						queryMol.getAtom(qsNeighbours[a]) 
					).getOrder();
					
					novelHs.addBond( dmces.get(n), translatedNeighbour, neighbourOrder );
				}
			}
			
		}
		
		/*
		do {
			for( int n = 0; n < dmces.size(); n++ ) {
				if( dmces.get(n) < 0 ) {
					long[] neighbours = getAtomNeighboursIndices(queryMol, n);
					System.out.println( "neighbours of " + n + " is " + Arrays.toString(neighbours) );
					
					allMapped = false;
				}
			}
			allMapped = true;
			break;  // just in-case it fails to break
		} while( ! allMapped );
		*/
		
		return novelHs;
	}
	
	
	
	protected IAtomContainer hs;
	protected IAtomContainer queryMol;
	private HashMap<Integer, ArrayList<Integer>> hAtomsMap, qAtomsMap;  // atom indices of atomic numbers in hs and query
	private HashMap<Integer, List<IBond>> hAtomBonds, qAtomBonds;  // atom indices corresponding to its neighbouring bonds - computational time saver
	private AllRingsFinder ringFinder;
	private int queryRings, hsRings;
	protected Random random;
	private SmilesGenerator sg;
	
}
