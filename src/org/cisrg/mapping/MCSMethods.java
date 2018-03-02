package org.cisrg.mapping;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.cisrg.ambit.SmartsHelper;
import org.cisrg.mapping.hyperstructures.DefaultHyperstructureAtomMatcher;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.cdk.smsd.algorithm.matchers.DefaultMCSPlusAtomMatcher;

//import cern.colt.matrix.impl.DenseDoubleMatrix2D;


/**
 * 
 * @author Edmund Duesbury
 * @date October 2014
 * 
 */

public abstract class MCSMethods {

	
	protected MCSMethods() {
		sGenerator = new SmilesGenerator().aromatic();
		sParser = new SmilesParser( DefaultChemObjectBuilder.getInstance() );
		smaH = new SmartsHelper( DefaultChemObjectBuilder.getInstance() );
	}

	
	public List<List<Integer>> getBestAtomIndexMatches() {
		return mcsAtomIsomorphisms;
	}
	
	
	
	
	public List<Map<IBond, IBond>> getBestBondMatches() {
		return mcsBondIsomorphisms;
	}
	
	
	/**
	 * 
	 * @return a list of target-to-query 2-element arrays (molecule bond indices)
	 */
	public List<List<int[]>> getBestBondIndexMatches() {
		return mcsBondIndexIsomorphisms;
	}
	


	public abstract void search( IAtomContainer graph1, IAtomContainer graph2 ) throws CDKException;

	
	
	
	public int execute() {
		
		
		
		mcsBondIsomorphisms = new ArrayList< Map<IBond, IBond> >();
		mcsAtomIsomorphisms = new ArrayList<List<Integer>>();
		mcsBondIndexIsomorphisms = new ArrayList<List<int[]>>();
		
		modProdTime = 0;
		modProdEdgeDensity = 0.0;
		
		if( verbose ) {
			System.out.println( "Starting MCS search using " + this.getClass() );
		}
		
		mcsStartTime = System.currentTimeMillis();
		
		try {
			search( hsMol, qMol );
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		mcsExecTime = (System.currentTimeMillis() - mcsStartTime);
		
		
		
		if( verbose ) {
			System.out.println( "mcs search time - " + mcsExecTime );
			System.out.println( "main mol atoms - " + hsMol.getAtomCount() );
		}
		
		if( mcsBondIsomorphisms.size() > 0 ) { 
			Map<IBond, IBond> commonSubgraphMap = mcsBondIsomorphisms.get(0);
			mcsSize = commonSubgraphMap.size();
			
			if( commonSubgraph == null ) {
				//IAtomContainer commonSubgraph = ConvenienceTools.createCommonSubgraphDeepCopy(hsMol, qMol, commonSubgraphMap);
				commonSubgraph = ConvenienceTools.createCommonSubgraph(hsMol, qMol, commonSubgraphMap);
			}
			
			IAtomContainerSet fragments = ConvenienceTools.partitionIntoMolecules(commonSubgraph);
			//IAtomContainerSet fragments = ConnectivityChecker.partitionIntoMolecules(commonSubgraph);
			if( fragments.isEmpty() ) {
				fragmentSizes = new int[]{0};
			} else {
				fragmentSizes = new int[ fragments.getAtomContainerCount() ];
				
				for( int n = 0; n < fragmentSizes.length; n++ ) {
					fragmentSizes[n] = fragments.getAtomContainer(n).getAtomCount();
				}
			}
			
			// remove evil listeners
			for( IAtom at : hsMol.atoms() ) {
	        	for( IAtomContainer ac : fragments.atomContainers() )
	        		at.removeListener(ac);
	        }
			
			
			try {
				//ConvenienceTools.initializeMolecule(commonSubgraph);
				
				//if( commonSubgraph instanceof IQueryAtomContainer )
				//	mcsSMARTS = smaH.toSmarts( (QueryAtomContainer) commonSubgraph );
				//else
				if( commonSubgraph.isEmpty() ) {
					mcsSMARTS = "";
				} else {
					if( commonSubgraph instanceof IQueryAtomContainer ) {
						mcsSMARTS = smaH.toSmarts( (QueryAtomContainer) commonSubgraph) ;
					} else {
						//ConvenienceTools.calculateImplicitHydrogens( commonSubgraph );  // yeah - unfortunately SMILES generator errors on this sometimes
						mcsSMARTS = sGenerator.create( commonSubgraph );
					}
				}
					//System.out.println( "MCS = " + mcsSMARTS + " " + commonSubgraph.getBondCount() + " " + commonSubgraph.getAtomCount() );
					//System.out.println( "mapping chromosome - " + mcsAtomIsomorphisms.get(0)  );
			} catch (CDKException e1) {
				e1.printStackTrace();
			}
		} else {  // handle nulls
			mcsBondIsomorphisms.add(new HashMap<IBond, IBond>(0) );
			mcsBondIndexIsomorphisms.add( new ArrayList<int[]>(0) );
			mcsAtomIsomorphisms.add(new ArrayList<Integer>(0) );
			fragmentSizes = new int[1];
			mcsSize = 0;
			mcsSMARTS = "";
		}
		
		
		return 0;
	}
	
	
	/**
	 * Converts atom map to (obsolete) Genetic Algorithm chromosome
	 * 
	 * @param commonAtomMap
	 * @return
	 */
	/*protected List<Integer> atomMapToChromosome( Map<IAtom, IAtom> commonAtomMap ) {
		
		// query to hyperstructure mapping
		ArrayList<Integer> mapping = new ArrayList<Integer>(qMol.getAtomCount());  

		for( int n = 0; n < qMol.getAtomCount(); n++ ) {
			mapping.add( (n * -1) - 1 );  // dummy assignments
		}

		for( Entry<IAtom, IAtom> e : commonAtomMap.entrySet() ) {
			int hsAtom1Index = hsMol.getAtomNumber( e.getValue() );
			int qAtom1Index = qMol.getAtomNumber( e.getKey() );

			mapping.set( qAtom1Index, hsAtom1Index );
		}

		return mapping;
	}*/
	
	
	
	/**
	 * Maps common atoms from a bond map
	 * 
	 * for a given atom, the corresponding atom that is the most frequent in the bonds it is found in, is chosen
	 * 
	 * // FIXME  not optimal, doesn't map all atoms
	 * 
	 * @param bondMap
	 */
	protected Map<IAtom, IAtom> createAtomMapFromBondMap( IAtomContainer mol1, IAtomContainer mol2, Map<IBond, IBond> bondMap ) {
		
		Map<IAtom, IAtom> atomMap = new HashMap<IAtom, IAtom>();
		
		double[][] matchMatrix = new double[ mol1.getAtomCount() ][ mol2.getAtomCount() ];
		for( double[] row : matchMatrix ) {
			for( int r = 0; r < row.length; r++ ) {
				row[r] = 0;
			}
		}
		
		Map<IAtom, Integer> indices1 = new HashMap<IAtom, Integer>();
		Map<IAtom, Integer> indices2 = new HashMap<IAtom, Integer>();
		
		for( int i = 0; i < mol1.getAtomCount(); i++ ) {
			indices1.put( mol1.getAtom(i), i );
		}
		for( int j = 0; j < mol2.getAtomCount(); j++ ) {
			indices2.put( mol2.getAtom(j), j );
		}
		
		// count frequencies
		for( Entry<IBond, IBond> e : bondMap.entrySet() ) {
			for( IAtom bond1At : e.getKey().atoms() ) {
				for( IAtom bond2At : e.getValue().atoms() ) {
					
					if( atomsMatch(mol1, bond1At, mol2, bond2At) ) {
						matchMatrix[ indices1.get(bond1At) ][ indices2.get(bond2At) ]--;
					} 
				}
			}
		}
		
		//System.out.println( new DenseDoubleMatrix2D( matchMatrix ) );
		
		/*// choose the atom with the highest frequency
		for( int r = 0; r < mol1.getAtomCount(); r++ ) {
			int[] row = matchMatrix[r];
			
			int max = 0;
			int maxIndex = 0;
			
			for( int c = 0; c < mol2.getAtomCount(); c++ ) {
				if( row[c] > max ) {
					max = row[c];
					maxIndex = c;
				}
			}
			
			if( max > 0 )
				atomMap.put( mol1.getAtom(r) , mol2.getAtom(maxIndex) );
		}*/
		
		List<int[]> bestMatches = ConvenienceTools.linearAssignment( matchMatrix );
		
		for( int n = 0, ns = bestMatches.size(); n < ns ; n++ ) {
			if( bestMatches.get(n)[1] >= 0 && matchMatrix[ n ][ bestMatches.get(n)[1] ] < 0 )
				atomMap.put( mol1.getAtom(n) , mol2.getAtom( bestMatches.get(n)[1] ) );
		}
		
		return atomMap;
	}
	

	
	// query to hyperstructure mapping
	
	public static List<Integer> atomMapToChromosome( IAtomContainer mol1, IAtomContainer mol2,  Map<IAtom, IAtom> commonAtomMap ) {
				ArrayList<Integer> mapping = new ArrayList<Integer>(mol1.getAtomCount());  
				
				for( int n = 0; n < mol1.getAtomCount(); n++ ) {
					mapping.add( (n * -1) - 1 );  // dummy assignments
				}
				
				for( Entry<IAtom, IAtom> e : commonAtomMap.entrySet() ) {
						int hsAtom1Index = mol1.getAtomNumber( e.getKey() );
						int qAtom1Index = mol2.getAtomNumber( e.getValue() );

						//mapping.set( qAtom1Index, hsAtom1Index );
						mapping.set( hsAtom1Index, qAtom1Index );
				}
				
				return mapping;
	}
	
	
	protected static boolean atomsMatch( IAtomContainer ac1, IAtom at1, IAtomContainer ac2, IAtom at2 ) {
		//Atom Matchers
        AtomMatcher atomMatcher1 = null;
        
        // What matches the first atom in the first bond? 
        if( at1 instanceof IQueryAtom )
        	atomMatcher1 = new DefaultHyperstructureAtomMatcher( (IQueryAtom) at1, (IQueryAtomContainer) ac1 );
        else
        	atomMatcher1 = new DefaultMCSPlusAtomMatcher(ac1, at1, true);
        
        return atomMatcher1.matches(ac2, at2);
	}
	
	
	public void setMainMol( IAtomContainer m1 ) {
		hsMol = m1;
	}
	
	public void setQueryMol( IAtomContainer m2 ) {
		qMol = m2;
	}
	
	
	public IAtomContainer getMainMol() {
		return hsMol;
	}
	
	
	public IAtomContainer getQueryMol() {
		return qMol;
	}
	
	public void setMatchBonds(boolean mb) {
		matchBonds = mb;
	}
	
	public long getModularProductConstructionTime() {
		return modProdTime;
	}
	
	public long getMCSSearchTime() {
		return mcsExecTime;
	}
	
	public int getModularProductNodeCount() {
		if( modProd != null )
			return modProd.getNodes().size();
		
		return 0;
	}
	
	public double getModularProductEdgeDensity() {
		return modProdEdgeDensity;
	}
	
	/**
	 * Milliseconds in which the algorithm stops after (if applicable) starting the "expansion" process
	 * 
	 * @param to
	 */
	public void setExpansionTimeLimit( int to ) {
		expansionTimeLimitMs = to;
	}



	protected List<Map<IBond, IBond>> mcsBondIsomorphisms;
	protected List<List<int[]>> mcsBondIndexIsomorphisms;
	protected List<List<Integer>> mcsAtomIsomorphisms;
	protected IAtomContainer hsMol, qMol;
	
	protected boolean verbose = false;
	protected boolean matchBonds = true;
	protected IAtomContainer commonSubgraph = null;
	public int mcsSize = -1;
	public int[] fragmentSizes = null;
	public String mcsSMARTS = null; 
	
	/**
	 * terminate algorithm after this time, during an "expansion stage"
	 * 
	 * Expansion would be traversing the search tree in a backtracking algorithm (e.g. expanding a clique).  
	 * Could also just be extension of MCS (like in kCombu build-up stage)
	 * 
	 * Some algorithms will not use a time limit due to a lack of apparent "expansion stage" (like spectral algorithms)
	 */
	protected long expansionTimeLimitMs = 10000;  
	protected long mcsStartTime;
	protected long mcsExecTime;
	protected int numberOfSteps;
	protected long modProdTime;
	protected double modProdEdgeDensity;
	
	/* XXX
	 * It'd make more sense to put the modular product stuff in the clique detection things only, but
	 * I wanted the modular product "statistics" methods to be accessible from this class.
	 * 
	 * Also, a modular product isn't necessarily used just for clique detection
	 */
	protected GenerateCompatibilityGraphEdges modProd;
	
	protected SmilesParser sParser = null;
	protected SmilesGenerator sGenerator = null;
	protected SmartsHelper smaH = null;
	
	
	/**
	 * method to implement if you are for any reason in need of calculating a "fitness value" for the mapping.
	 * 
	 * Most useful when working with genetic algorithms or any other random search procedure.
	 * 
	 * @param plugins
	 */
	//public void setFitnessMeasure( GAPlugins plugins );
	
	//public GAPlugins getFitnessMeasure();
	//public GAChromosome getBestMatch();
}

