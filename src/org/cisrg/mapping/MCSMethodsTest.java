package org.cisrg.mapping;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;



import org.apache.xmlbeans.impl.piccolo.io.FileFormatException;
import org.cisrg.ambit.SmartsParser;
import org.cisrg.knime.ExtendedAlgorithm;
import org.cisrg.knime.ExtendedIsomorphism;
import org.cisrg.mapping.ChemAxonMCS.ChemAxonMCSOptions;
import org.cisrg.mapping.CliqueDetection.ModularProductOptions;
import org.junit.Before;
import org.junit.Test;
//import org.knime.cisrg.hyperstructures.BitSetCollectionAL;
import org.cisrg.BitSetExtended;
import org.cisrg.mapping.GenerateCompatibilityGraphEdges;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;


//import ambit2.smarts.SmartsHelper;

public class MCSMethodsTest {

	@Before
	public void setUp() throws Exception {
	}
	
	
	protected boolean hasCALicense() {
		  try {
			Class<?> exists = Class.forName( "chemaxon.license.LicenseHandler", false, this.getClass().getClassLoader() )  ;
			
			if( exists.desiredAssertionStatus() ) {
				System.out.println( chemaxon.license.LicenseHandler.MCES );
				
				if( chemaxon.license.LicenseHandler.MCES != null ) {
					return true;
				}
				
			}
			//LicenseHandler lh = new LicenseHandler();
			
			return false;
		} catch (ClassNotFoundException e) {
			return false;
		}
	}

	
	public int upperBoundTest( IAtomContainer mol1, IAtomContainer mol2 ) {
		
		int upperBound = 0;
		
		GenerateCompatibilityGraphEdges modProd = null;
		try {
			modProd = new GenerateCompatibilityGraphEdges(mol1, mol2, true, true, false, -1, null, true) ;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		List<int[]> nodes = modProd.getNodes();
		
		
		// idea 1 - calculate modal average of degrees in all neighbourhood subgraphs in modular product
		// neighbourhood degree - |intersection of this neighbourhood and adjacency list of each node|
		for( int n=0; n < nodes.size(); n++ ) {
			Collection<Integer> adj = new HashSet<>( modProd.getAdjacencyList().get(n) );
			List<Integer> degrees = new ArrayList<Integer>();
			
			for( Integer neighbour : adj ) {
				Collection<Integer> adj2 = new HashSet<>( modProd.getAdjacencyList().get(neighbour) );
				adj2.retainAll(adj);
				degrees.add( adj2.size() );
				
				if( adj2.size() > upperBound )
					upperBound = adj2.size();
			}
			
			System.out.println( degrees );
		}
		System.out.println( "avg degree - " + upperBound );
		
		return upperBound;
		
	}

	
	
	
	@Test
	public void testExample1() {
		
		
		
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mos_mapping_test.sdf";
		//String inputFileName2 = compoundPath + "/mos_mapping_test20.smi";
		//String inputFileName2 = compoundPath + "/other_mcs/graphene_ring_mapping.smi";
		String inputFileName2 = compoundPath + "/Franco_92a.smi";
		//String inputFileName2 = compoundPath + "/zhu_graphs/zhu.sdf";
		//String inputFileName2 =  compoundPath + "/non_planar.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/ring_mapping_test18.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/aid466_actives_5.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mddr_norings_5.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mddr/mddr_78374_MaxMin_10.sdf";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mddr/mddr_renin_random_10.sdf";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/reaction_mapping/JW1.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/symmetric/Raymond_02.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/thesis/media/misc/5HT_reuptake_inhibitors_AndrewMadden.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/symmetric/Ersmark_02.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/symmetric/Bone_Babine.smi";
		//String inputFileName2 = compoundPath + "/symmetric/Libby_01.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/other_mcs/norings_vs_CHEMBL243.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/other_mcs/norings_vs_CHEMBL2094108.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/other_mcs/non_planar_2.sdf";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mod_prod_simple_test.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/chembl751606.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/rascal_test_dir/duff.smi";
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/ergosterol_analogues.smi";  // 1 & 4 are hard
		//String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/rascal_test_dir/smaller.smi";  // same size molecules but swap molecules to see if answer is optimal
		//String inputFileName2 = compoundPath + "/bond_freq_test1.smi";
		//String inputFileName2 = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/hyperstructure_query_test1.smi";
		//String inputFileName2 = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/aromatic_nonaromatic_rings.smi";
		//String inputFileName2 = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/single_triple_bonds.smi";
		
		SmilesGenerator sGenerator = new SmilesGenerator().aromatic() ;
		//SmilesParser sParser = new SmilesParser( DefaultChemObjectBuilder.getInstance() );
		//SmartsHelper smaH = new SmartsHelper( DefaultChemObjectBuilder.getInstance() ) ;
		SmartsParser smaP = new SmartsParser();
		
		 // get molecules
        ArrayList<IAtomContainer> compounds2 = null;
		try {
			compounds2 = ConvenienceTools.getQueryMolecules( new File( inputFileName2 ), null, false );
			System.out.println( "compound file size - " + compounds2.size() );
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		
		IAtomContainer compound1 = compounds2.get(1);
		IAtomContainer compound2 = compounds2.get(0); 
		//Collections.reverse( compounds2 );
		
		/*int[][] bala1 = new int[compound1.getBondCount()][]; 
		List<List<Integer>> bal1 = ConvenienceTools.createBondAdjacencyList(compound1) ;
		for( int n = 0; n < bala1.length; n++ ) {
			List<Integer> rowList = bal1.get(n);
			bala1[n] = new int[rowList.size()];
			
			for( int b=0 ; b<rowList.size(); b++ ) {
				bala1[n][b] = rowList.get(b);
			}
		}
		Canon.symmetry(compound1, bala1);
		*/
		
		
		
		//upperBoundTest( compound1, compound2 );
		
		try {
			System.out.println( "Compound 1 is " + sGenerator.create( compound1 ) ); 
			System.out.println( "Compound 2 is " + sGenerator.create( compound2 ) );
		} catch (CDKException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		
		long searchTime = System.currentTimeMillis();
		

		//ExtendedIsomorphism exI = new ExtendedIsomorphism (ExtendedAlgorithm.kCombu_dMCES , true);
		//ExtendedIsomorphism exI = new ExtendedIsomorphism (ExtendedAlgorithm.CP_dMCES , true);
		//ExtendedIsomorphism exI = new ExtendedIsomorphism (ExtendedAlgorithm.Depolli_dMCES , true);
		//ExtendedIsomorphism exI = new ExtendedIsomorphism (ExtendedAlgorithm.RASCAL_dMCES , true);

		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.DEFAULT  ;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.MCSPlus ;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.VFLibMCS  ;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.CDKMCS;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.fMCS;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.ChemAxon_dMCES;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.kCombu_dMCES;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.ZhuAERdMCES;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.BK_cMCES;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.CP_dMCES;
		//ExtendedAlgorithm algorithm = ExtendedAlgorithm.RASCAL_dMCES;
		ExtendedAlgorithm algorithm = ExtendedAlgorithm.Depolli_dMCES;
		
		
		boolean rHeuristics = false;
		boolean ringHeuristics = false;
		int topoDistLim = 0;
		int timeLimit = 10000;
		
		
		/*try {
			ConvenienceTools.writeDIMACSGraph( 
					new GenerateCompatibilityGraphEdges(compound1, compound2, true, rHeuristics, rHeuristics, topoDistLim, null, true) 
			);
			
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			
		}*/
		
		
		
		ConvenienceTools.countRings(compound1);
		ConvenienceTools.countRings(compound2);
		
		
		try {
			ConvenienceTools.calculateAromaticity( compound1 );
			ConvenienceTools.calculateAromaticity( compound2 );
		} catch (CDKException e) {
			e.printStackTrace();
		}
		
		
		ExtendedIsomorphism exI = new ExtendedIsomorphism(compound1, compound2, algorithm, true, false, true, rHeuristics, ringHeuristics, topoDistLim, timeLimit);
		//exI.setTimeLimit(5000);

		//exI.setMatchBonds(false);
		//exI.setTopologicalDistanceLimit(2);
		//exI.setUseRaymondHeuristics(true);
/*
		
		MCSMethods tempMapper = new DepolliCliqueDetection( new ModularProductOptions(rHeuristics, ringHeuristics, topoDistLim) );
		try {
			GenerateCompatibilityGraphEdges tempMp = new GenerateCompatibilityGraphEdges(compound1, compound2, true, rHeuristics, ringHeuristics, topoDistLim, null, true);
			GenerateCompatibilityGraphEdges tempMp2 = new GenerateCompatibilityGraphEdges(compound1, compound2, true, rHeuristics, ringHeuristics, topoDistLim, null, true);
			//List<Collection<Integer>> mpEdges = tempMapper.modProd.getAdjacencyList() ;
			List<Collection<Integer>> mpEdges = tempMp.getAdjacencyList() ;
			int nE = 0;
			if( mpEdges.get(0) instanceof BitSetCollectionAL ) {
				for( Collection<Integer> row : mpEdges ) {
					BitSetCollectionAL bsAL = (BitSetCollectionAL) row;
					OpenBitSet bs = (OpenBitSet) bsAL.getBitSet();
					
					System.out.println( bs.size() + " " + bs.cardinality() );
					
					bs.flip(0, mpEdges.size() );
					bs.clear(nE);
					
					System.out.println( bs.size() + " " + bs.cardinality() );
					++nE;
				}
			}
			
			
			
			int[][] adjList = ConvenienceTools.listOfListsToMatrix(mpEdges);
			
			HashMap<Integer, Collection<Integer>> cover1Source = new HashMap<Integer, Collection<Integer>>();
			cover1Source.clear();
			for( int n = 0; n < adjList.length; n++ ) {
				Collection<Integer> set = new HashSet<Integer>();
				for( int s : adjList[n] ) { set.add(s); }
				cover1Source.put(n, set);
			}
			
			List<Integer> cover = ZhuSpectralMCES.findCoverIsolationAlgorithm( cover1Source, adjList);
			
			cover1Source.clear();
			for( int n = 0; n < adjList.length; n++ ) {
				Collection<Integer> set = new HashSet<Integer>();
				for( int s : adjList[n] ) { set.add(s); }
				cover1Source.put(n, set);
			}
			List<Integer> cover2 = ZhuSpectralMCES.findCoverSRA(cover1Source);
			
			cover1Source.clear();
			for( int n = 0; n < adjList.length; n++ ) {
				Collection<Integer> set = new HashSet<Integer>();
				for( int s : adjList[n] ) { set.add(s); }
				cover1Source.put(n, set);
			}
			List<Integer> cover3 = ZhuSpectralMCES.findGreedyCover(cover1Source);
			
			cover1Source.clear();
			for( int n = 0; n < adjList.length; n++ ) {
				Collection<Integer> set = new HashSet<Integer>();
				for( int s : adjList[n] ) { set.add(s); }
				cover1Source.put(n, set);
			}
			List<Integer> cover4 = ZhuSpectralMCES.findRandomCover(cover1Source, adjList);
			System.out.println( cover.size() + " " + ZhuSpectralMCES.isVertexCover(adjList, cover) );
			System.out.println( cover4.size() + " " + ZhuSpectralMCES.isVertexCover(adjList, cover4) );
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		*/
		
		//try {
			//exI.init(compound1, compound2, true, false);
			searchTime = System.currentTimeMillis() - searchTime;
		//} catch (CDKException e1) {
			// TODO Auto-generated catch block
			//e1.printStackTrace();
		//}
		
		System.out.println( "mapping size = " + exI.getAllBondMaps().get(0).size() );
		System.out.println( "mp size = " + exI.getModularProductNodeCount() );
		System.out.println( "mp density = " + exI.getModularProductEdgeDensity() );

		
		
		// delta-Y exchanges must be avoided
		int maxMaps = Math.min( 100, exI.getAllBondMaps().size() );
		for( int c = 0; c < maxMaps; c++ ) {
			//Map<IBond, IBond> commonSubgraphMap = exI.getAllBondMaps().get(c);  
			Map<IBond, IBond> commonSubgraphMap = new HashMap<IBond, IBond>();
			Map<IBond, IBond> reverseSubgraphMap = new HashMap<IBond, IBond>();
			
			
			//findCliques.cliqueMax.remove( findCliques.cliqueMax.toArray(new int[0][])[32] );
			
			
			
			//findMCS.deltaYExchangeOccured( findMCS.bestCliques.get(c) );
			
			ArrayList<Integer> degrees1 = new ArrayList<Integer>( commonSubgraphMap.size() );
			ArrayList<Integer> degrees2 = new ArrayList<Integer>( commonSubgraphMap.size() );
			
			
			for( int[] match : exI.getAllBondIndexMaps().get(c) ) {
				/*IBond a = match.getKey();
				IBond b = match.getValue();*/
				IBond a = compound1.getBond( match[0] );
				IBond b = compound2.getBond( match[1] );
				
				commonSubgraphMap.put(a, b);
				reverseSubgraphMap.put(b, a);
			}
			
			System.out.println( c + " mcs atom size - " + exI.getAllAtomMapping().get(c).getCount() );
			System.out.println( c + " mcs bond size - " + exI.getAllBondMaps().get(c).size() );
			System.out.println( c + " mcs atom size - " + exI.getAllAtomIndexMapping().get(c).size() );
			System.out.println( c + " mcs bond size - " + commonSubgraphMap.size() );
			
			IAtomContainer sGraph = ConvenienceTools.createCommonSubgraph( compound1, compound2, commonSubgraphMap );
			IAtomContainer sGraph2 = ConvenienceTools.createCommonSubgraph( compound2, compound1, reverseSubgraphMap );
			
			
			// I'm using degree sequences but TBH One could probably just compare degree frequencies instead - avoids costly sorting
			// must be degrees of atoms, not bonds
			for( IAtom b : sGraph.atoms() ) {
				degrees1.add( sGraph.getConnectedAtomsCount( b ) );
			}
			
			for( IAtom b : sGraph2.atoms() ) {
				degrees2.add( sGraph2.getConnectedAtomsCount( b ) );
			}
			
			Collections.sort( degrees1 );
			Collections.sort( degrees2 );
			
			
			System.out.println( " 1 degrees - " + degrees1 );
			System.out.println( " 2 degrees - " + degrees2 );
			
			try {
				System.out.println(" subgraph SMILES - " + sGenerator.create(sGraph) );
				System.out.println(" subgraph SMILES 2 - " + sGenerator.create(sGraph2) );
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			/*// path distance checks (strong ring heuristic)
			List<Integer> pathDistances = new ArrayList<Integer>();
			for( Entry<Integer, Integer> e1 : exI.getAllMapping().get(c).entrySet() )  {
				for( Entry<Integer, Integer> e2 : exI.getAllMapping().get(c).entrySet() )  {
					
				}
			}*/
		}
		
		System.out.println( "mapping time = " + searchTime );
		
		
		
		if( hasCALicense() ) {
		
		// ChemAxon comparison
		long caMapperTime = System.currentTimeMillis();
		ChemAxonMCSOptions caOpts = new ChemAxonMCSOptions();
		caOpts.SMARTSHandling = true;
		caOpts.ringEnforcement = false;
		caOpts.matchBonds = true;
		caOpts.connectedMode = true;
		//caOpts.verbose = true;
		//ChemAxonMCS caMapper = new ChemAxonMCS( new MolHandler(compound1, true, false).getMolecule(), new MolHandler(compound2, true, false).getMolecule(), caOpts );
		
		
		
		
		ChemAxonMCS caMapper = new ChemAxonMCS( compound1, compound2, caOpts );
		
		
		System.out.println( "test MCS" );
		//caMapper.setMainMol(compound1);
		//caMapper.setQueryMol(compound2);
		//caMapper.getOptions().verbose = true;
		caMapper.execute();
		System.out.println( "test MCS2" );
		caMapperTime = System.currentTimeMillis() - caMapperTime;
		
		System.out.println( "ChemAxon benchmark size = " + caMapper.getBestBondMatches().get(0).size() );
		System.out.println( "ChemAxon benchmark time = " + caMapperTime );
		//System.out.println( "ChemAxon benchmark time = " + caMapper.mcsTimeTaken );
		System.out.println( "ChemAxon MCS SMILES - " + caMapper.mcsSMARTS );
		System.out.println( "ChemAxon MCS bond indices - "  );
		for( int[] pair : caMapper.getBestBondIndexMatches().get(0) ) { 
			System.out.print( "[" + pair[0] + ", " + pair[1] + "]  " );
		}
		System.out.println();
		
		
			//IAtomContainer caSg = sParser.parseSmiles( caMapper.mcsSMARTS  );
			IAtomContainer caSg = smaP.parse( caMapper.mcsSMARTS  );
			System.out.println( "ChemAxon, is subgraph of 1 - " + ConvenienceTools.isSubgraph(caSg, compound1) );
			System.out.println( "ChemAxon, is subgraph of 2 - " + ConvenienceTools.isSubgraph(caSg, compound2) );
		
		}
	}
	
	
	
	@Test
	public void testExample2() {
		
		List<int[]> nodes = new ArrayList<int[]>(6);
		
		
		for( int n = 0; n < 6; n++ ) {
		  nodes.add( new int[]{ n,n } );
		}

		List<Collection<Integer>> edges = new ArrayList<Collection<Integer>>(6);
		edges.add( new BitSetExtended<>( Arrays.asList( new Integer[]{ 1,2,3,4,5 } ) ) );
		edges.add( new BitSetExtended<>( Arrays.asList( new Integer[]{ 0, 2,3,4 } ) )  );
		edges.add( new BitSetExtended<>(Arrays.asList( new Integer[]{ 0,1, 3 } ) ) );
		edges.add( new BitSetExtended<>(Arrays.asList( new Integer[]{ 0,1,2 } ) ) );
		edges.add( new BitSetExtended<>(Arrays.asList( new Integer[]{ 0,1 } ) ) );
		edges.add( new BitSetExtended<>(Arrays.asList( new Integer[]{ 0 } ) ) );
		
		//readDIMACSFile(nodes, edges, "/home/edmund/commsys.ijs.si/~matjaz/maxclique/DIMACS/DIMACS_subset/C125.9.clq");
		//readDIMACSFile(nodes, edges, "/home/edmund/commsys.ijs.si/~matjaz/maxclique/DIMACS/DIMACS_subset/san1000.clq");
		//readDIMACSFile(nodes, edges, "/home/edmund/commsys.ijs.si/~matjaz/maxclique/DIMACS/DIMACS_subset/p_hat1000-1.clq");
		//readDIMACSFile(nodes, edges, "/opt/source_MaxCliquePara_v2.2/DIMACS_subset/tempGraph4.clq");
		
		// bitset edges
		for( int n = 0; n < edges.size(); n++ ) {
			edges.set(n, new BitSetExtended<>(edges.get(n)) );
		}
		
		
		try {
			Thread.sleep(1200);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		
		
		
		int[][] adjList = ConvenienceTools.listOfListsToMatrix(edges);
		int counter = 0;
		List<Collection<Integer>> edgesInverse = new ArrayList<Collection<Integer>>(6);
		for( Collection<Integer> edgeList : edges ) {
			BitSetExtended  <Integer> bsEL = (BitSetExtended<Integer>) edgeList;
			BitSet bs = bsEL.getBitSet();
			BitSet bsInv = (BitSet) bs.clone();
			bsInv.flip(0, adjList.length);
			bsInv.clear( counter++ );
			
			//if( bsInv.cardinality() > 0 )
			edgesInverse.add( new BitSetExtended<Integer>(bsInv) );
		}
		
		int[][] adjListInv = ConvenienceTools.listOfListsToMatrix(edgesInverse);
		
		HashMap<Integer, Collection<Integer>> cover1Source = new HashMap<Integer, Collection<Integer>>();
		cover1Source.clear();
		for( int n = 0; n < adjList.length; n++ ) {
			Collection<Integer> set = new HashSet<Integer>();
			for( int s : adjList[n] ) { set.add(s); }
			cover1Source.put(n, set);
		}
		
		
		long beforeCoverTime = System.currentTimeMillis();
		//List<Integer> cover = findCoverLocalSearch( adjList );
		List<Integer> cover = findCoverIsolationAlgorithm2( adjListInv );
		//List<Integer> cover = findCoverVSA( adjListInv );
		//List<Integer> cover = ZhuSpectralMCES.findCoverSRA( cover1Source );
		long afterCoverTime = System.currentTimeMillis() - beforeCoverTime;
		
		
		
		
		// estimate clique by removing cover
		Set<Integer> estimatedClique = new BitSetExtended<Integer>( adjList.length );
		for( int n = 0; n < adjList.length; n++ ) {
			if( !  cover.contains(n) )
				estimatedClique.add(n);
		}
		
		/*for( Collection<Integer> row : edges ) {
			row.retainAll(estimatedClique);
		}
		for( Integer e : estimatedClique ) {
			System.out.print( edges.get(e).size() + " " );
		}*/
		
		System.out.println( "cover - " + cover + " " + ZhuSpectralMCES.isVertexCover(adjListInv, cover) + " " + cover.size() );
		System.out.println( "cover algorithm time (ms) - " + afterCoverTime );
		System.out.println( "cover-derived clique - " + estimatedClique + " " + estimatedClique.size() );
		
		
		ModularProductOptions mpOpts = new ModularProductOptions(false, false, -1);
		DepolliCliqueDetection findCliques = new DepolliCliqueDetection(mpOpts);
		//findCliques.verbose = true;
		 
		findCliques.setExpansionTimeLimit(10000);
		findCliques.mcsStartTime = System.currentTimeMillis();
		findCliques.bestCliques = new ArrayList<List<Integer>>();
		findCliques.findCliquesDP( nodes, edges, null, null );
		
		
		System.out.println( "clique " + findCliques.bestCliques.get(0) + " " + findCliques.bestCliques.get(0).size() );
		
		
		
		
	}
	
	
	

	/**
	 * XXX Approximate minimum vertex cover.  Isolation Algorithm by Ugurlu 2012 - "New heuristic algorithm for unweighted minimum vertex cover"
	 * 
	 * @param workBase  A starting adjacency list in the form of a Map (in which the algorithm performs removals on)
	 * @return
	 */
		public static List<Integer> findCoverIsolationAlgorithm( int[][] adjList ) {
	        // C <-- {}
			
			List<Set<Integer>> adjMatrix = new ArrayList<>( adjList.length );
			List<Set<Integer>> adjMatrixCopy = new ArrayList<>( adjList.length );  // for the redundancy check
			
			
			long bitsetTime = System.currentTimeMillis();
			
			// adjacency list to bitset matrix
			for( int n = 0; n < adjList.length; n++ ) {
				Set<Integer> rowBs = new BitSetExtended<>(adjList.length);
				for( int elem : adjList[n] ) {
					rowBs.add( elem );
				}
				adjMatrix.add(rowBs);
				adjMatrixCopy.add( new BitSetExtended<>(rowBs) );
			}

			Set<Integer> cover = new BitSetExtended<Integer>( adjList.length );
			
			System.out.println( "Time to construct bitset matrix from adjacency list - " + (System.currentTimeMillis() - bitsetTime) );
			
			

			int origSize = adjList.length;
			
			int minDegree = origSize;
			int minDegreeAt = 0;
			
			boolean empty = false;
			
	        // while G' != {} (i.e. stop process when there're no edges left)
	        while ( ! empty ) { 

	        	// finding of vertex with minimum degree
	        	minDegree = origSize;
	        	for( int n = 0; n < adjList.length; n++ ) {
	        		Set<Integer> rowBs = adjMatrix.get(n);
	        				
	        		if( rowBs.size() < minDegree && rowBs.size() > 0 ) {
	        			minDegree = rowBs.size();
	        			minDegreeAt = n;
	        		}
	        	}
	            
	            // add all neighbours of the new vertex to the solution
	        	Set<Integer> neighbours = new BitSetExtended<>( adjMatrix.get( minDegreeAt ) );
	            cover.addAll( neighbours );
	            
	        
	            
	            
	            // remove from G' every edge incident on the neighbours (as well as original v)
	            for( int amr = 0; amr < adjMatrix.size(); amr++ ) {
	            	Set<Integer> rowBs = adjMatrix.get(amr);
	            	rowBs.removeAll( cover );
	            	rowBs.remove(minDegreeAt);
	            	
	            	if( cover.contains(amr) )
	            		rowBs.clear();
	            }
	            adjMatrix.get(minDegreeAt).clear();
	            
	            
	            // check if adjMatrix is empty (has no set things)
	            empty = true;
	            for( Set<Integer> rowBs : adjMatrix ) {
	            	if( ! rowBs.isEmpty() ) {
	            		empty = false;
	            		break;
	            	}
	            }
	        
	            
	            //System.out.println( "min degree info - " + minDegreeAt + " " + minDegree + " workbase info (neighbours) - "  + neighbours + " " + adjMatrix );
	            //System.out.println( "cover info - " + cover.size() + " " + cover );
	            
	            
	        }
	        
	       
	        //cover.add(2);
	        //cover.add(5);
	         
	        
	        System.out.println( "cover info (pre-redundancy check) - " + cover.size() + " " + cover );
	        
	        //boolean test = cover.containsAll(cover);

	        // removal of redundant vertices from cover
	        Set<Integer> keys = new BitSetExtended<Integer>( cover );
	        Collection<Integer> toRemove = new ArrayList<>();
	        for( Integer v : keys ) {
	        	boolean redundant = false;
	        	
	        	//Set<Integer> cNeighbours = 
	        	
	        	// a vertex is redundant if all its neighbouring vertices in the graph, are also in the cover
	        	if( cover.containsAll(adjMatrixCopy.get(v)) ) {
	        		redundant = true;
	        	}
	        	 
	        	if( redundant ) {
	        		toRemove.add(v);
	        		//keys = new HashSet<>(cover);
	        	}
	        	
	        	/*cover.remove(v);
	        	if( ZhuSpectralMCES.isVertexCover(adjList, cover) ) {
	        		toRemove.add(v);
	        	}
	        	cover.add(v);*/
	        }
	        
	        
	        
	        cover.removeAll(toRemove);
	        System.out.println( "cover info (post-redundancy) - " + cover.size() + " " + cover + " " + toRemove );
	        

	        return new ArrayList<Integer>( cover );
	    }
		
		
		
		
		/*
		 * trying out alternative methods of vertex selection
		 */
		public static List<Integer> findCoverIsolationAlgorithm2( int[][] adjList ) {
	        // C <-- {}
			
			List<Set<Integer>> adjMatrix = new ArrayList<>( adjList.length );
			List<Set<Integer>> adjMatrixCopy = new ArrayList<>( adjList.length );  // for the redundancy check
			
			
			long bitsetTime = System.currentTimeMillis();
			
			// adjacency list to bitset matrix
			for( int n = 0; n < adjList.length; n++ ) {
				Set<Integer> rowBs = new BitSetExtended<>(adjList.length);
				for( int elem : adjList[n] ) {
					rowBs.add( elem );
				}
				adjMatrix.add(rowBs);
				adjMatrixCopy.add( new BitSetExtended<>(rowBs) );
			}

			
			List<Set<Integer>> covers = new ArrayList<>();
			
			System.out.println( "Time to construct bitset matrix from adjacency list - " + (System.currentTimeMillis() - bitsetTime) );
			
			
			
			List<Integer> minDegreeIndices = new ArrayList<>();

        	// finding of vertex with minimum degree
        	int minDegree = adjList.length;
        	for( int n = 0; n < adjList.length; n++ ) {
        		Set<Integer> rowBs = adjMatrix.get(n);

        				
        		if( rowBs.size() < minDegree && rowBs.size() > 0 ) {
        			minDegree = rowBs.size();
        			//minDegreeAt = n;
        			minDegreeIndices.clear();
        			minDegreeIndices.add( n );
        		} else if( rowBs.size() == minDegree && rowBs.size() > 0 ) {
        			minDegreeIndices.add( n );
        		}
        	}
			

			int origSize = adjList.length;
			
			//minDegree = origSize;
			int minDegreeAt = 0;
			
			
			
			System.out.println( "minimum degrees (starting) - " + minDegreeIndices );
			
			
			for( int mda : minDegreeIndices ) {
				
				adjMatrix = new ArrayList<>( adjList.length );
				
				// adjacency list to bitset matrix
				for( int n = 0; n < adjList.length; n++ ) {
					Set<Integer> rowBs = new BitSetExtended<>(adjList.length);
					for( int elem : adjList[n] ) {
						rowBs.add( elem );
					}
					adjMatrix.add(rowBs);
				}
				
				Set<Integer> cover = new BitSetExtended<Integer>( adjList.length );
				minDegreeAt = mda;
				
				boolean empty = false;
				
		        // while G' != {} (i.e. stop process when there're no edges left)
		        while ( ! empty ) { 
		        	
		        	/*List<Integer> minDegreeIndices = new ArrayList<>();
	
		        	// finding of vertex with minimum degree
		        	minDegree = origSize;
		        	//int minSupport = adjList.length * adjList.length;
		        	for( int n = 0; n < adjList.length; n++ ) {
		        		Set<Integer> rowBs = adjMatrix.get(n);
		        		
		        		int support = 0;
		        		for( int e : rowBs ) {
		        			support += adjMatrix.get(e).size();
		        		}
		        				
		        		if( rowBs.size() < minDegree && rowBs.size() > 0 ) {
		        			minDegree = rowBs.size();
		        			//minDegreeAt = n;
		        			minDegreeIndices.clear();
		        			minDegreeIndices.add( n );
		        		} else if( rowBs.size() == minDegree && rowBs.size() > 0 ) {
		        			minDegreeIndices.add( n );
		        		}
		        		if( support < minSupport && rowBs.size() > 0 ) {
		        			minSupport = support;
		        			minDegreeAt = n;
		        		}
		        	}
		        	
		        	minDegreeAt = minDegreeIndices.get( Math.min( minDegreeIndices.size() - 1, minDegreeIndices.size() - 1 ) );*/
		        	
		        	
		            
		            // add all neighbours of the new vertex to the solution
		        	Set<Integer> neighbours = new BitSetExtended<>( adjMatrix.get( minDegreeAt ) );
		            cover.addAll( neighbours );
		            
		        
		            
		            
		            // remove from G' every edge incident on the neighbours (as well as original v)
		            for( int amr = 0; amr < adjMatrix.size(); amr++ ) {
		            	Set<Integer> rowBs = adjMatrix.get(amr);
		            	rowBs.removeAll( cover );
		            	
		            	if( cover.contains(amr) )
		            		rowBs.clear();
		            }
		            
		            
		            // check if adjMatrix is empty (has no set things)
		            empty = true;
		            for( Set<Integer> rowBs : adjMatrix ) {
		            	if( ! rowBs.isEmpty() ) {
		            		empty = false;
		            		break;
		            	}
		            }
		            
		            
		            // finding of vertex with minimum degree
		        	minDegree = origSize;
		        	//int minSupport = adjList.length * adjList.length;
		        	for( int n = 0; n < adjList.length; n++ ) {
		        		Set<Integer> rowBs = adjMatrix.get(n);
		        				
		        		if( rowBs.size() < minDegree && rowBs.size() > 0 ) {
		        			minDegree = rowBs.size();
		        			minDegreeAt = n;
		        		}
		        	}
		        
		            
		            //System.out.println( "min degree info - " + minDegreeAt + " " + minDegree + " workbase info (neighbours) - "  + neighbours + " " + adjMatrix );
		            //System.out.println( "cover info - " + cover.size() + " " + cover );
		            
		        }    
		        
		        covers.add( cover );
	        }
	        
	       
	        //cover.add(2);
	        //cover.add(5);
			
			//System.out.println( covers );
			int minCSize = covers.get(0).size();
			int minC = 0;
			for( int c = 0; c < covers.size(); c++ ) {
				Set<Integer> currentCover = covers.get(c);
				System.out.println( "cover size - " + currentCover.size() + " " + currentCover );
				
				if( currentCover.size() < minCSize ) {
					minCSize = currentCover.size();
					minC = c;
				}
			}
	         
			Set<Integer> cover = covers.get(minC);
			
	        System.out.println( "cover info (pre-redundancy check) - " + cover.size() + " " + cover );
	        
	        //boolean test = cover.containsAll(cover);

	        // removal of redundant vertices from cover
	        Set<Integer> keys = new BitSetExtended<Integer>( cover );
	        Collection<Integer> toRemove = new ArrayList<>();
	        for( Integer v : keys ) {
	        	boolean redundant = false;
	        	
	        	//Set<Integer> cNeighbours = 
	        	
	        	// a vertex is redundant if all its neighbouring vertices in the graph, are also in the cover
	        	if( cover.containsAll(adjMatrixCopy.get(v)) ) {
	        		redundant = true;
	        	}
	        	 
	        	if( redundant ) {
	        		toRemove.add(v);
	        		//keys = new HashSet<>(cover);
	        	}
	        	
	        	/*cover.remove(v);
	        	if( ZhuSpectralMCES.isVertexCover(adjList, cover) ) {
	        		toRemove.add(v);
	        	}
	        	cover.add(v);*/
	        }
	        
	        
	        
	        cover.removeAll(toRemove);
	        System.out.println( "cover info (post-redundancy) - " + cover.size() + " " + cover + " " + toRemove );
	        

	        return new ArrayList<Integer>( cover );
	    }
		
		
		
		
		
		
		/** XXX	Using greedy minimal cover rather than random minimal cover
		 * 
		 * Uses vertex support algorithm by Balaji et al (2010) "Optimization of Unweighted Minimum Vertex Cover"
		 * 
		 * @param adjList  an adjacency list of the graph
		 * @return
		 */
			public static List<Integer> findCoverVSA( int[][] adjList ) {
		        // C <-- {}
				
				 
				
				List<Set<Integer>> adjMatrix = new ArrayList<>( adjList.length );
				//List<Set<Integer>> adjMatrixCopy = new ArrayList<>( adjList.length );  // for the redundancy check
				
				
				long bitsetTime = System.currentTimeMillis();
				
				// adjacency list to bitset matrix
				for( int n = 0; n < adjList.length; n++ ) {
					Set<Integer> rowBs = new BitSetExtended<>(adjList.length);
					for( int elem : adjList[n] ) {
						rowBs.add( elem );
					}
					adjMatrix.add(rowBs);
					//adjMatrixCopy.add( new BitSetCollectionAL<>(rowBs) );
				}

				Set<Integer> cover = new BitSetExtended<Integer>( adjList.length );
				
				System.out.println( "Time to construct bitset matrix from adjacency list - " + (System.currentTimeMillis() - bitsetTime) );

				
				
				
				boolean empty = false;
				
				
		        // while G' != {}
		        while ( ! empty  ) { 
		             
		        	int maxSupport = -1;  // max support
		        	List<Integer> potentialIndices = new ArrayList<>();
		        	
		        	
		        	// list all vertices of maximum support
		        	for( int n = 0; n < adjList.length; n++ ) {
		        		Set<Integer> rowBs = adjMatrix.get(n);
		        		
		        		int support = 0;
		        		for( int e : rowBs ) {
		        			support += adjMatrix.get(e).size();
		        		}
		        		 
		        		
		        		if( support > maxSupport && rowBs.size() > 0 ) {
		        			maxSupport = support;
		        			potentialIndices.clear();
		        			potentialIndices.add(n);
		        		} else if( support == maxSupport && rowBs.size() > 0 ) {
		        			potentialIndices.add(n);
		        		}
		        	}
		        	
		        	// choose of max support vertices, the vertex with the maximum degree
		        	int maxDegreeAt = -1;
		        	int maxDegree = -1;  // max support
		        	
		        	for( Integer v : potentialIndices ) {
		        		int degree = adjMatrix.get(v).size();
		        		
		        		if( degree > maxDegree ) {
		        			maxDegree = degree;
		        			maxDegreeAt = v;
		        		}
		        	}
		        	
		         
		            
		        	
		            //System.out.println( maxDegreeAt.getProperty(idProperty) + " " + punchingBag.getBondCount() );
		            
		            // C <-- C U {v}
		            cover.add( maxDegreeAt );

		            // remove from G' every edge incident on v, and v itself
		            for( int amr = 0; amr < adjMatrix.size(); amr++ ) {
		            	Set<Integer> rowBs = adjMatrix.get(amr);
		            	rowBs.remove(maxDegreeAt);
		            	
		            }
		            adjMatrix.get(maxDegreeAt).clear();
		            
		           

		         // check if adjMatrix is empty (has no set things)
		            empty = true;
		            for( Set<Integer> rowBs : adjMatrix ) {
		            	if( ! rowBs.isEmpty() ) {
		            		empty = false;
		            		break;
		            	}
		            }
		            
		            
		          System.out.println( "max degree info - " + maxDegreeAt + " " + maxDegree + " workbase info (neighbours) - " + " " + adjMatrix );
		            System.out.println( "cover info - " + cover.size() + " " + cover );
		            
		        }
		        
		        //if( verbose )	System.out.println( "cover verification = " + isVertexCover(g, cover) );

		        return new ArrayList<Integer>( cover );
		    }
			
		
		
		public static List<Integer> findCoverLocalSearch( int[][] adjList ) {
	        // C <-- {}
			
			List<Set<Integer>> adjMatrix = new ArrayList<>( adjList.length );
			List<Set<Integer>> adjMatrixCopy = new ArrayList<>( adjList.length );  // for the redundancy check
			
			for( int n = 0; n < adjList.length; n++ ) {
				Set<Integer> rowBs = new BitSetExtended<>(adjList.length);
				for( int elem : adjList[n] ) {
					rowBs.add( elem );
				}
				adjMatrix.add(rowBs);
				adjMatrixCopy.add( new BitSetExtended<>(rowBs) );
			}


			Set<Integer> d = new BitSetExtended<Integer>( adjList.length );
			Set<Integer> cover = new BitSetExtended<Integer>( adjList.length );
			
			

			int origSize = adjList.length;
			 
			
			//int maxDegree = 0;
			int maxSupport = 0;
			int maxSupportAt = 0;
			
			boolean empty = false;
			
	        // while G' != {} (i.e. stop process when there're no edges left)
	        while ( ! empty ) { 

	        	// finding of vertex with minimum degree
	        	
	        	maxSupport = 0;
	        	for( int n = 0; n < adjList.length; n++ ) {
	        		Set<Integer> rowBs = adjMatrix.get(n);
	        				
	        		int support = rowBs.size();
	        		
	        		for( int e : rowBs ) {
	        			support += adjMatrix.get(e).size();
	        		}
	        		
	        		if( support > maxSupport ) {
	        			maxSupport = support;
	        			maxSupportAt = n;
	        		}
	        		
	        	}
	        	
	        	d.add( maxSupportAt );
	            
	        

	            
	            // remove all adjacent edges of added vertex from search space
	            for( int amr = 0; amr < adjMatrix.size(); amr++ ) {
	            	Set<Integer> rowBs = adjMatrix.get(amr);
	            	rowBs.remove( maxSupportAt );
	            }
	            adjMatrix.get( maxSupportAt ).clear();
	            
	            // check if adjMatrix is empty (has no set things)
	            empty = true;
	            for( Set<Integer> rowBs : adjMatrix ) {
	            	if( rowBs.size() > 0 ) {
	            		empty = false;
	            		break;
	            	}
	            }
	        
	            //workBase.remove( minDegreeAt );
	            
	            System.out.println( "max degree info - " + maxSupportAt + " " + maxSupport + " workbase info (neighbours) - " + " " + adjMatrix );
	            //System.out.println( "cover info - " + cover.size() + " " + cover );
	            
	            
	        }
	        
	        System.out.println( "cover (pre-redundancy check) info - " + d.size() + " " + d );
	        
	        
	        
	        
	        

	        // removal of redundant vertices from cover
	        redundancy: while( ! d.isEmpty() ) {
	        	
	        	// create adjacency matrix for the subgraph D
		        List<Set<Integer>> adjMatrixD = new ArrayList<>( adjList.length );  // for the redundancy check
				
				for( int n = 0; n < adjList.length; n++ ) {
					Set<Integer> rowBs = new BitSetExtended<>( adjMatrixCopy.get(n) );
					rowBs.retainAll( d );
					adjMatrixD.add( rowBs );
				}
	        	
	        	// find vertex with minimum support
	        	int minSupport = adjList.length * adjList.length;  // worst case value
	        	int w = 0;  // index of minimum support
	        	
	        	for( int n = 0; n < adjList.length; n++ ) {
	        		Set<Integer> rowBs = adjMatrixD.get(n);
	        		
	        		if( ! d.contains(n) )
	        			continue;
	        				
	        		int support = rowBs.size();
	        		
	        		for( int e : rowBs ) {
	        			support += adjMatrixD.get(e).size();
	        		}
	        		
	        		if( support < minSupport ) {
	        			minSupport = support;
	        			w = n;
	        		}
	        	}
	        	
	        	
	        	// check adj matrix to see if all edges adjacent to w is a subset of edges of all other vertices - implies redundancy
	        	boolean isSubset = true;
	        	
	        	Set<Integer> wAdj = adjMatrixD.get(w);
	        	for( Set<Integer> oAdj : adjMatrixD ) {
	        		if( oAdj.size() > 0 ) {
	        			if( ! oAdj.containsAll(wAdj) ) {
	        				isSubset = false;
	        				break;
	        			}
	        		}
	        	}
	        	
	        	
	        	if( isSubset ) {
	        		boolean subset2 = true;
	        		
	        		// if we remove w from d - for each v in D-w, are the edges a subset of in the original graph?
	        		for( int v : d ) {
	        			if( v == w )
	        				continue;
	        			
	        			
	        		}
	        		
	        		cover.add(w);
	        		d.remove(w);
	        		
	        	} else {
	        		cover.add(w);
	        		d.remove(w);
	        		  
	        	}
	        	
	        	System.out.println( "min support info " + w + " " + minSupport + " " + d );
	        	System.out.println( "cover info - " + cover.size() + " " + cover + " " + d );
	        	
	        }
	        
	        //cover.removeAll(toRemove);
	        System.out.println( "cover info - " + cover.size() + " " + cover );
	        

	        return new ArrayList<Integer>( cover );
	    }
		
		
		
		private void readDIMACSFile( List<int[]> nodes, List<Collection<Integer>> edges, String filePath ) {
			File graphFile = new File( filePath ); 
			 
			
			BufferedReader br = null;
			try {
				br = new BufferedReader( new FileReader(graphFile) );
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			Pattern nodeCounter = Pattern.compile("^p edge\\s+(\\d+).*");
			Pattern edgeInformationFinder = Pattern.compile("^e (\\d+) (\\d+)");
			
			nodes.clear();
			edges.clear();
			
			//List<int[]> nodes = null;
			//List<Collection<Integer>> edges = null;
			
			// read file & get info
			int nodeNumber = 1000;
			String graphLine = null;
			try {
				while( (graphLine = br.readLine()) != null ) {
					
					Matcher ncMatcher = nodeCounter.matcher(graphLine);
					Matcher eiMatcher = edgeInformationFinder.matcher(graphLine);
					
					if( ncMatcher.find() ) {
						
						String edgeNoString = ncMatcher.group(1);
						nodeNumber = Integer.parseInt(edgeNoString);
						
						
						//nodes = new ArrayList<int[]>(nodeNumber);
						//edges = new ArrayList<Collection<Integer>>(nodeNumber);
						
						for( int n = 0; n < nodeNumber; n++ ) {
							  nodes.add( new int[]{ n,n } );
							  edges.add( new ArrayList<Integer>(nodeNumber/2) );
						}
					}
					
					if( eiMatcher.find() ) {

						String firstNode = eiMatcher.group(1);
						String secondNode = eiMatcher.group(2);
						
						int node1 = Integer.parseInt(firstNode) - 1;
						int node2 = Integer.parseInt(secondNode) - 1;
						
						edges.get(node1).add(node2);
						edges.get(node2).add(node1);
					}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			
			//System.out.println( edges );
		}
		
		
		
		
		
	
	//private static String compoundPath = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/";
	private static String compoundPath = "/home/edmund/Documents/workspace/cisrg/data/input/";
	//private static String compoundPath = "/home/u054444/workspace/cisrg/data/input/";

}
