package test.java.org.cisrg.test;

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.ExtendedAlgorithm;
import org.cisrg.mapping.SimilarityComparator;
import org.cisrg.executable.MCSCommandLine;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import junit.framework.TestCase;



public class SimilarityComparatorTest extends TestCase {

	// Test all the MCSs shown in my thesis, to get the correct sizes and types (test for degree sequence orders of MCSs)
	
	
	public SimilarityComparatorTest() {
	    super( );
	    
	    algorithm = ExtendedAlgorithm.Depolli_dMCES;
	    
	    simHub = new SimilarityComparator(
				null, bondWeightFlag, ghostSubstructures, algorithm, 
				raymondHeuristics, ringHeuristics, topoDistLim, timeLimit, false
		);
	}
	
	
	public static List<Integer> getSortedDegreeListFromBondMap( IAtomContainer cp1, IAtomContainer cp2, Map<IBond, IBond> bondMap ) {
		ArrayList<Integer> degrees = new ArrayList<Integer>(cp1.getAtomCount());
		IAtomContainer sGraph = ConvenienceTools.createCommonSubgraph( cp1, cp2, bondMap );
		
		// I'm using degree sequences but TBH One could probably just compare degree frequencies instead - avoids costly sorting
		// must be degrees of atoms, not bonds
		for( IAtom b : sGraph.atoms() ) {
			degrees.add( sGraph.getConnectedBondsCount( b ) );
		}
		Collections.sort(degrees);
		
		return(degrees);
	}
	
	
	private void molGenericTest( String molFilePath, List<Integer> refFragmentSizes, List<Integer> refDegreeSeq , int refMCSSize ) {

        ArrayList<IAtomContainer> compounds = null;
        //System.out.println("Working Directory = " + System.getProperty("user.dir"));
		try {
			compounds = ConvenienceTools.getQueryMolecules( new File( molFilePath ), null, true );
			
			
			//System.out.println( "compound file size - " + compounds2.size() );
		} catch (IOException e) {
			e.printStackTrace();
		}  
		
		
		
		simHub.calculateSimilarity( compounds.get(0), compounds.get(1) );
		
		List<Integer> fragmentSizes = Arrays.stream(simHub.fragmentSizes).boxed().toList();
		int mcsSize = simHub.bondMaps.get(0).size();
		List<Integer> degrees = getSortedDegreeListFromBondMap( compounds.get(0), compounds.get(1), simHub.bondMaps.get(0) );
		
		// sort fragment sizes
		List<Integer> fragmentSizesSorted = new ArrayList<Integer>(fragmentSizes);
		Collections.sort(fragmentSizesSorted);
		
		System.out.println( "reference vs db fragment sizes: " + refFragmentSizes + " " + fragmentSizesSorted );
		System.out.println( "mcs Size (bonds):" + mcsSize );
		System.out.println( "dbMol degree sequence: " + degrees  );
		System.out.println( "MCS time (ms): " + simHub.mcsExecTime  );
		System.out.println( "MCS SMARTS: " + simHub.mcsSMARTS  );
		System.out.println( "Tanimoto similarity (from MCS): " + simHub.tanimoto  );
		System.out.println( "\n" );
		
		//String outputImagePath = "./test_mols.png" ;
		//MCSCommandLine.depictMolecules( compounds.get(0), compounds.get(1), 0, 0, outputImagePath, simHub ); 
	
		
		assertIterableEquals(refFragmentSizes, fragmentSizesSorted, "Fragment sizes not equal nor in order to reference");
		assertIterableEquals(refDegreeSeq, degrees, "Degree sequence not equal nor in order to reference");
		//assertEquals( refMCSSize, mcsSize, "MCS size does not equal reference" );
		Assertions.assertEquals(refMCSSize, mcsSize, "MCS size does not equal reference");
		
		
		
		
	}
	
	@Test
	public void testSymmetric1() {
		// get molecules
		String inputFilePath = dataPath + "/symmetric/Raymond_01.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{16, 16} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3} );
		int refMCSSize = 32;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testSymmetric2() {
		// get molecules
		String inputFilePath = dataPath + "/symmetric/Raymond_02.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{2, 3, 3, 3, 14, 14} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3} );
		int refMCSSize = 35;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testSymmetric3() {
		// get molecules
		String inputFilePath = dataPath + "/symmetric/Ersmark_01.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{34} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3} );
		int refMCSSize = 33;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testSymmetric4() {
		// get molecules
		String inputFilePath = dataPath + "/symmetric/Ersmark_02.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{2, 6, 7, 20} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3} );
		int refMCSSize = 33;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	
	@Test
	public void testSymmetric5() {
		// get molecules
		String inputFilePath = dataPath + "/symmetric/Bone_Babine.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{6, 6, 32} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3} );
		int refMCSSize = 45;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testSymmetric6() {
		// get molecules
		String inputFilePath = dataPath + "/symmetric/Libby_01.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{2, 2, 46} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3} );
		int refMCSSize = 54;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testNoRings1() {
		// get molecules
		String inputFilePath = dataPath + "/other_mcs/norings_vs_CHEMBL243.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{3, 3, 4, 4, 4, 6, 6} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3} );
		int refMCSSize = 23;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testNoRings2() {
		// get molecules
		String inputFilePath = dataPath + "/other_mcs/norings_vs_CHEMBL5373.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{2, 3, 3, 5, 6} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2} );
		int refMCSSize = 14;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testNoRings3() {
		// get molecules
		String inputFilePath = dataPath + "/other_mcs/norings_vs_CHEMBL2094108.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{2, 3, 4, 4, 4} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2} );
		int refMCSSize = 12;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testNonPlanar1() {
		// get molecules
		String inputFilePath = dataPath + "/other_mcs/non_planar_1.sdf";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{2, 3, 3, 8} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3} );
		int refMCSSize = 12;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testNonPlanar2() {
		// get molecules
		String inputFilePath = dataPath + "/other_mcs/non_planar_2.sdf";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{16} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4} );
		int refMCSSize = 16;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	
	@Test
	public void testManyRingsMapping() {
		// get molecules
		String inputFilePath = dataPath + "/other_mcs/graphene_ring_mapping.smi";
		List<Integer> refFragmentSizes = Arrays.asList( new Integer[]{22} );
		List<Integer> refDegreeSeq = Arrays.asList( new Integer[]{1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3} );
		int refMCSSize = 25;
		
		molGenericTest(inputFilePath, refFragmentSizes, refDegreeSeq, refMCSSize);
	}
	 
	
	private String dataPath = "./data/input" ;
	
	private ExtendedAlgorithm algorithm ;
	
	boolean bondWeightFlag = false;
	boolean bondTypeFlag = true;
	boolean matchRingFlag = true;
	boolean matchAtomTypeFlag = true;
	boolean ghostSubstructures = false;
	boolean raymondHeuristics = true;
	boolean ringHeuristics = true;
	
	int topoDistLim = 10;
	int timeLimit = 20000;
	
	private SimilarityComparator simHub ;
}
