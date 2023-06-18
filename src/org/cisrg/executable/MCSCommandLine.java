/**
 * 
 */
package org.cisrg.executable;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.ExtendedAlgorithm;
import org.cisrg.mapping.SimilarityComparator;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

/**
 * @author edmund
 *
 *
 * Takes two smiles strings, a ref and a database molecule, and outputs MCS in smiles
 * 
 * Output MCS as image of both molecules (highlighted, CDK style)
 */


@Command(name = "checksum", mixinStandardHelpOptions = true, version = "checksum 4.0",
         description = "Prints the checksum (SHA-256 by default) of a file to STDOUT.")
public class MCSCommandLine implements Callable<Integer> {
	
	@Parameters(index = "0", description = "Reference Molecule(s) in .smi (SMILES) or .sdf format")
    private File refFile;
	
	@Parameters(index = "1", description = "Database Molecule(s) in .smi or .sdf format")
    private File dbFile;

    @Option(names = {"-a", "--algorithmName"}, description = "Depolli_dMCES, ...")
    private String algorithmName = "Depolli_dMCES";
    
    @Option(names = {"-i", "--outputImagePath"}, defaultValue = Option.NULL_VALUE, description = "outputs a (PNG) image to the specified file path")
    private String outputImagePath;
    
    @Option(names = {"-d", "--topologicalDistanceLimit"}, defaultValue = "2", description = "(Integer) maximum distance between 2 bond pairs (between ref and db molecules) to consider when buildng MCS.  Larger allows for less geometrically constrained results, but increases time requirements.  Defaults to 2 if not set.  To disable, set to -1 or a number larger than the number of bonds in your largest input molecule. ")
    private int topologicalDistanceLimit;
    
    @Option(names = {"-r", "--ringHeuristics"}, defaultValue = "false", description = "(true/false) Enable to impose constraints on ring matching to hasten MCS calculations.  May reduce MCS size in some comparisons. ")
    private boolean ringHeuristics = false;
    
    @Option(names = {"-R", "--RaymondHeuristics"}, defaultValue = "false", description = "(true/false) Enable to impose constraints from John Raymond's RASCAL publication to hasten MCS calculations.  May reduce MCS size in some comparisons. ")
    private boolean raymondHeuristics = false;
    
    @Option(names = {"-t", "--timeLimit"}, defaultValue = "10000", description = "(Integer) Maximum time allowed per MCS calculation.")
    private int timeLimit;


    
    
    private ArrayList<IAtomContainer> setupMolecules( String inputFileName ) {

    	// get molecules
        ArrayList<IAtomContainer> compounds = null;
		try {
			compounds = ConvenienceTools.getQueryMolecules( new File( inputFileName ), null, true );
			//System.out.println( "compound file size - " + compounds.size() );
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return compounds;
    }
    
    public static void depictMolecules( IAtomContainer refMol, IAtomContainer dbMol, int rIndex, int dIndex, String outImgPath, SimilarityComparator sh ) {
    	
    	String fileSuffix = "_r" + rIndex + "_d" + dIndex + ".png";
    	String nonExtension = outImgPath.substring(0,  outImgPath.lastIndexOf(".") );
    	String currentFilePath = nonExtension + fileSuffix ;
    	
    	
		DepictionGenerator dg = new DepictionGenerator().withSize(600, 400);
		
		ArrayList<IAtomContainer> molPair = new ArrayList<IAtomContainer>(2);
		molPair.add(refMol);
		molPair.add(dbMol);
		
		dg = dg.withHighlight( sh.bondMaps.get(0).keySet() , Color.RED).withOuterGlowHighlight(2.0)
			   .withHighlight( sh.bondMaps.get(0).values() , Color.RED).withOuterGlowHighlight(2.0) ;
				
		try {
			dg.depict(molPair).writeTo( currentFilePath );
			
			System.out.println("Reference and Database molecule MCS written to " + currentFilePath );
		} catch (IOException | CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

    @Override
    public Integer call() throws Exception { // your business logic goes here...
        
    	BufferedReader reader;
    	
    	ArrayList<IAtomContainer> refMols = setupMolecules( refFile.getAbsolutePath() );
    	ArrayList<IAtomContainer> dbMols = setupMolecules( dbFile.getAbsolutePath() );

    	ExtendedAlgorithm algorithm = ExtendedAlgorithm.valueOf(algorithmName);
    	
    	
    	boolean bondWeightFlag = false;  // assign weights to matched bonds in MCS
		boolean bondTypeFlag = true;  // set True to make bond labels match, otherwise bond labels are ignored
		boolean ghostSubstructures = false;  // not currently used
		boolean detailedGhostInfo = false;  // not currently used

		  // set true to enable ring matching constraints, also to reduce search space
		//int topologicalDistanceLimit = 4;  // see CLI instructions -d option for info
		//int timeLimit = 10000;  // milliseconds
		//boolean matchRingFlag = true;  /
		//boolean matchAtomTypeFlag = true;
		
		//IAtomContainer refMol = refMols.get(0);
		//IAtomContainer dbMol = dbMols.get(1);
		
    	/*
		ExtendedIsomorphism exI = new ExtendedIsomorphism(
					refMol, dbMol, algorithm, bondTypeFlag, matchRingFlag, 
					matchAtomTypeFlag, raymondHeuristics, ringHeuristics, topoDistLim, timeLimit
		);
		*/
		
		SimilarityComparator simHub = new SimilarityComparator(
				null, bondWeightFlag, ghostSubstructures, algorithm, 
				raymondHeuristics, ringHeuristics, topologicalDistanceLimit, timeLimit, detailedGhostInfo
		);
		
		
		for( int r = 0; r < refMols.size(); r++ ) {
			IAtomContainer refMol = refMols.get(r);
			
			for( int d = 0; d < dbMols.size(); d++ ) {
				IAtomContainer dbMol = dbMols.get(d);
				
				simHub.calculateSimilarity(refMol, dbMol );

				System.out.println( simHub.mcsSMARTS );
				System.out.println( "MCS calculation time (ms): " + simHub.mcsExecTime );
				

				if( outputImagePath != null )
					depictMolecules( refMol, dbMol, r, d, outputImagePath, simHub );

				
				System.out.println( "" );  // final new line to delineate pairs
			}
		}
		
		 
		
    	return 0;
    }

	/**
	 * @param args
	 */
	// this example implements Callable, so parsing, error handling and handling user
    // requests for usage help or version help can be done with one line of code.
    public static void main(String... args) {
        int exitCode = new CommandLine(new MCSCommandLine()).execute(args);
        System.exit(exitCode);
    }

}
