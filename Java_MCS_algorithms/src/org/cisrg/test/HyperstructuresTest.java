package org.cisrg.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.io.StringBufferInputStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.xmlbeans.impl.piccolo.io.FileFormatException;
import org.cisrg.ambit.SmartsBondExpression;
import org.cisrg.ambit.SmartsHelper;
import org.cisrg.ambit.SmartsParser;
import org.cisrg.hyperstructures.CDKSMARTSHyperstructureFitness;
import org.cisrg.mapping.ChemAxonMCS;
import org.cisrg.mapping.CliqueDetection;
import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.DepolliCliqueDetection;
import org.cisrg.mapping.FMCS;
import org.cisrg.mapping.MCSMethods;
import org.cisrg.mapping.ChemAxonMCS.ChemAxonMCSOptions;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowIterator;
import org.knime.core.data.RowKey;
import org.knime.core.node.BufferedDataContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.knime.type.CDKValue;
import org.openscience.cdk.smiles.FixBondOrdersTool;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smsd.tools.MolHandler;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

//import chemaxon.formats.MolImporter;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;
//import chemaxon.struc.Molecule;


public class HyperstructuresTest {

	@Before
	public void setUp() throws Exception {
		
		//mapper = new ChemAxonMCS( hs, compounds.get(1), false, false );
		sGenerator = new SmilesGenerator();
		sGenerator.setUseAromaticityFlag(true);
		sGenerator = sGenerator.aromatic();
    	sp = new SmilesParser( DefaultChemObjectBuilder.getInstance() );
    	smaH = new SmartsHelper( DefaultChemObjectBuilder.getInstance() );
    	smaP = new SmartsParser();
    	uit = new UniversalIsomorphismTester();
    	
    	//sp.setPreservingAromaticity(true);
	}
	
	
	

	@Test
	public void testRingEnforcementConstruction() {
		
    	
		String inputFileName = molDir + "mos_mapping_test4.smi";

 		
        // get molecules
        ArrayList<IAtomContainer> compounds = null;
		try {
			compounds = ConvenienceTools.getQueryMolecules( new File( inputFileName ), null );
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

        /*
        while( rows.hasNext() ) {
        	DataRow row = rows.next();
        	
        	if( row.getCell(molColIndex) instanceof CDKValue ) {
        		
        		IAtomContainer mol = ((CDKValue) row.getCell(molColIndex)).getAtomContainer();
        		//compounds.add( fbot.kekuliseAromaticRings( mol ) );
        		//ConvenienceTools.initializeMolecule(mol);
        		//mol = AtomContainerManipulator.removeHydrogens( mol );
        		
        		compounds.add( sp.parseSmiles( sGenerator.createSMILES(mol) ) );  // remove hydrogens - no need to match these
        		//System.out.println( "SMILES 1 " + sGenerator.createSMILES(mol) );
        	}
        }
*/
 

        
        // sort molecules
        //ConvenienceTools.sortBySimilarityMCS(compounds);
        //Collections.sort(compounds, GATest.BONDASCENDING);
	


		//String hsSmiles = sGenerator.createSMILES( compounds.get(0) );
		System.out.println( "SMILES 2 " + sGenerator.createSMILES( compounds.get(1) ) );
		
		
		IQueryAtomContainer hs = null;
		MCSMethods mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		//mcsOpts.matchAtoms = false;
		mcsOpts.ringEnforcement = true;
		mcsOpts.SMARTSHandling = false;
		mapper = new ChemAxonMCS( hs, compounds.get(1), mcsOpts );
		
		hs = createHyperstructure(mapper, inputFileName, mcsOpts, 0);
		
		
		/*
		for( int m = 0; m < compounds.size(); m++ ) {
			try {
				ConvenienceTools.initializeMolecule(compounds.get(m));
				compounds.set(m, AtomContainerManipulator.removeHydrogens( compounds.get(m) ) ); // remove hydrogens - no need to match these
				
				//queries.set(m, fbot.kekuliseAromaticRings( queries.get(m) ) );
			} catch (CDKException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		*/
		
	
		//hs = (IQueryAtomContainer) ConvenienceTools.calculateHyperstructure(compounds, mapper, 3, 3, false, false, true, false);

		String hsSMARTS = smaH.toSmarts( (QueryAtomContainer) hs);
		System.out.println( "hsSMARTS - " + hsSMARTS );
		
		
		
		
		// more complex
		hs = null;
		inputFileName = molDir + "mos_mapping_test.smi";

		// get molecules
		try {
			compounds = ConvenienceTools.getQueryMolecules(new File(
					inputFileName), null);
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		// sort molecules
		// ConvenienceTools.sortBySimilarityMCS(compounds);
		// Collections.sort(compounds, GATest.BONDASCENDING);

		mapper = new ChemAxonMCS(hs, compounds.get(1), mcsOpts);

/*		hs = (IQueryAtomContainer) ConvenienceTools.calculateHyperstructure(
				compounds, mapper, 24, 30, false, false, true, false);*/
		hs = createHyperstructure(mapper, inputFileName, mcsOpts, 0);
		ConvenienceTools.countRings(hs);
		try {
			System.out.println("mol 3: " + sGenerator.create( compounds.get(2) ) );
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		String hsSMARTS2 = smaH.toSmarts((QueryAtomContainer) hs);
		System.out.println("hsSMARTS2 - " + hsSMARTS2);
				
		
		//Assert.assertEquals("O!@;=C1!@;-N!@;-C2(!@;-C!@;-C!@;-C!@;-1)@;-C@;-C@;-C@;-C@;-C@;-2", hsSMARTS);
		Assert.assertEquals("O!@;=C2!@;-N!@;-C1(@;-C@;-C@;-C@;-C@;-C@;-1)!@;-C!@;-C!@;-C!@;-2", hsSMARTS);
		String correctHS1 = "C!@;-[#6]3c1@:c2!@;-,:N!@;-C6(!@;-O!@;-c(@:c@:[#6]1)(@:c@:2)@:[#6]3!@;-C4@;-C@;-C(!@;-C!@C)@;-C@;-C5@;-N@;-C@C@;-C@;-4@;-5)-C(@;-C@;-C@;-S~6)!@;-O";
		String correctHS2 = "C!@;-[#6,#6]4c3(:c2:c:c5(!@;-O!@;-C1(-C(@;-C@;-C@;-S~1)!@;-O)!@;-N!@;-2):c(:c:3)!@;-,:C6(!@;-[#6,#6]45)@;-C@;-C(!@;-C!@C)@;-C@;-C7@;-N@;-C@C@;-C@;-6@;-7)!@;-,:C";
		//Assert.assertEquals("C!@;-[#6]4c3:c2:c:c(!@;-O!@;-C1(-C(@;-C@;-C@;-S~1)!@;-O)!@;-N!@;-2)(:c:[#6]3)[#6]4!@;-C5@;-C@;-C(!@;-C!@C)@;-C@;-C6@;-N@;-C@C@;-C@;-5@;-6", hsSMARTS2);
		Assert.assertTrue( hsSMARTS2.equals(correctHS1) || hsSMARTS2.equals(correctHS2) );
		
		//fail("Not yet implemented");
	}
	
	
	@Test
	public void testNonRingEnforcementConstruction() {
		
		IQueryAtomContainer hs = null;
		MCSMethods mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		//mcsOpts.matchAtoms = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		//mcsOpts.verbose = true;
		
		//mapper = new ChemAxonMCS( hs, null, mcsOpts );
		mapper = setUpMapper(mcsOpts);
		
		// Simple example 
		String inputFileName = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mos_mapping_test4.smi";

        // get molecules
        //ArrayList<IAtomContainer> compounds = null;
		

        // sort molecules
        //ConvenienceTools.sortBySimilarityMCS(compounds);
        //Collections.sort(compounds, GATest.BONDASCENDING);
	
		
		
		//mapper.execute();
		System.out.println( "MCS String1 : " + mapper.mcsSMARTS );
		
		//hs = (IQueryAtomContainer) ConvenienceTools.calculateHyperstructure(compounds, mapper, 3, 3, false, false, false, false);
		
		hs = createHyperstructure(mapper, inputFileName, mcsOpts, 2);
		
		
		ConvenienceTools.countRings( hs );
		//hsf.setRingBonds();
		System.out.println( "MCS String: " + mapper.mcsSMARTS );
		System.out.println( "ringed: " + smaH.toSmarts( (QueryAtomContainer) hs ) );

		String hsSMARTS1 = smaH.toSmarts( (QueryAtomContainer) hs);
		System.out.println( "hsSMARTS example 1 - " + hsSMARTS1 );
		
		
		// more complex
		inputFileName = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/mos_mapping_test.smi";

		hs = (IQueryAtomContainer) createHyperstructure(mapper, inputFileName, mcsOpts, 3);
		ConvenienceTools.countRings( hs );
		System.out.println( "ringed: " + smaH.toSmarts( (QueryAtomContainer) hs ) );

		String hsSMARTS2 = smaH.toSmarts( (QueryAtomContainer) hs);
		System.out.println( "hsSMARTS2 - " + hsSMARTS2 );
		
		
		try {
			System.out.println( uit.isIsomorph( sp.parseSmiles("O=C2NC1(CCCCC1)-C-C-C-2"), smaP.parse(hsSMARTS1)  ) );
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		// more complex
				inputFileName = molDir + "chembl751606.smi";

				hs = (IQueryAtomContainer) createHyperstructure(mapper, inputFileName, mcsOpts, 5);
				ConvenienceTools.countRings( hs );
				System.out.println( "ringed: " + smaH.toSmarts( (QueryAtomContainer) hs ) );

				String hsSMARTS3 = smaH.toSmarts( (QueryAtomContainer) hs);
				System.out.println( "hsSMARTS3 - " + hsSMARTS3 );
				
				
		
		Assert.assertEquals("O=C2NC1(CCCCC1)-C-C-C-2", hsSMARTS1);
		Assert.assertEquals("C[#6]3[#6](C1CC(C~C)CC2NC~CC12)c45:c:c(:c3[#6]:c:4)N-C6(O5)C(CCS~6)-O", hsSMARTS2);
		//fail("Not yet implemented");
	}
	
	
	@Test
	public void testRingNonRingAttributes() {
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		mcsOpts.bondFrequencies = true;
		
		mapper = setUpMapper( mcsOpts );
		mcsOpts.verbose = true;
		
		//String inputFileName2 = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/mos_mapping_test.smi";
		//String inputFileName1 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/bond_freq_test1.smi";
		String inputFileName1 = molDir + "ghost_example2.smi";
		String inputFileName2 = molDir + "mos_mapping_test.smi";
		String inputFileName3 = molDir + "mos_mapping_test6.smi";
		String inputFileName4 = molDir + "mddr/mddr_06233_MaxMin_10_2.smi";
		
		//IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName2, mcsOpts, 0 );


		int compoundLimit = 0;
		SmartsHelper gSmaH = new SmartsHelper(DefaultChemObjectBuilder.getInstance(), false);
		IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName4, mcsOpts, compoundLimit );
		
		mcsOpts.ringEnforcement = true;
		IQueryAtomContainer hs2re = createHyperstructure( mapper, inputFileName4, mcsOpts, compoundLimit );
		
		System.out.println( hs2.getAtomCount() + " " + hs2re.getAtomCount() + " " + hs2.getBondCount() + " " + hs2re.getBondCount() );
		System.out.println( gSmaH.toSmarts((QueryAtomContainer) hs2) );
		System.out.println( gSmaH.toSmarts((QueryAtomContainer) hs2re) );
		
		try {
			//getStandardMolecule(hs2);
			for( IAtom at : hs2.atoms() ) {
				at.setImplicitHydrogenCount(0);
			}
			
			for( IAtom at : hs2re.atoms() ) {
				at.setImplicitHydrogenCount(0);
			}
			
			 
			
			String hs2SMILES = sGenerator.create( hs2 ).replaceAll(":", "");
			String hs2reSMILES = sGenerator.create( hs2re ).replaceAll(":", "");
			
			IAtomContainer hs2s = sp.parseSmiles(hs2SMILES);
			IAtomContainer hs2res = sp.parseSmiles(hs2reSMILES);
			
			System.out.println( hs2SMILES );
			System.out.println( hs2reSMILES );
			
			System.out.println( uit.isIsomorph(hs2s, hs2res) );
		} catch (CDKException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		

	}
	
	
	
	
	
	
	@Test
	public void testRingEnforcementMapping() {
		
		ChemAxonMCS mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = true;
		mcsOpts.SMARTSHandling = false;
		
		mapper = new ChemAxonMCS( null, null, mcsOpts );
		
		// Simple example 
		String inputFileName = molDir + "mos_mapping_test4.smi";
		IQueryAtomContainer hs1 = createHyperstructure( mapper, inputFileName, mcsOpts, 2 );
		
		
		String inputFileName2 = molDir + "mos_mapping_test.smi";
		IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName2, mcsOpts, 0 );
		
		 // get molecules
        ArrayList<IAtomContainer> compounds2 = null;
		try {
			compounds2 = ConvenienceTools.getQueryMolecules( new File( inputFileName2 ), null );
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String inputFileName3 = molDir + "chembl751606.smi";
		IQueryAtomContainer hs3 = createHyperstructure( mapper, inputFileName3, mcsOpts, 0 );
		
		 // get molecules
        ArrayList<IAtomContainer> compounds3 = null;
		try {
			compounds3 = ConvenienceTools.getQueryMolecules( new File( inputFileName3 ), null );
			/*
			for( IAtomContainer cdkMol : compounds3 ) {
				try {
					AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cdkMol);
					cdkMol = new FixBondOrdersTool().kekuliseAromaticRings(cdkMol);
					cdkMol = AtomContainerManipulator.removeHydrogens( cdkMol ); // remove hydrogens - no need to match these
					getStandardMolecule(cdkMol);
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			*/
			ConvenienceTools.sortBySimilarityMCS(compounds3);
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		String inputFileName4 = molDir + "ergosterol_analogues.smi";
		IQueryAtomContainer hs4 = createHyperstructure( mapper, inputFileName4, mcsOpts, 0 );
		
		 // get molecules
        ArrayList<IAtomContainer> compounds4 = null;
		try {
			compounds4 = ConvenienceTools.getQueryMolecules( new File( inputFileName4 ), null );
			ConvenienceTools.sortBySimilarityMCS(compounds4);
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Now for the mappings
		
		mcsOpts.matchBonds = true;
		mcsOpts.SMARTSHandling = true;
		
		try {
			IQueryAtomContainer hsSM = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs1) );
			IAtomContainer cGraph = sp.parseSmiles("C(C)1NC(=O)CCC1");
			ConvenienceTools.countRings(cGraph);
			/*
			//calculateMapping( hs,  cGraph, mcsOpts );
			mapper.setMainMol(hsSM);
			mapper.setQueryMol(cGraph);
			mapper.execute();
			//System.out.println( mapper.mcsSize );
			System.out.println( "MCS = " + mapper.mcsSMARTS );
			Assert.assertEquals( 5, mapper.mcsSize);
			Assert.assertEquals( 1, mapper.fragmentSizes[0]);
			*/
			
			// should be 5 out of 8 bonds matching
			Assert.assertTrue( ! isSubstructure2(hsSM, cGraph) );
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		

		IQueryAtomContainer hsSM = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs2) );
		for( int cGraph = 0; cGraph < compounds2.size(); cGraph++ ) {
			//mcsOpts.SMARTSHandling = true;
				System.out.println( mcsOpts.SMARTSHandling );
				
			/*
				mapper.setMainMol(hsSM);
				mapper.setQueryMol( compounds2.get(cGraph) );
				mapper.execute();
				//System.out.println( mapper.mcsSize );
				
				Assert.assertEquals( compounds2.get(cGraph).getBondCount(), mapper.mcsSize);
				Assert.assertEquals( compounds2.get(cGraph).getBondCount(), mapper.fragmentSizes[0]);
			*/
				
				System.out.println( "subgraph test, compound " + cGraph + ": ");
				
				try {
					isSubstructure( sGenerator.create( compounds2.get(cGraph) ), smaH.toSmarts( (QueryAtomContainer) hsSM) );
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				ConvenienceTools.countRings( compounds2.get(cGraph) );
				Assert.assertTrue( isSubstructure2(hsSM, compounds2.get(cGraph)) );
		}
		
		/*
		hsSM = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs3) );
		//ConvenienceTools.correctAtomBondTypes( (IQueryAtomContainer) hsSM );
		//mcsOpts.verbose = true;
		for( int cGraph = 0; cGraph < compounds3.size(); cGraph++ ) {
			//mcsOpts.SMARTSHandling = true;
				//System.out.println( mcsOpts.SMARTSHandling );
			
				mapper.setMainMol(hsSM);
				mapper.setQueryMol( compounds3.get(cGraph) );
				mapper.execute();
				
				System.out.println( "MCS = " + mapper.mcsSMARTS );
				
				Assert.assertEquals( compounds3.get(cGraph).getBondCount(), mapper.mcsSize);
				Assert.assertEquals( compounds3.get(cGraph).getBondCount(), mapper.fragmentSizes[0]);
			
		}
		*/
		
		hsSM = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs4) );
		//ConvenienceTools.correctAtomBondTypes( (IQueryAtomContainer) hsSM );
		//mcsOpts.verbose = true;
		for( int cGraph = 0; cGraph < compounds4.size(); cGraph++ ) {
			//mcsOpts.SMARTSHandling = true;
				//System.out.println( mcsOpts.SMARTSHandling );
			
				/*
				mapper.setMainMol(hsSM);
				mapper.setQueryMol( compounds4.get(cGraph) );
				mapper.execute();
				
				System.out.println( "MCS = " + mapper.mcsSMARTS );
				
				Assert.assertEquals( compounds4.get(cGraph).getBondCount(), mapper.mcsSize);
				Assert.assertEquals( compounds4.get(cGraph).getBondCount(), mapper.fragmentSizes[0]);
			*/
				Assert.assertTrue( isSubstructure2(hsSM, compounds2.get(cGraph)) );
		}
	}
	
	
	
	
	@Test
	public void testNonRingEnforcementMapping() {
		
		System.out.println("Begin testNonRingEnforcementMapping");
		MCSMethods mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		//mcsOpts.matchAtoms = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		//mcsOpts.verbose = true;
		
		//mapper = new ChemAxonMCS( null, null, mcsOpts );
		mapper = setUpMapper(mcsOpts);
		
		// Simple example 
		String inputFileName = molDir + "/mos_mapping_test4.smi";
		IQueryAtomContainer hs1 = createHyperstructure( mapper, inputFileName, mcsOpts, 2 );
		
		String hsSMARTS1 = smaH.toSmarts( (QueryAtomContainer) hs1);
		System.out.println( "hsSMARTS - " + hsSMARTS1 );
		
		
		String inputFileName2 = molDir + "/mos_mapping_test.smi";
		IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName2, mcsOpts, 0 );
		hs2.setNotification(false);
		
		 // get molecules
        ArrayList<IAtomContainer> compounds2 = null;
		try {
			compounds2 = ConvenienceTools.getQueryMolecules( new File( inputFileName2 ), null );
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println( "compoudns2 " + compounds2.size() );
		//String inputFileName3 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mddr/mddr_78374_MaxMin_10.sdf";
		String inputFileName3 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mddr_maxmin_10_1.smi";
		//String inputFileName3 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mddr_mw300_5.smi";
		//String inputFileName3 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/aid466_actives_5.smi";
		//String inputFileName3 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/mos_mapping_test2.smi";
		//String inputFileName3 = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/ergosterol_analogues.smi";
		IQueryAtomContainer hs3 = createHyperstructure( mapper, inputFileName3, mcsOpts, 0 );
		
		 // get molecules
        ArrayList<IAtomContainer> compounds3 = null;
		try {
			compounds3 = ConvenienceTools.getQueryMolecules( new File( inputFileName3 ), null );
			//compounds3 = new ArrayList<IAtomContainer>( compounds3.subList(2, 3) );
			
			for( int c = 0; c < compounds3.size(); c++ ) {
				IAtomContainer cdkMol = compounds3.get(c);
				try {
					AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cdkMol);
					//cdkMol = new FixBondOrdersTool().kekuliseAromaticRings(cdkMol);
					//getStandardMolecule(cdkMol);
					ConvenienceTools.initializeMolecule(cdkMol);
					//cdkMol = AtomContainerManipulator.removeHydrogens( cdkMol );
					compounds3.set(c, cdkMol);
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			ConvenienceTools.sortBySimilarityMCS(compounds3);
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		
		// Now for the mappings
		
		mcsOpts.matchBonds = true;
		//mcsOpts.matchAtoms = false;
		mcsOpts.SMARTSHandling = true;
		
		try {
			IQueryAtomContainer hsSM = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs1) );
			//IAtomContainer cGraph = sp.parseSmiles("C(C)1NC(=O)CCC1");
			IAtomContainer cGraph = sp.parseSmiles("O=CNC1CCCCC1");
			/*
			//calculateMapping( hs,  cGraph, mcsOpts );
			mapper.setMainMol(hsSM);
			mapper.setQueryMol(cGraph);
			mapper.execute();
			//System.out.println( mapper.mcsSize );
			
			Assert.assertEquals( 8, mapper.mcsSize);
			Assert.assertEquals( 8, mapper.fragmentSizes[0]); 
			*/
			try {
				isSubstructure( sGenerator.create(cGraph), smaH.toSmarts( (QueryAtomContainer) hsSM) );
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println( "hs1 - " + smaH.toSmarts( (QueryAtomContainer) hs1) );
			Assert.assertTrue( isSubstructure2(hsSM, cGraph) );
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		

		IQueryAtomContainer hsSM = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs2) );
		for( int cGraph = 0; cGraph < compounds2.size(); cGraph++ ) {
			//mcsOpts.SMARTSHandling = true;
				//System.out.println( mcsOpts.SMARTSHandling );
			/*
				mapper.setMainMol(hsSM);
				mapper.setQueryMol( compounds2.get(cGraph) );
				mapper.execute();
				//System.out.println( mapper.mcsSize );
				
				Assert.assertEquals( compounds2.get(cGraph).getBondCount(), mapper.mcsSize);
				Assert.assertEquals( compounds2.get(cGraph).getBondCount(), mapper.fragmentSizes[0]);
				*/
				Assert.assertTrue( isSubstructure2(hsSM, compounds2.get(cGraph)) );
		}
		
		hsSM = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs3) );
		//hsSM = hs3;
		//ConvenienceTools.correctAtomBondTypes( (IQueryAtomContainer) hsSM );
		//mcsOpts.matchBonds = false;
		//mcsOpts.verbose = true;
		for( int cGraph = 0; cGraph < compounds3.size(); cGraph++ ) {
			//mcsOpts.SMARTSHandling = true;
				//System.out.println( mcsOpts.SMARTSHandling );
			
				mapper.setMainMol(hsSM);
				mapper.setQueryMol( compounds3.get(cGraph) );
				mapper.execute();
				
				try {
					isSubstructure( sGenerator.create(compounds3.get(cGraph)), smaH.toSmarts( (QueryAtomContainer) hs3) );
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println( "HS = " + smaH.toSmarts( (QueryAtomContainer) hsSM) );
				/*
				System.out.println( "MCS = " + mapper.mcsSMARTS );
				System.out.println( mapper.mcsSize + " " + mapper.fragmentSizes[0] + " out of " + compounds3.get(cGraph).getBondCount() + " bonds" );
				Assert.assertEquals( compounds3.get(cGraph).getBondCount(), mapper.mcsSize);
				Assert.assertEquals( compounds3.get(cGraph).getBondCount(), mapper.fragmentSizes[0]);
				*/
			
		}
		
		ConvenienceTools.printHyperstructureStats(hs3, compounds3.toArray( new IAtomContainer[0] ));
	}
	
	
	
	
	@Test
	public void testSingleAromaticCreation() {
		
		MCSMethods mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		
		mapper = setUpMapper( mcsOpts );
		
		// Simple example 
		String inputFileName = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/aromatic_nonaromatic_rings.smi";
		
		IQueryAtomContainer hs = createHyperstructure( mapper, inputFileName, mcsOpts, 0 );
		
		String hsSMARTS = smaH.toSmarts( (QueryAtomContainer) hs);
		System.out.println( "hsSMARTS benzene - " + hsSMARTS );
		
		
		Assert.assertEquals("[#6]1[#6][#6][#6][#6][#6]1", hsSMARTS );
		
	}
	
	
	@Test
	public void testAnyBondCreation() {
		
		MCSMethods mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		
		//mapper = setUpMapper( mcsOpts );
		mapper = setUpMapper( mcsOpts );
		
		// Simple example 
		String inputFileName = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/single_triple_bonds.smi";
		IQueryAtomContainer hs = createHyperstructure( mapper, inputFileName, mcsOpts, 0 );
		
		String hsSMARTS = smaH.toSmarts( (QueryAtomContainer) hs );
		System.out.println( "hsSMARTS benzene - " + hsSMARTS );
		
		
		Assert.assertEquals("C~C", hsSMARTS );
		
	}
	
	
	
	/*
	 * Seeing what happens when using different MCS type (i.e tdMCS)
	 */
	
	@Test
	public void testTopologyHS() {
		
		MCSMethods mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		
		//mapper = setUpMapper( mcsOpts );
		mapper = setUpMapper( mcsOpts );
		
		
		String inputFileName = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/ring_mapping_test1.smi";
		//String inputFileName = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/ergosterol_analogues.smi";
		IQueryAtomContainer hs = createHyperstructure( mapper, inputFileName, mcsOpts, 2 );
		
		String hsSMARTS = smaH.toSmarts( (QueryAtomContainer) hs );
		System.out.println( "hsSMARTS ringed - " + hsSMARTS );
		
		
		Assert.assertEquals("C~C", hsSMARTS );
		
	}
	
	
	/*
	 * Presumed methodology to test that the correct bond frequencies are assigned to hs edges:
	 * 
	 * - For each pair of molecules, generate MCES
	 * - copy/maintain previous bond weights
	 * - create hyperstructure with bond weights
	 * - use MCES to check that the bonds in the MCES have their weights incremented by one
	 * - thus the others NOT in the MCES should be the same
	 * 
	 */
	@Test
	public void testBondFrequency() {
		
		MCSMethods mapper = null;
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		mcsOpts.bondFrequencies = true;
		
		mapper = setUpMapper( mcsOpts );
		
		//IQueryAtomContainer deleteThis = smaP.parse("C1|1|C23~C|3|C~|3|C-|1|3(|2|C|1|1-|1|2)-|1|O");
		//IQueryAtomContainer deleteThis = smaP.parse("C|1|C3|3|C(|3|C1|3|C|3|C(|2|C~|2|C)|3|C|3|C2|3|N|3|C~|3|C|3|C|3|1|3|2)|2|c45:|3|c:|3|c(:|3|c|3|3|2|c:|1|c:|1|4)|2|N-|2|C6(|1|O|1|5)|2|C(|1|C|1|C|1|S~|3|6)-|1|O");
		//System.out.println( smaH.toSmarts( (QueryAtomContainer) deleteThis ) );
		
		//String inputFileName2 = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/mos_mapping_test.smi";
		String inputFileName2 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/bond_freq_test1.smi";
		//IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName2, mcsOpts, 0 );

		// get molecules
		ArrayList<IAtomContainer> compounds2 = null;
		try {
			compounds2 = ConvenienceTools.getQueryMolecules( new File( inputFileName2 ), null );
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

		
		String startingSMILES = null;
		try {
			startingSMILES = sGenerator.create( compounds2.get(0) );
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// before we start, account for the special case where an explicit hydrogen is added to aromatic nitrogens (don't want explicit hydrogens)
		startingSMILES = startingSMILES.replaceAll("\\[nH\\]", "n");
		
		IQueryAtomContainer hs = smaP.parse( startingSMILES  );
		ConvenienceTools.correctAtomBondTypes(hs);

		
		
		
		// Now for the mappings


	/*	for( int cGraph = 1; cGraph < compounds2.size(); cGraph++ ) {
			
			ArrayList<IAtomContainer> tempC = new ArrayList<IAtomContainer>();
			System.out.println("cGraph = " + cGraph);
			
			//tempC.add( hs );
			tempC.add( compounds2.get(0) );
			tempC.add( compounds2.get(cGraph) );
			
			ConvenienceTools.calculateHyperstructure(tempC, mapper, 33, 41, false, false, false);
			
			CDKSMARTSHyperstructureFitness hsf = new CDKSMARTSHyperstructureFitness(hs, compounds2.get(cGraph), false );
			hsf.getAllBondMatchHyperstructure(false) ;
			CDKSMARTSHyperstructureFitness.configureBondOrders(hs);
			
			mapper.setMainMol( hs );
			mapper.setQueryMol( compounds2.get(cGraph) );
			mapper.execute();
			
			ArrayList<ArrayList<Integer>> matches = mapper.getBestMatches();
			ArrayList<Integer> bm = matches.get(0) ;
			Map<IBond, IBond> bondMap = mapper.getBestBondMatches().get(0);
			hs = hsf.createHyperstructure( hsf.atomMapFromChromosome(bm), bondMap );
			hsf.restoreHyperstructureBonds(false);
			
			System.out.println( "hs : " + smaH.toSmarts( (QueryAtomContainer) hs) );


		}
*/
		
		Map<IBond, IBond> bondMapCopy;
		
		
		// first test - using the MCES, we should see that the bonds in the MCES should have their weights incremented by one per iteration
		for( int cGraph = 1; cGraph < compounds2.size(); cGraph++ ) {
			System.out.println("cGraph = " + cGraph);
			//IQueryAtomContainer hsCDKDupe = new QueryAtomContainer( hs, hs.getBuilder() );
			CDKSMARTSHyperstructureFitness hsf = new CDKSMARTSHyperstructureFitness(hs, compounds2.get(cGraph), false, true, cGraph );
			//hsf.getAllBondMatchHyperstructure(false) ;
			CDKSMARTSHyperstructureFitness.configureBondOrders(hs);
			
			mapper.setMainMol( hs );
			mapper.setQueryMol( compounds2.get(cGraph) );
			mapper.execute();
			
			List<List<Integer>> matches = mapper.getBestAtomIndexMatches();
			List<Integer> bm = matches.get(0) ;
			Map<IBond, IBond> bondMap = mapper.getBestBondMatches().get(0);
			
			// Register current hs bond frequency values here
			Map<IBond, Integer> prevFreqs = new HashMap<IBond, Integer>( bondMap.size() );
			for( IBond mB : hs.bonds() ) {
				prevFreqs.put(mB, (Integer) mB.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ) ); 
			}
			
			
			
			
			bondMapCopy = new LinkedHashMap<IBond, IBond>( bondMap );  // shallow copy for MCS reference
			hs = hsf.createHyperstructure( hsf.atomMapFromChromosome(bm), bondMap );
			//hsf.restoreHyperstructureBonds(false);
			
			System.out.println( "hs : " + smaH.toSmarts( (QueryAtomContainer) hs) );

			// see if weights in MCES are incremented properly
			for( IBond mB : hs.bonds() ) {
				// look for MCES bonds
				if( bondMapCopy.containsValue( mB ) ) {
					System.out.println( "in MCS, old vs new: " + prevFreqs.get(mB) + " - " + mB.getProperty("_bfreq") + "  " + mB );
					Assert.assertEquals( prevFreqs.get(mB) + 1, mB.getProperty("_bfreq") );
				} else {  // ones not in MCES should not have incremented
					System.out.println( "not in MCS, old vs new: " + prevFreqs.get(mB) + " - " + mB.getProperty("_bfreq") + "  " + mB );
					
					/* a flaw with this testing - as some bonds in the hyperstructure are changed, we lose track of them.
					 * They are ignored for the moment
					 */
					if( prevFreqs.containsKey(mB))
						Assert.assertEquals( prevFreqs.get(mB), mB.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ) );
				}
			}
		}
		
		
		
		// now to see if the hyperstructure created using the convenience method here, has the same bond frequencies
		IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName2, mcsOpts, 0 );
		for( int b = 0; b < hs.getBondCount(); b++ ) {
			IBond b1 = hs.getBond(b);
			IBond b2 = hs2.getBond(b);
			System.out.println( b1.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ) + " - " + b2.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType )  );
			System.out.println( b1.getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType ) + " - " + b2.getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType )  );
			Assert.assertEquals(b1.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ), b2.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ));
		}
		
		String hsSMARTS = smaH.toSmarts( (QueryAtomContainer) hs);
		System.out.println( hsSMARTS );
		
		
		
		
		
		
		// it'd be useful to check that the SMARTS is also correct
		// only thing I can think of doing ATM is to perform a "frequency comparison" on the bond frequencies - number of 1s, 2s, 3s etc
		// (because the mapping correspondence is broken, so we have to resort to frequencies)
		int[] hsFreqs = new int[ compounds2.size() + 1 ];
		int[] hsSMARTSFreqs = new int[ compounds2.size() + 1 ];
		
		// initialise these arrays
		for( int n = 0; n < compounds2.size(); n++ ) {
			hsFreqs[n] = 0;
			hsSMARTSFreqs[n] = 0;
		}
		
		
		// go through current hs frequencies
		for( IBond b : hs.bonds() ) {
			int freq = b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType );
			hsFreqs[ freq ]++;
		}
		
		// go through SMARTS hs frequencies
		 Pattern findFrequencies = Pattern.compile("\\|(\\d+)(.*?)\\|");
		 Matcher m = findFrequencies.matcher( hsSMARTS );
		 
		 while( m.find() ) {
			 //boolean b = m.matches();
			 int freq = Integer.parseInt(m.group(1));
			 hsSMARTSFreqs[ freq ]++;
			 //System.out.println( m.groupCount() + " " + " " + m.group(1) );
		 }
		 
		 // now to compare
		for( int n = 0; n < compounds2.size() + 1; n++ ) {
				System.out.println( "freq " + n + ": " + hsFreqs[n] + "," + hsSMARTSFreqs[n] );
				Assert.assertEquals(hsFreqs[n], hsSMARTSFreqs[n]);
		}
		
		
		
		// topology type test
		
		//mcsOpts.bondFrequencies = false;
		mcsOpts.ringEnforcement = true;
		String inputFileName3 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/chembl751606.smi";
		//String inputFileName3 = molDir + "mos_mapping_test.smi";
		//String inputFileName3 = molDir + "mos_topology_difference_test.smi";
		IQueryAtomContainer hsEnforced = createHyperstructure(mapper, inputFileName3, mcsOpts, 0);
		
		mcsOpts.ringEnforcement = false;
		IQueryAtomContainer hsNonenforced = createHyperstructure(mapper, inputFileName3, mcsOpts, 0);
		//ConvenienceTools.countRings(hs3);
		
		/*for( IBond hb : hs3.bonds() ) {
			if( hb instanceof SmartsBondExpression ) {
				SmartsBondExpression sb = (SmartsBondExpression) hb;
				System.err.println( "topo: " +  sb + " , " +  sb.getProperty(CDKSMARTSHyperstructureFitness.topologyType) + " , " + sb.getProperty(CDKSMARTSHyperstructureFitness.bondFrequencyType) );
			}
		}*/
		
		String enforcedSMARTS = smaH.toSmarts( (QueryAtomContainer) hsEnforced );
		String nonenforcedSMARTS = smaH.toSmarts( (QueryAtomContainer) hsNonenforced );
		System.out.println( "topology type SMARTS = " + enforcedSMARTS );
		System.out.println( "topology type SMARTS (no ring enforcement) = " + nonenforcedSMARTS );
		
		/*List<Integer> enforcedFreqs = new ArrayList<Integer>();
		List<Integer> nonenforcedFreqs = new ArrayList<Integer>();
		List<String> enforcedTopologies = new ArrayList<String>();
		List<String> nonenforcedTopologies = new ArrayList<String>();*/
		
		// check that SMARTS is also correct
		// frequency correspondence check (because the mapping correspondence is broken, so we have to resort to frequencies)
		int[] enforcedFreqs = new int[ 6 ];
		int[] nonenforcedFreqs = new int[ 6 ];		
		int[] enforcedTopologies = new int[ 3 ];
		int[] nonenforcedTopologies = new int[ 3 ];	
		
		//Pattern findFrequencies = Pattern.compile("\\|\\d+\\l?\\|");
		Matcher enforcedMatches = findFrequencies.matcher(enforcedSMARTS);
		
		while( enforcedMatches.find() ) {
					 int freq = Integer.parseInt(enforcedMatches.group(1));
					 String topType = enforcedMatches.group(2);
					 
					 //enforcedFreqs.add( freq );
					 //enforcedTopologies.add( enforcedMatches.group(2) );
					 
					 enforcedFreqs[freq]++;
					 
					 if( topType.equals("r") ) {
						 enforcedTopologies[0]++;
					 } else if( topType.equals("c") ) {
						 enforcedTopologies[1]++;
					 } else {
						 enforcedTopologies[2]++;
					 }
		 }
		
		Matcher nonenforcedMatches = findFrequencies.matcher(nonenforcedSMARTS);
		
		while( nonenforcedMatches.find() ) {
					 int freq = Integer.parseInt(nonenforcedMatches.group(1));
					 String topType = nonenforcedMatches.group(2);
					 
					 //nonenforcedFreqs.add( freq );
					 //nonenforcedTopologies.add( nonenforcedMatches.group(2) );
					 
					 nonenforcedFreqs[freq]++;
					 
					 if( topType.equals("r") ) {
						 nonenforcedTopologies[0]++;
					 } else if( topType.equals("c") ) {
						 nonenforcedTopologies[1]++;
					 } else {
						 nonenforcedTopologies[2]++;
					 }
		 }
		
		/*Assert.assertEquals( enforcedFreqs.size(), nonenforcedFreqs.size() );
		for( int f = 0; f < enforcedFreqs.size(); f++ ) {
			System.out.println( "enforced-nonenforced = " + enforcedFreqs.get(f) + enforcedTopologies.get(f) + "-" + nonenforcedFreqs.get(f) + nonenforcedTopologies.get(f) );
		}*/
		
		
		
		// compare ringEnforcement weights with non-ringEnforcement weights
		for( int f = 0; f < compounds2.size() + 1; f++ ) {
			System.out.println( "freq " + f + ": " + enforcedFreqs[f] + "," + nonenforcedFreqs[f] );
			Assert.assertEquals(enforcedFreqs[f], nonenforcedFreqs[f]);
		}
		
		
		// see if ringEnforcement topology types match non-ringEnforcement topology types
		for( int f = 0; f < 3; f++ ) {
			System.out.println( "topology " + f + ": " + enforcedTopologies[f] + "," + nonenforcedTopologies[f] );
			Assert.assertEquals(enforcedTopologies[f], nonenforcedTopologies[f]);
		}
		
		
		
		
		
		
		
		// special case mapping test to see if the frequencies are in the "right places"
		IAtomContainer testMol = null;
		try {
			testMol = sp.parseSmiles("C1=CC=C1");
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		mapper.setMainMol( hs2 );
		mapper.setQueryMol( testMol );
		mapper.execute();
		
		List<List<Integer>> matches = mapper.getBestAtomIndexMatches();
		List<Integer> bm = matches.get(0) ;
		Map<IBond, IBond> bondMap = mapper.getBestBondMatches().get(0);
		
		for( IBond b : bondMap.values() ) {
			System.out.println( "freq " + b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ) );
			//Assert.assertEquals(hsFreqs[n], hsSMARTSFreqs[n]);
		}
	}
	
	
	
	
	
	/**
	 * Presence of Ghost substructures - looking for a specific number for given hyperstructures
	 * 
	 * Also test for MCS ghost substructure thing.
	 * 
	 * TBH this is rather hard to test.  Other than manually inspecting the output I can't think of what to test for
	 */
	@Test
	public void testGhostFinding() {
		
		
		ChemAxonMCS.ChemAxonMCSOptions mcsOpts = new ChemAxonMCS.ChemAxonMCSOptions();
		mcsOpts.connectedMode = false;
		mcsOpts.matchBonds = false;
		mcsOpts.ringEnforcement = false;
		mcsOpts.SMARTSHandling = false;
		mcsOpts.bondFrequencies = true;
		
		mapper = setUpMapper( mcsOpts );
		mcsOpts.verbose = false;
		
		//String inputFileName2 = "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/mos_mapping_test.smi";
		//String inputFileName1 = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/bond_freq_test1.smi";
		String inputFileName1 = molDir + "ghost_example2.smi";
		String inputFileName2 = molDir + "mos_mapping_test.smi";
		String inputFileName3 = molDir + "mos_mapping_test6.smi";
		String inputFileName4 = molDir + "mddr/mddr_06233_MaxMin_10_2.smi";
		
		//IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName2, mcsOpts, 0 );


		SmartsHelper gSmaH = new SmartsHelper(DefaultChemObjectBuilder.getInstance(), false);
		IQueryAtomContainer hs2 = createHyperstructure( mapper, inputFileName2, mcsOpts, 0 );
		
		
		
		IQueryAtomContainer hs1 = createHyperstructure( mapper, inputFileName1, mcsOpts, 0 );
		try {
			//System.out.println( smaH.toSmarts(test) );
			HyperstructureConstructor.getStandardMolecule(hs1);
			System.out.println( gSmaH.toSmarts( (QueryAtomContainer) hs1));
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		
		for( int n=0; n < hs1.getAtomCount(); n++ ) {

			IAtomContainer[] radii = ConvenienceTools.getNeighbourhoodGraphEdgeInduced(hs1, n, 3);
			
			for( int r = 1; r < radii.length; r++ ) {
				IAtomContainer test = radii[r];
	
				System.out.print( n + " " + r + " - " + gSmaH.toSmarts( (QueryAtomContainer) test) );
					
				for( IBond b : test.bonds() ) {
					List<Integer> bondOrigins = (List<Integer>) b.getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType );
					System.out.print( " " + bondOrigins + "" + b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType )  );
					
					// This caused nasty problems before - the bond origins were causing issues
					Assert.assertEquals(bondOrigins.size(), b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ));
				}
				System.out.println( "  ghost - " + ConvenienceTools.isGhostSubstructure(test) );
			}
		}
		
		
		// SMARTS test
		String hs2SMARTS = smaH.toSmarts( (QueryAtomContainer) hs2 );
		IQueryAtomContainer hs2Redux = smaP.parse(hs2SMARTS);
		System.out.println( "hs2 smarts test - " + smaH.toSmarts( (QueryAtomContainer) hs2Redux ) );
		ConvenienceTools.correctAtomBondTypes(hs2Redux);
		
		// see if Ghost info is preserved from SMARTS output and re-parsing
		hs2 = hs2Redux;
		
		/*for( IBond b : hs2.bonds() ) {
			System.out.println( "hs2 " + b.getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType ) + "" + b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType )  );
		}*/
		
		HashMap<Integer, Boolean> uniqueSubstructuresHs2 = new HashMap<Integer, Boolean>();
		
		for( int n=0; n < hs2.getAtomCount(); n++ ) {

			IAtomContainer[] radii = ConvenienceTools.getNeighbourhoodGraphEdgeInduced(hs2, n, 3);
			
			for( int r = 1; r < radii.length; r++ ) {
				IAtomContainer test = radii[r];
	
				System.out.print( n + " " + r + " - " + gSmaH.toSmarts( (QueryAtomContainer) test) );
					
				for( IBond b : test.bonds() ) {
					List<Integer> bondOrigins = (List<Integer>) b.getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType );
					System.out.print( " " + bondOrigins + "" + b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ) + "  " + b.getProperty( ConvenienceTools.origBondIndicesProperty ) + ", "  );
					
					// This caused nasty problems before - the bond origins were causing issues
					Assert.assertEquals(bondOrigins.size(), b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ));
					
					boolean isGhost = ConvenienceTools.isGhostSubstructure(test);
					uniqueSubstructuresHs2.put( (Integer) b.getProperty( ConvenienceTools.origBondIndicesProperty ), isGhost ) ;
				}
				System.out.println( "  ghost - " + ConvenienceTools.isGhostSubstructure(test)  );
			}
		}
		
		
		
		int hsGhostCount = 0;
		for( boolean b : uniqueSubstructuresHs2.values() ) {
			if(b)
				hsGhostCount++;
		}
		//mcsGhostCount = uniqueGhosts.size();
		System.out.println( "Unique substructures in hs2 - " + uniqueSubstructuresHs2.size() + ", Unique Ghosts in MCS - " + hsGhostCount );
		
		
		
		
		
		// a test involving MCS ghost substructures
		
		class MCSGhostTest {
			
			public MCSGhostTest( String querySMILES, IAtomContainer hs, MCSMethods m, int r ) {
				try {
					queryMol = sp.parseSmiles( querySMILES );
					HyperstructureConstructor.getStandardMolecule(queryMol);
					
					System.out.println( "query smiles - " + sGenerator.create(queryMol) );
					System.out.println( "hs smarts test (again) - " + smaH.toSmarts( (QueryAtomContainer) hs ) );
				} catch (InvalidSmilesException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				catch (CDKException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				
				this.hs = hs;
				this.mapper = m;
				this.radius = r;
				
				//uniqueSubstructures = new HashMap<Integer, Boolean>();
			}
			
			public void test() {
				mcsGhostCount = 0;
				
				mapper.setMainMol( hs );
				mapper.setQueryMol( queryMol );
				
				mapper.execute();
				
				List<List<Integer>> matches = mapper.getBestAtomIndexMatches();
				//List<Integer> bm = matches.get(0) ;
				Map<IBond, IBond> bondMap = mapper.getBestBondMatches().get(0);
				
				//IAtomContainer mcs = ConvenienceTools.createCommonSubgraph( queryMol, hs, bondMap );
				IAtomContainer mcs = ConvenienceTools.createCommonSubgraph( hs, queryMol, bondMap );
			
					//System.out.println( "mcs - " + sGenerator.create(mcs) );
				System.out.println( "mcs - " + smaH.toSmarts( (QueryAtomContainer) mcs ) );
				
				
				
				//long mcsGhostTimePrior = System.currentTimeMillis();
				
				//HashSet<Integer> uniqueGhosts = new HashSet<Integer>();
				uniqueSubstructures = new HashMap<Integer, Boolean>();
				
				for( int n=0; n < mcs.getAtomCount(); n++ ) {

					IAtomContainer[] radii = ConvenienceTools.getNeighbourhoodGraphEdgeInduced(mcs, n, radius);
					
					for( int r = 1; r < radii.length; r++ ) {
						IAtomContainer test = radii[r];
			
						/*try {
							System.out.print( n + " " + r + " - " +  sGenerator.create(test) );
						} catch (CDKException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}*/
						System.out.print( n + " " + r + " - " + smaH.toSmarts( (QueryAtomContainer) test) );
							
						for( IBond b : test.bonds() ) {
							List<Integer> bondOrigins = (List<Integer>) b.getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType );
							System.out.print( " " + b.hashCode() + " " + bondOrigins + " " + b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType )  );
							
							// This caused nasty problems before - the bond origins were causing issues
							Assert.assertEquals(bondOrigins.size(), b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ));
						}
						
						boolean isGhost = ConvenienceTools.isGhostSubstructure(test);
						
						System.out.println( "  ghost - " + isGhost + "  " + test.getProperty( ConvenienceTools.origBondIndicesProperty ) );
						
						
						uniqueSubstructures.put( (Integer) test.getProperty( ConvenienceTools.origBondIndicesProperty ), isGhost );
					}
				}
				
				
				for( boolean b : uniqueSubstructures.values() ) {
					if(b)
						mcsGhostCount++;
				}
				//mcsGhostCount = uniqueGhosts.size();
			}
			
			
			IAtomContainer queryMol, hs;
			MCSMethods mapper;
			int radius = 4;
			
			public HashMap<Integer, Boolean> uniqueSubstructures;
			public int mcsGhostCount = 0;
		};
		
		
		
		long mcsGhostTimePrior = System.currentTimeMillis();
		MCSGhostTest mcsGhostTest2 = new MCSGhostTest("C1CSC(C1)NC1CCCC2NC=CC12", hs2, mapper, 4) ;
		mcsGhostTest2.test();
		
		System.out.println( "Unique substructures in MCS - " + mcsGhostTest2.uniqueSubstructures.size() + ", Unique Ghosts in MCS - " + mcsGhostTest2.mcsGhostCount + " obtained at milliseconds - " + (System.currentTimeMillis() - mcsGhostTimePrior) );
		Assert.assertEquals(4, mcsGhostTest2.mcsGhostCount);
		
		
		
		
		
		
		// Now for a topology enforcement example
		mcsOpts.ringEnforcement = true;
		
		IQueryAtomContainer hs3 = createHyperstructure( mapper, inputFileName3, mcsOpts, 0 );
	
		MCSGhostTest mcsGhostTest3 = new MCSGhostTest("C=CNC(=C)CN(=O)=O", hs3, mapper, 4) ;
		mcsGhostTest3.test();
		
		System.out.println( "Unique substructures in MCS3 - " + mcsGhostTest3.uniqueSubstructures.size() + ", Unique Ghosts in MCS - " + mcsGhostTest3.mcsGhostCount   );

		
		
		
		
		// More complex
		mcsOpts.ringEnforcement = false;
		
		
				
		IQueryAtomContainer hs4 = createHyperstructure( mapper, inputFileName4, mcsOpts, 0 );
		hs4 = smaP.parse( smaH.toSmarts( (QueryAtomContainer) hs4 ) );
		ConvenienceTools.correctAtomBondTypes(hs4);
		
		mcsOpts.SMARTSHandling = true;  // MUST be on otherwise the degenerate bonds etc aren't used in MCS mapping!
		mcsOpts.matchBonds = true;
		// CC(=CCCC(=CCCC(=CCCC(=CCc1cc(O)ccc1O)C)C)C)C
		// CN(C)CCCN1c2ccccc2CCc3ccc(Cl)cc13
		MCSGhostTest mcsGhostTest4 = new MCSGhostTest("CN(C)CCCN1c2ccccc2CCc3ccc(Cl)cc13", hs4, mapper, 4) ;
		mcsGhostTest4.test();
		
		//System.out.println( "HS SMARTS - " + gSmaH.toSmarts( (QueryAtomContainer) hs4) );
		System.out.println( "Unique substructures in MCS4 - " + mcsGhostTest4.uniqueSubstructures.size() + ", Unique Ghosts in MCS - " + mcsGhostTest4.mcsGhostCount + " " + mapper.getBestBondMatches().get(0).size()  );
				
	}
	
	
	
	
	private MCSMethods setUpMapper( ChemAxonMCS.ChemAxonMCSOptions mcsOpts ) {
		
		boolean method1 = false;
		MCSMethods mapper = null;
		
		if( method1 ) {
			
			mapper = new ChemAxonMCS( null, null, mcsOpts );
			
		} else {
			CliqueDetection.ModularProductOptions mpOpts = new CliqueDetection.ModularProductOptions(false, false, 1);
			mapper = new DepolliCliqueDetection(mpOpts);
			//mapper = new FMCS();
			mapper.setMatchBonds( mcsOpts.matchBonds );
			
		}
		
		
		return mapper;
	}
	
	
	private IQueryAtomContainer createHyperstructure( MCSMethods mapper, String inputFileName, ChemAxonMCS.ChemAxonMCSOptions mcsOpts, int maxMols ) {
		// get molecules
        ArrayList<IAtomContainer> compounds = null;
		try {
			compounds = ConvenienceTools.getQueryMolecules( new File( inputFileName ), null );
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		/*
		for( int m=0; m < compounds.size(); m++ ) {
			IAtomContainer  mol = compounds.get(m);
			try {
				compounds.set( m, sp.parseSmiles( sGenerator.create(mol) ) );
			} catch (InvalidSmilesException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		*/
        // sort molecules
        //ConvenienceTools.sortBySimilarityMCS(compounds);
        //Collections.sort(compounds, ConvenienceTools.BONDASCENDING);
		
		
		ConvenienceTools.sortBySimilarityMCS( compounds );
		//Collections.reverse(compounds);
		
		if( maxMols > 0 ) {
			compounds = new ArrayList<IAtomContainer>( compounds.subList(0, maxMols) );
		}
		
		System.out.println( "there are " + compounds.size() + " compounds" );
		
		//IQueryAtomContainer hs = (IQueryAtomContainer) ConvenienceTools.calculateHyperstructure(compounds, mapper, 30, 32, false, mcsOpts.matchBonds, mcsOpts.ringEnforcement, mcsOpts.bondFrequencies );
		IQueryAtomContainer hs = (IQueryAtomContainer) HyperstructureConstructor.createHyperstructure(compounds, mapper, mcsOpts.ringEnforcement, mcsOpts.bondFrequencies );
		

		
		return hs;
	}
	
	
	 private boolean isSubstructure(String querySMILES, String targetSMARTS) {
		 Molecule query = null, target = null;
		 
         try {
             MolSearch s = new MolSearch();
 
             
            try {
 				MolImporter mim1 = new MolImporter( new StringBufferInputStream( querySMILES ), "smiles" );
 				MolImporter mim2 = new MolImporter( new StringBufferInputStream( targetSMARTS ), "smarts" );
 				query = mim1.read();
 				target = mim2.read();
 			} catch (IOException e2) {
 				// TODO Auto-generated catch block
 				e2.printStackTrace();
 			}
             
             s.setQuery(query);
             s.setTarget(target);
 
             // search all matching substructures and print hits
             int[] hitAtoms = null;
 
             //s.findAllHits();
             s.findAll();
             hitAtoms = s.findFirst();
             if (hitAtoms == null)
                 System.out.println("No hits");
             else {
                 System.out.println( "There are " + hitAtoms.length + " matched atoms, out of " + query.getAtomCount() + " total atoms in query" );
             }// end else
         } catch (SearchException e) {
             e.printStackTrace();
             System.exit(1);
         }// end catch
         
         return false;
     }// end method
	 
	 
	 private boolean isSubstructure2( IQueryAtomContainer hs, IAtomContainer query ) {
		 
		 ChemAxonMCSOptions caOpts = new ChemAxonMCSOptions();
		 caOpts.ringEnforcement = true;
		 ChemAxonMCS mapper = new ChemAxonMCS( hs, query, caOpts );
		 
		 CliqueDetection.ModularProductOptions mpOpts = new CliqueDetection.ModularProductOptions(false, false, 20);
		 //MCSMethods mapper = new DepolliCliqueDetection(mpOpts);
		 //MCSMethods mapper = new FMCS();
		 mapper.setMatchBonds( true );
		 mapper.setMainMol(hs);
		 mapper.setQueryMol(query);
		
		 
		 ConvenienceTools.correctAtomBondTypes(hs);
		 ConvenienceTools.correctAtomBondTypes(query);
			
		 mapper.execute();
		 
		 System.out.println( "HS = " + smaH.toSmarts( (QueryAtomContainer) hs ) );
		 try {
			System.out.println( "query = " + sGenerator.create( query ) );
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		 System.out.println( "MCS = " + mapper.mcsSMARTS ); 
		 System.out.println( mapper.mcsSize + " bonds identified out of " + query.getBondCount() );
         
         return ( mapper.mcsSize == query.getBondCount() && mapper.fragmentSizes[0] == query.getAtomCount() );
     }// end method


	
	
	
	//private String molDir = "/home/edmund/Documents/work/PhD/workspace/cisrg/data/input/";
	private String molDir = "/home/edmund/Documents/workspace/cisrg/data/input/";
	
	MCSMethods mapper;
	SmilesGenerator sGenerator;
	SmilesParser sp;
	SmartsHelper smaH;
	SmartsParser smaP;
	UniversalIsomorphismTester uit;
}
