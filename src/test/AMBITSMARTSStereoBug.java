package test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import org.cisrg.ambit.SmartsHelper;
import org.cisrg.ambit.SmartsParser;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class AMBITSMARTSStereoBug {

	public static void main( String[] argv ) {
		
		File inputFile = new File( "/usr/users/people/edmund/Documents/PhD/workspace/cisrg/data/input/ergosterol_analogues.smi" );
		
		ArrayList<IAtomContainer> queries = new ArrayList<IAtomContainer>(5);
		
		IteratingSMILESReader isr = null;
		try {
			isr = new IteratingSMILESReader( new FileInputStream(inputFile), DefaultChemObjectBuilder.getInstance() );
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		while( isr.hasNext() ) {
			try {
				IAtomContainer molecule = isr.next();

			    try {
			        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
			        //CDKHueckelAromaticityDetector.detectAromaticity(molecule);
			       // aromaticity.apply(molecule);
			    } catch (CDKException e) {
			        //logger.debug(e.toString());
			        throw new CDKException(e.toString(), e);
			    }
				
				queries.add( molecule );

			} catch( Exception e2 ) {
				e2.printStackTrace();
			}
		}
		
		try {
			isr.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		SmartsParser smartsParser = new SmartsParser( );
		SmilesGenerator sGenerator = new SmilesGenerator().isomeric().aromatic();
		
		try {
			for( IAtomContainer mol : queries ) {
				String refSMILES = sGenerator.create(mol);
				
				System.out.println( refSMILES ); 
				
				IQueryAtomContainer ref = smartsParser.parse( refSMILES );
			}
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
}
