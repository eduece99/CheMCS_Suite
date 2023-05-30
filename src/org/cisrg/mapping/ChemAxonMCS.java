package org.cisrg.mapping;



//import org.knime.cisrg.hyperstructures.GAPlugins;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;






public class ChemAxonMCS extends MCSMethods {

	
	public static class ChemAxonMCSOptions {
		
		public boolean matchBonds = true;
		public boolean matchAtoms = true;
		public boolean connectedMode = false;
		public boolean ringEnforcement = false;
		public boolean SMARTSHandling = true;
		public boolean bondFrequencies = false;
		public boolean verbose = false;
		
		
	}
	
	
	public ChemAxonMCS(IAtomContainer h, IAtomContainer q, ChemAxonMCSOptions opts) {
		super(  );
	}
	
	 
	@Override
	public void search( IAtomContainer graph1, IAtomContainer graph2 ) throws CDKException {
		 
		throw new CDKException("Error - this algorithm has not been released here for IP reasons" );
		
	}
	




	
}
