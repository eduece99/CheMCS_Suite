package org.cisrg.mapping;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.xmlbeans.impl.piccolo.io.FileFormatException;
import org.cisrg.BitSetExtended;
import org.cisrg.ambit.AromaticSymbolQueryAtom;
import org.cisrg.ambit.DoubleNonAromaticBond;
import org.cisrg.ambit.RingQueryBond;
import org.cisrg.ambit.SmartsAtomExpression;
import org.cisrg.ambit.SmartsBondExpression;
import org.cisrg.ambit.SmartsHelper;
import org.cisrg.ambit.SmartsParser;
import org.cisrg.ambit.TripleBondAromaticityNotSpecified;
import org.cisrg.hyperstructures.CDKSMARTSHyperstructureFitness;
import org.cisrg.mapping.ChemAxonMCS.ChemAxonMCSOptions;
import org.cisrg.mapping.VentoFoggia.VentoFoggia2;
import org.cisrg.utilities.HungarianAlgorithm;
//import org.knime.cisrg.hyperstructures.CDKHybridHyperstructureFitness;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IElectronContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.iterator.DefaultIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.isomorphism.matchers.CTFileQueryBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.SymbolQueryAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AnyAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AnyOrderQueryBond;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticOrSingleQueryBond;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticQueryBond;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticSymbolAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AtomicNumberAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.LogicalOperatorAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.RingBond;
import org.openscience.cdk.isomorphism.matchers.smarts.SMARTSAtom;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smsd.algorithm.matchers.DefaultBondMatcher;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;
import org.openscience.smsd.algorithm.matchers.DefaultAtomMatcher;

import ambit2.smarts.SmartsExpressionToken;

public class ConvenienceTools {


	/**
	 * Sort the input set of molecules by similarity, starting with the largest molecule.  The second molecule would thus be the one which is most
	 * similar to the first, then the 3rd would be most similar to the 2nd etc.
	 * 
	 * @param molecules
	 */
	public static void sortBySimilarity( ArrayList<IAtomContainer> molecules ) {
		
		ExtendedFingerprinter ef = new ExtendedFingerprinter();
		final String fpId = "_fp";
		final String rankId = "_rank";
		
		IAtomContainer biggest = Collections.max( molecules, ConvenienceTools.ATOMASCENDING );  // This, returns the maximum size 
		//System.out.println("biggest - " + biggest.getAtomCount() );
		biggest.setProperty(rankId, 1);
		int count = 1;
		
		for( IAtomContainer mol : molecules ) {
			try {
				mol.setProperty( fpId, ef.getBitFingerprint(mol) );
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		IAtomContainer current = biggest;
		int currentId = -1;
		while( count < molecules.size() ) {
			
			double maxSim = 0;
			for( int m = 0; m < molecules.size(); m++ ) {
				if( molecules.get(m).getProperty(rankId) == null ) {
					double sim = Tanimoto.calculate( 
							(IBitFingerprint) current.getProperty(fpId), 
							(IBitFingerprint) molecules.get(m).getProperty(fpId) 
					);
				
					if( sim > maxSim ) {
						maxSim = sim;
						currentId = m;
						current = molecules.get(currentId);
					}
				}
			}
			
			molecules.get(currentId).setProperty(rankId, ++count);
			
		}
		
		// descending order of rank
		Collections.sort( molecules, new Comparator<IAtomContainer>() {
	
			public int compare(IAtomContainer o1, IAtomContainer o2) {
	
				if( (Integer) o2.getProperty(rankId) > (Integer) o1.getProperty(rankId) )
					return -1;
	
				if( (Integer) o2.getProperty(rankId) < (Integer) o1.getProperty(rankId) )
					return 1;
	
				return 0;
			}
			
		});
		
		/*
		for( int m = 0; m < molecules.size(); m++ ) {
			double sim = Tanimoto.calculate( 
					(IBitFingerprint) molecules.get(0).getProperty(fpId), 
					(IBitFingerprint) molecules.get(m).getProperty(fpId) 
			);
			
			//logger.debug( m + " " + sim );
		}
		*/
	}
	
	
	
public static void sortBySimilarityMCS( ArrayList<IAtomContainer> molecules ) {
		
		SmilesGenerator sGenerator = new SmilesGenerator().aromatic();
		SmartsParser smartsParser = new SmartsParser();
		
		ChemAxonMCSOptions mcsOpts = new ChemAxonMCSOptions();
		mcsOpts.matchBonds = false;  // set to false to make topology-enforcement and non-enforced ones comparable
		mcsOpts.ringEnforcement = false;
		mcsOpts.connectedMode = false;
		//mcsOpts.verbose = true;
		MCSMethods mcs = new ChemAxonMCS(null, null, mcsOpts);
		//MCSMethods mcs = new KawabataBuildupMCS(-1, false, false);
		//mcs.setMatchBonds(false);
		
		//mcs.setSmartsHandling(true);  // needed for aromatic representation
		final String rankId = "_rank";
		
		if( molecules.size() <= 2 ) {
			Collections.sort( molecules, ConvenienceTools.ATOMASCENDING );
			return;
		}
		
		IAtomContainer biggest = Collections.max( molecules, ConvenienceTools.ATOMASCENDING );  // This, returns the maximum size 
		//System.out.println("biggest - " + biggest.getAtomCount() );
		biggest.setProperty(rankId, 1);
		int count = 1;
		
		
		IAtomContainer current = biggest;
		int currentId = -1;
		while( count < molecules.size() ) {
			
			double maxSim = 0;
			
			for( int m = 0; m < molecules.size(); m++ ) {
				if( molecules.get(m).getProperty(rankId) == null ) {
					
					/*String refSMILES = null;
					try {
						refSMILES = sGenerator.create( current );
					} catch (CDKException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					System.out.println( m + " | " + refSMILES );
					IQueryAtomContainer ref = smartsParser.parse( refSMILES );*/
					IAtomContainer ref = current;
					
					correctAtomBondTypes(ref);
					mcs.setMainMol(ref);
					mcs.setQueryMol( molecules.get(m) );
					mcs.execute();
					
					//double sim = (double) mcs.mcsSize / ( (double) ref.getBondCount() + (double) molecules.get(m).getBondCount() - (double) mcs.mcsSize );
					
					
					double sim = mcs.mcsSize;
					
					//System.out.println( count + " | " +  refSMILES + " " + sGenerator.createSMILES(molecules.get(m)) + " " + ref.getBondCount() + " " + molecules.get(m).getBondCount() + " " + mcs.mcsSize + " " + sim );
					if( sim > maxSim ) {
						maxSim = sim;
						currentId = m;
					}
				}
			}
			
			current = molecules.get(currentId);
			molecules.get(currentId).setProperty(rankId, ++count);
			
		}
		
		// ascending order of rank
		Collections.sort( molecules, new Comparator<IAtomContainer>() {
	
			public int compare(IAtomContainer o1, IAtomContainer o2) {
	
				if( (Integer) o2.getProperty(rankId) > (Integer) o1.getProperty(rankId) )
					return -1;
	
				if( (Integer) o2.getProperty(rankId) < (Integer) o1.getProperty(rankId) )
					return 1;
	
				return 0;
			}
			
		});
		
		for( int m = 0; m < molecules.size(); m++ ) {
			String mSMILES = null;
			try {
				mSMILES = sGenerator.create(molecules.get(m));
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println( m + " | " + molecules.get(m).getProperty(rankId) + " " + mSMILES  );
			//mSMILES = mSMILES.replaceAll("\\[nH\\]", "n");
			//System.out.println( m + " | " + molecules.get(m).getProperty(rankId) + " " + mSMILES  );
		}
		
		/*
		for( int m = 0; m < molecules.size(); m++ ) {
			double sim = Tanimoto.calculate( 
					(IBitFingerprint) molecules.get(0).getProperty(fpId), 
					(IBitFingerprint) molecules.get(m).getProperty(fpId) 
			);
			
			//logger.debug( m + " " + sim );
		}
		*/
	}
	

	public static Comparator<IAtomContainer> BONDASCENDING = new Comparator<IAtomContainer>() {
	
		@Override
		public int compare(IAtomContainer o1, IAtomContainer o2) {
	
			if( o2.getBondCount() > o1.getBondCount() )
				return -1;
	
			if( o2.getBondCount() < o1.getBondCount() )
				return 1;
	
			return 0;
		}
	};
	public static Comparator<IAtomContainer> ATOMASCENDING = new Comparator<IAtomContainer>() {
	
		@Override
		public int compare(IAtomContainer o1, IAtomContainer o2) {
	
			if( o2.getAtomCount() > o1.getAtomCount() )
				return -1;
	
			if( o2.getAtomCount() < o1.getAtomCount() )
				return 1;
	
			return 0;
		}
	};
	
	
	public static ArrayList<IAtomContainer> getQueryMolecules( File inputFile, Comparator<IAtomContainer> comparator, boolean processCompounds ) throws FileFormatException {
		
		
		
		DefaultIteratingChemObjectReader<IAtomContainer> isr = null;
		ArrayList<IAtomContainer> queries = new ArrayList<IAtomContainer>(20);
	
		/*
		// Atom symbol to atomic number hash hack
		HashMap<String, Integer> symbolToNumber = new HashMap<String, Integer>();
		symbolToNumber.put("H", 1);
		symbolToNumber.put("C", 6);
		symbolToNumber.put("N", 7);
		symbolToNumber.put("O", 8);
		symbolToNumber.put("P", 15);
		symbolToNumber.put("S", 16);
		symbolToNumber.put("F", 9);
		symbolToNumber.put("Cl", 17);
		symbolToNumber.put("Br", 35);
		symbolToNumber.put("I", 53);
		*/
		
		// open connection to input SMILES file
		try {
			
			if( ! inputFile.exists() || ! inputFile.isFile() || ! inputFile.canRead() ) {
				throw new FileNotFoundException("Error - Molecule file does not exist");
			}
			
			if( inputFile.getName().matches(".*sdf$") ) {
				isr = new IteratingSDFReader( new FileInputStream(inputFile), DefaultChemObjectBuilder.getInstance() );
			} else { 
				isr = new IteratingSMILESReader( new FileInputStream(inputFile), DefaultChemObjectBuilder.getInstance() );
			}
			
		} catch (FileNotFoundException e3) {
			// TODO Auto-generated catch block
			e3.printStackTrace();
		}
		
		if( ! isr.hasNext() ) {
			//System.err.println( "Error - molecule file contains no readable molecules" );
			throw new FileFormatException("Error - molecule file contains no readable molecules " + isr);
		}
		
		// load the query molecules
		while( isr.hasNext() ) {
			try {
				//queries[0] = new AtomContainer( sp.parseSmiles("ONC1C(CO)CCC1") );  // hyperstructure
				//queries[1] = new AtomContainer( sp.parseSmiles("ONC1CCC(CCC)C1") );  // query
				IAtomContainer molecule = isr.next();
				
				// Sometimes new lines are read as molecules so these must be ignored
				if( molecule.getAtomCount() <= 0 )
					continue;
				
				
				if( processCompounds ) {
					 // check for aromaticity
				    // a good model for writing SMILES
				    Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(),
				                                              Cycles.all());
				    
				    try {
				        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
				        //CDKHueckelAromaticityDetector.detectAromaticity(molecule);
				        calculateImplicitHydrogens(molecule);
				        correctAtomBondTypes(molecule);
				        countRings(molecule);
				        aromaticity.apply(molecule);
				    } catch (CDKException e) {
				        //logger.debug(e.toString());
				        throw new CDKException(e.toString(), e);
				    }
				    
				}
					
					queries.add( molecule );
					//System.err.println("lolcat");
					//AtomTypeTools att = new AtomTypeTools();
					//att.assignAtomTypePropertiesToAtom(rdmol1);
					//att.assignAtomTypePropertiesToAtom(rdmol2);
					
				
			
	
	
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
		
		if( queries == null || queries.size() == 0 )
			throw new FileFormatException("Error - no molecules read in file");
	
		// sort in ascending order of number of bonds - Nathan Brown found this to give best compressions
		if( comparator != null )
			Collections.sort(queries, comparator);
	
		//System.out.println( queries.get(0).getBondCount() + " " + queries.get(2).getBondCount() );
	
	
		
	
		return queries;
	}
	
	public static ArrayList<IAtomContainer> getQueryMolecules( File inputFile, Comparator<IAtomContainer> comparator ) throws FileFormatException {
		return getQueryMolecules( inputFile, comparator, true );
	}
	
	
	/**
	 * Returns a set of connected components identified in the AtomContainer supplied
	 * 
	 * E.g - break a chain between 2 rings in a molecule.  This will yield 2 connected components, though CDK won't know this so you need
	 * to use this method to identify them.
	 * 
	 * FIXME  I've made this method synchronized.  Too high a risk of having the same AtomContainer object being used by this method in parallel.
	 * 
	 * @param atomContainer
	 * @return
	 */
	public static IAtomContainerSet partitionIntoMolecules(IAtomContainer atomContainer) {
		
		
		IAtomContainer newContainer = null;
		IAtomContainerSet molecules = null;
		
		if( atomContainer instanceof IQueryAtomContainer ) {
				//newContainer = new QueryAtomContainer2( atomContainer, atomContainer.getBuilder() ); 
				newContainer = atomContainer.getBuilder().newInstance(IAtomContainer.class, atomContainer);
				molecules = new AtomContainerSet();
			
		} else {
			newContainer = atomContainer.getBuilder().newInstance(IAtomContainer.class);
			molecules = atomContainer.getBuilder().newInstance(IAtomContainerSet.class);
		}
		
		newContainer.setNotification(false);  // avoid listeners from causing issues
		
		IAtomContainer molecule;
		
		List<IAtom> sphere = new ArrayList<IAtom>();

        for (IAtom atom : atomContainer.atoms()) {
            atom.setFlag(CDKConstants.VISITED, false);
            newContainer.addAtom(atom);
            newContainer.removeListener(newContainer);
            newContainer.removeListener(atomContainer);
        }

        for (IBond bond : atomContainer.bonds()) {
            bond.setFlag(CDKConstants.VISITED, false);
            newContainer.addBond(bond);
        }

        for (IElectronContainer eContainer : atomContainer.electronContainers()) {
            eContainer.setFlag(CDKConstants.VISITED, false);
            newContainer.addElectronContainer(eContainer);
        }

        while(newContainer.getAtomCount() > 0) {
			IAtom atom = newContainer.getAtom(0);
			molecule = null;
			/*
			if( atomContainer instanceof IQueryAtomContainer ) {
				molecule = new QueryAtomContainer( atomContainer );
			} else {
				molecule = atomContainer.getBuilder().newInstance(IAtomContainer.class);
			}
			*/
			//TempACFix nC = new TempACFix( atomContainer );
			molecule = atomContainer.getBuilder().newInstance(IAtomContainer.class);
			
			sphere.clear();
			sphere.add(atom);
			atom.setFlag(CDKConstants.VISITED, true);
			PathTools.breadthFirstSearch(newContainer, sphere, molecule);
			molecules.addAtomContainer(molecule);
			newContainer.remove(molecule);
		}
        
      /*  for( IAtom at : atomContainer.atoms() ) {
        	for( IAtomContainer ac : molecules.atomContainers() )
        		at.removeListener(ac);
        }*/
        
        
        //System.out.println( molecules.getAtomContainerCount() );
		return molecules;
	}
	
	
	public static void initializeMolecule( IAtomContainer atomContainer ) throws CDKException {
	    // Code copied from
	    // org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
		
		
	    Map<String, Integer> valencesTable = new HashMap<String, Integer>();
	    valencesTable.put("H", 1);
	    valencesTable.put("Li", 1);
	    valencesTable.put("Be", 2);
	    valencesTable.put("B", 3);
	    valencesTable.put("C", 4);
	    valencesTable.put("N", 5);
	    valencesTable.put("O", 6);
	    valencesTable.put("F", 7);
	    valencesTable.put("Na", 1);
	    valencesTable.put("Mg", 2);
	    valencesTable.put("Al", 3);
	    valencesTable.put("Si", 4);
	    valencesTable.put("P", 5);
	    valencesTable.put("S", 6);
	    valencesTable.put("Cl", 7);
	    valencesTable.put("K", 1);
	    valencesTable.put("Ca", 2);
	    valencesTable.put("Ga", 3);
	    valencesTable.put("Ge", 4);
	    valencesTable.put("As", 5);
	    valencesTable.put("Se", 6);
	    valencesTable.put("Br", 7);
	    valencesTable.put("Rb", 1);
	    valencesTable.put("Sr", 2);
	    valencesTable.put("In", 3);
	    valencesTable.put("Sn", 4);
	    valencesTable.put("Sb", 5);
	    valencesTable.put("Te", 6);
	    valencesTable.put("I", 7);
	    valencesTable.put("Cs", 1);
	    valencesTable.put("Ba", 2);
	    valencesTable.put("Tl", 3);
	    valencesTable.put("Pb", 4);
	    valencesTable.put("Bi", 5);
	    valencesTable.put("Po", 6);
	    valencesTable.put("At", 7);
	    valencesTable.put("Fr", 1);
	    valencesTable.put("Ra", 2);
	    valencesTable.put("Cu", 2);
	    valencesTable.put("Mn", 2);
	    valencesTable.put("Co", 2);
	    
		// correct atomic numbers - CDK bug
	    for (int c = 0; c < atomContainer.getAtomCount(); c++) {
				
				atomContainer.getAtom(c).setAtomicNumber( 
						PeriodicTable.getAtomicNumber( atomContainer.getAtom(c).getSymbol() ) 
				);

		}
	    
	    correctAtomBondTypes( atomContainer );

			// set atom ids in the query for the fitness function to use
		for (IAtom at : atomContainer.atoms()) {
				at.setID( atomContainer.getAtomNumber(at) + "");
		}
	
	    // do all ring perception
	    AllRingsFinder arf = new AllRingsFinder();
	    IRingSet allRings;
	    try {
	        allRings = arf.findAllRings(atomContainer);
	    } catch (CDKException e) {
	        //logger.debug(e.toString());
	        throw new CDKException(e.toString(), e);
	    }
	
	    // sets SSSR information
	    SSSRFinder finder = new SSSRFinder(atomContainer);
	    IRingSet sssr = finder.findEssentialRings();
	
	    for (IAtom atom : atomContainer.atoms()) {
	
	        // add a property to each ring atom that will be an array of
	        // Integers, indicating what size ring the given atom belongs to
	        // Add SSSR ring counts
	        if (allRings.contains(atom)) { // it's in a ring
	            atom.setFlag(CDKConstants.ISINRING, true);
	            
	            // lets find which ring sets it is a part of
	            List<Integer> ringsizes = new ArrayList<Integer>();
	            IRingSet currentRings = allRings.getRings(atom);
	            int min = 0;
	            for (int i = 0; i < currentRings.getAtomContainerCount(); i++) {
	                int size = currentRings.getAtomContainer(i).getAtomCount();
	                if (min > size) min = size;
	                ringsizes.add(size);
	            }
	            atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
	            atom.setProperty(CDKConstants.SMALLEST_RINGS, sssr.getRings(atom));
	        } else {
	            atom.setFlag(CDKConstants.ISINRING, false);
	        }
	
	        // determine how many rings bonds each atom is a part of
	        int hCount;
	        if (atom.getImplicitHydrogenCount() == CDKConstants.UNSET) hCount = 0;
	        else hCount = atom.getImplicitHydrogenCount();
	        
	        if( hCount < 0 )
	        	hCount = 0;
	
	        List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);
	        int total = hCount + connectedAtoms.size();
	        for (IAtom connectedAtom : connectedAtoms) {
	            if (connectedAtom.getSymbol().equals("H")) {
	                hCount++;
	            }
	        }
	        atom.setImplicitHydrogenCount(hCount);
	        atom.setProperty(CDKConstants.TOTAL_CONNECTIONS, total);
	        atom.setProperty(CDKConstants.TOTAL_H_COUNT, hCount);
	        /*
	        if (valencesTable.get(atom.getSymbol()) != null) {
	            int formalCharge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
	            atom.setValency(valencesTable.get(atom.getSymbol()) - formalCharge);
	        }
	        */
	    }
	
	    for (IBond bond : atomContainer.bonds()) {
	        if (allRings.getRings(bond).getAtomContainerCount() > 0) {
	            bond.setFlag(CDKConstants.ISINRING, true);
	        }
	    }
	
	    for (IAtom atom : atomContainer.atoms()) {
	        List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);
	
	        int counter = 0;
	        IAtom any;
	        for (IAtom connectedAtom : connectedAtoms) {
	            any = connectedAtom;
	            if (any.getFlag(CDKConstants.ISINRING)) {
	                counter++;
	            }
	        }
	        atom.setProperty(CDKConstants.RING_CONNECTIONS, counter);
	    }
	
	    // check for aromaticity
	    // a good model for writing SMILES
	    calculateAromaticity( atomContainer );
	    
	    // remove charges
	    for( IAtom at : atomContainer.atoms() ) {
	    	at.setFormalCharge(0);
	    }

	}
	
	public static void calculateImplicitHydrogens( IAtomContainer mol ) {
		
		final String[] HYVALENCE_EL={"C","N","O","S","P", "F"};
        final int[]   HYVALENCE_VAL={ 4,  3,  2,  2,  3,   1 };
        for (int n=0; n < mol.getAtomCount() ;n++)
        {
        	IAtom atom=mol.getAtom(n);
        	String el=atom.getSymbol();
        	int hy=0;
        	
        	// FIXME  I don't care at the moment what this is - I just want the SMILES generators to stop complaining!
        	atom.setImplicitHydrogenCount( 0 );
        	
        	if( el != null ) {
        		for (int i=0;i<HYVALENCE_EL.length;i++) if (el.equals(HYVALENCE_EL[i])) {hy=HYVALENCE_VAL[i]; break;}
        	} else {
        			atom.setImplicitHydrogenCount( 0 );
        	}
        	
        	if (hy==0) continue;
        	int ch = atom.getFormalCharge() == null ? 0 : atom.getFormalCharge();
        	if (el.equals("C")) ch=-Math.abs(ch);
        	final int unpaired=0; // (not current available, maybe introduce later)
        	hy+=ch-unpaired;
        	
			// (needs to include actual H's) for (int i=0;i<bondAdj[n].length;i++) hy-=bondOrder[bondAdj[n][i]];
			for (IBond bond : mol.getConnectedBondsList(atom))
			{
				if (bond.getOrder()==IBond.Order.SINGLE) hy-=1;
				else if (bond.getOrder()==IBond.Order.DOUBLE) hy-=2;
				else if (bond.getOrder()==IBond.Order.TRIPLE) hy-=3;
				else if (bond.getOrder()==IBond.Order.QUADRUPLE) hy-=4;
				// (look for zero-bonds later on)
			}
			
			// treat aromatics as an extra bond (one less hydrogen)
			if( isAromatic(atom) )
				hy--;
			
        	atom.setImplicitHydrogenCount( Math.max(0,hy) );
        }
		
	}
	
	
	public static void calculateAromaticity( IAtomContainer atomContainer ) throws CDKException {
		
		// check for aromaticity
	    // a good model for writing SMILES
	    Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk() ,
	                                              Cycles.all() );
	    try {
	        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
	        //CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
	        aromaticity.apply(atomContainer);
	    } catch (CDKException e) {
	        //logger.debug(e.toString());
	        throw new CDKException(e.toString(), e);
	    }
	    
	}
	
	public static boolean isAromatic( IAtom at, IAtomContainer mol ) {
		
		if( at.getFlag(CDKConstants.ISAROMATIC) )
			return true;
		
		if( at instanceof AromaticSymbolAtom || at instanceof AromaticSymbolQueryAtom || at instanceof AromaticAtom )
			return true;
		
		if( mol != null ) {
			for( IBond b : mol.getConnectedBondsList(at) ) {
				if( b.getFlag(CDKConstants.ISAROMATIC) )
					return true;
				
				if( b instanceof AromaticQueryBond /*|| b instanceof SingleOrAromaticBond || b instanceof AromaticOrSingleQueryBond*/ )
					return true;
			}
		}
		
		return false;
		
	}
	
	public static boolean isAromatic( IBond b ) {
		
		if( b.getFlag(CDKConstants.ISAROMATIC) )
			return true;
		
		return ( b instanceof AromaticQueryBond /*|| b instanceof SingleOrAromaticBond || b instanceof AromaticOrSingleQueryBond*/ );
			
		
		/*if( isAromatic( b.getAtom(0) ) && isAromatic( b.getAtom(1) ) )
			return true;*/
		
		//return false;
		
	}
	
	public static boolean isAromatic( IAtom at ) {
		
		return isAromatic( at, null );
		
	}
	
	
	public static void correctAtomBondTypes( IAtomContainer atomContainer ) {
		
		IBond[] newBonds = new IBond[ atomContainer.getBondCount() ]; 
		
		/*
		// do all ring perception
		AllRingsFinder arf = new AllRingsFinder();
		//HanserRingFinder arf = new HanserRingFinder();
	    IRingSet allRings = null;
	    try {
	    	//arf.getRingSet(atomContainer);
	        //allRings = arf.getRingSet(atomContainer);
	        allRings = arf.findAllRings(atomContainer);
	    } catch (CDKException e) {
	    	e.printStackTrace();
	        //logger.debug(e.toString());
	    }
	    */
		//System.out.println("Correcting atom and bond types");
	    
		// set symbols for missing atoms
		for( IAtom at : atomContainer.atoms() ) {
			if( at instanceof AtomicNumberAtom ) {
				at.setSymbol( PeriodicTable.getSymbol( at.getAtomicNumber() ) );
			} else if( at instanceof SymbolQueryAtom ) {
				at.setAtomicNumber( PeriodicTable.getAtomicNumber( at.getSymbol() ) );
			} else if( at instanceof SMARTSAtom && at.getSymbol() != null) {
				at.setAtomicNumber( PeriodicTable.getAtomicNumber( at.getSymbol() ) );
			} else if( at instanceof SmartsAtomExpression ) {
				SmartsAtomExpression sAt = (SmartsAtomExpression) at;
				List<SmartsExpressionToken> tokens = sAt.tokens;
				
				for( SmartsExpressionToken t : tokens ) {
					if( t.type == 11 || t.type == 2 ) {
						at.setAtomicNumber( t.param );
						at.setSymbol( PeriodicTable.getSymbol( t.param ) );
						break;
					}
					
					if( t.type == 1 ) {
						at.setAtomicNumber( t.param );
						at.setSymbol( PeriodicTable.getSymbol( t.param ) );
						break;
					}
				}
				
			} else if( at instanceof LogicalOperatorAtom ) {
				//at.setSymbol( PeriodicTable.getSymbol( ((LogicalOperatorAtom) at).getLeft().getAtomicNumber() ) );
				at.setSymbol( ((LogicalOperatorAtom) at).getLeft().getSymbol() );
			} else if( at instanceof AromaticSymbolQueryAtom ) {
				at.setFlag(CDKConstants.ISAROMATIC, true);
			} else if( at instanceof AnyAtom ) {
				at.setSymbol("*");
			}
			
			// we're ignoring charges for the moment
			if( at.getFormalCharge() == null )
		    	at.setFormalCharge(0);
		}
		
		// correct to allow for aromaticity and missing bond orders
		int bCount = 0;
		for( IBond b : atomContainer.bonds() ) {
			IBond newBond = b;
			
			
			/*
			if( b instanceof AromaticOrSingleQueryBond || b.getOrder() == null ) {
				b.setOrder(Order.SINGLE);
				//newBond = new OrderQueryBond( (IQueryAtom) b.getAtom(0), (IQueryAtom) b.getAtom(1), b.getOrder(), atomContainer.getBuilder() );
				//atomContainer.removeBond(b);
				//atomContainer.addBond(newBond);
			}
			/*
			if( b instanceof RingBond ) {
				newBond = new AnyOrderQueryBondNotNull( (IQueryAtom) b.getAtom(0), (IQueryAtom) b.getAtom(1), b.getOrder(), atomContainer.getBuilder() );
			}
			*/
			/*
			if( b instanceof AnyOrderQueryBond ) {
				newBond = new AnyOrderQueryBondNotNull( (IQueryAtom) b.getAtom(0), (IQueryAtom) b.getAtom(1), b.getOrder(), atomContainer.getBuilder() );
			}
			*/
			newBonds[bCount] = newBond;
			bCount++;
			
			if( newBond.getOrder() == Order.UNSET || newBond.getOrder() == null ) {
				if( b instanceof AnyOrderQueryBond ) {
					b.setOrder( Order.TRIPLE );
				} else if ( b instanceof DoubleNonAromaticBond ) {
					b.setOrder( Order.DOUBLE );
				} else if ( b instanceof TripleBondAromaticityNotSpecified ) {
					b.setOrder( Order.TRIPLE );
				} else if ( b instanceof CTFileQueryBond ) {
					CTFileQueryBond ctb = (CTFileQueryBond) b;
					
					if( ctb.getType().equals( CTFileQueryBond.Type.DOUBLE ) || ctb.getType().equals( CTFileQueryBond.Type.DOUBLE_OR_AROMATIC ) ) {
						b.setOrder( Order.DOUBLE );
					} else if( ctb.getType().equals( CTFileQueryBond.Type.TRIPLE )  ) {
						b.setOrder( Order.TRIPLE );
					} else {
						b.setOrder( Order.SINGLE );
					}
					
				} else {
					b.setOrder( Order.SINGLE );
				}
			}
			
		    //if (allRings.getRings(b).getAtomContainerCount() > 0) {
		    //        b.setFlag(CDKConstants.ISINRING, true);
		    //}
		}
		
		atomContainer.setBonds(newBonds);
		
		// XXX remove all listeners - not needed, cause concurrency problems
		for( IBond nb : newBonds ) {
			nb.removeListener(atomContainer);
		}
	}
	
	

	
	
	public static int atomSymbolToNumber( IAtom at ) {
		
		/*
		if( symbolToNumber.size() == 0 ) {
			symbolToNumber.put("H", 1);
			symbolToNumber.put("C", 6);
			symbolToNumber.put("N", 7);
			symbolToNumber.put("O", 8);
			symbolToNumber.put("P", 15);
			symbolToNumber.put("S", 16);
			symbolToNumber.put("F", 9);
			symbolToNumber.put("Cl", 17);
			symbolToNumber.put("Br", 35);
			symbolToNumber.put("I", 53);
		}
		
		if( at instanceof AtomicNumberAtom ) {
			return at.getAtomicNumber();
		}
		
		return symbolToNumber.get( at.getSymbol() );
		*/
		
		return PeriodicTable.getAtomicNumber( at.getSymbol() );
	}
	
	/**
	 *  Gets the bondSymbol attribute of the Fingerprinter class
	 *
	 *@param  bond  Description of the Parameter
	 *@return       The bondSymbol value
	 */
	public static String getBondSymbol(IBond bond)
	{
		String bondSymbol = "";
		if (bond.getFlag(CDKConstants.ISAROMATIC))
		{
			bondSymbol = ":";
		} else if (bond.getOrder() == IBond.Order.SINGLE)
		{
			bondSymbol = "-";
		} else if (bond.getOrder() == IBond.Order.DOUBLE)
		{
			bondSymbol = "=";
		} else if (bond.getOrder() == IBond.Order.TRIPLE)
		{
			bondSymbol = "#";
		}
		return bondSymbol;
	}
	
	
	
	public static List<List<Integer>> createBondAdjacencyList( IAtomContainer mol ) {
		
		ArrayList<List<Integer>> adjList = new ArrayList<List<Integer>>( mol.getBondCount() );
		
		for( int a = 0; a < mol.getBondCount(); a++ ) {
			adjList.add(a, new ArrayList<Integer>(8)); // assume max 8 connections ( 2 * 4 )
		}
		
		for( int a = 0; a < mol.getBondCount(); a++ ) {
			for( int b = a; b < mol.getBondCount(); b++ ) {
				if( a != b && mol.getBond(a).isConnectedTo( mol.getBond(b) ) ) {
					
					adjList.get(a).add(b);
					adjList.get(b).add(a);
				}
			}
		}
		
		return adjList;
		
	}
	
	
	/**
	 * Not queue-based - relies on removing neighbours from a set then adding distances, per iteration (until the set is empty)
	 * The set starts as full (i.e all members of the graph)
	 * 
	 * @param adjList
	 * @return
	 */
	public static int[][] bfsShortestPathLengths( List<List<Integer>> adjList ) {
			
			int numBonds = adjList.size();
			
			int[][] pathDistances = new int[ numBonds ][ numBonds ];
			
			
			for( int n = 0; n < adjList.size(); n++ ) {
	
				Set<Integer> elements = new HashSet<Integer>();
				for( int r=0; r < adjList.size(); r++ ) {
					elements.add(r);
				}
				
				// store indices of atoms/bonds instead of the actual objects.  
				// Doing this with bitsets because we want the things uniquified 
				List<Set<Integer>> orderElems = new ArrayList<Set<Integer>>( adjList.size() );
				
				Set<Integer> initialContainer = new HashSet<Integer>();
				initialContainer.add( n );
				orderElems.add( 0, initialContainer );
				
				elements.remove(n);
				
				
				// neighbours of neighbours
				// take atoms of latest radius
				// find neighbouring atoms and bonds
				// remove the atoms and bonds from these neighbours to obtain a radius
				for( int o = 1; ! elements.isEmpty(); o++ ) {
					
					Set<Integer> radialAtomsSet = new HashSet<Integer>();
					
					for (int a : orderElems.get(o-1) ) {
						
						for( Integer radAt : adjList.get(a) ) {
							radialAtomsSet.add( radAt );
							
							
							if( elements.contains( radAt ))
								pathDistances[n][radAt] = o;
							
							elements.remove( radAt );
						}
						
						
					}
					
					
					
					radialAtomsSet.removeAll( orderElems.get( o-1 ) );  // remove all previously used atoms
					
					
					orderElems.add(o, radialAtomsSet);
				}
			}
			
			return pathDistances;
			
		}
	
	/*public static int[][] bfsShortestPathLengths( List<List<Integer>> adjList ) {
		
		int numBonds = adjList.size();
		
		int[][] pathDistances = new int[ numBonds ][ numBonds ];
		
		for( int n=0; n < numBonds; n++ ) {
			
			int distance = 0;
			
			//Set<Integer> discovered = new HashSet<Integer>();
			LinkedList<Integer> queue = new LinkedList<Integer>();
			
			//discovered.add( n );
			queue.addLast( n );
			
			while( ! queue.isEmpty() ) {
				
				int e = queue.removeFirst();
				++distance; 
				
				for( int neighbour : adjList.get(e) ) {
					if( pathDistances[n][neighbour] <= 0 && neighbour != n ) {
						queue.add(neighbour);
						//discovered.add(neighbour);
						
						pathDistances[n][neighbour] = distance;
						//pathDistances[neighbour][n] = distance;
					}
				}
				
			}
		}
		
		return pathDistances;
		
	}*/
	
	
	public static int[][] bondAdjacencyMatrix( IAtomContainer graph ) {
		
		int[][] adjMat = new int[ graph.getBondCount() ][ graph.getBondCount() ];
		
		for( int i = 0; i < graph.getBondCount(); i ++ ) {
			for( int j = 0; j < graph.getBondCount(); j ++ ) {
				
				adjMat[i][j] = 0;
				
			}
		}
		
		for( int i = 0; i < graph.getBondCount(); i ++ ) {
				List<Integer> connected = getConnectedBondsIndices( graph, i );
				
				for( int b : connected ) {
					adjMat[i][b] = 1;
				}

		}
		
		return adjMat;
		
	}
	
	
	public static Map<Integer, Collection<Integer>> adjListToMap( int[][] adjList ) {
	
		Map<Integer, Collection<Integer>> adjMap = new HashMap<Integer, Collection<Integer>>();
	
		for( int i = 0; i < adjList.length; i++ ) {
	
				Set<Integer> set = new HashSet<Integer>();
				for( int a : adjList[i] ) { 
					set.add(a); 
				}
				
				adjMap.put( i, set );
		}
	
		return adjMap;
	}
	
	
	public static int[][] listOfListsToMatrix( Collection<Collection<Integer>> lol ) {
		
		int[][] matrix = new int[ lol.size() ][];
		
		int i = 0;
		for( Collection<Integer> list : lol ) {
			
			
			int[] row = new int[ list.size() ];
			
			int n = 0;
			for( int elem : list ) {
				row[n++] = elem;
			}
			
			matrix[i++] = row;
		}
		
		return matrix;
	}
	
	
	
	public static List<Integer> getConnectedBondsIndices( IAtomContainer graph, int bondIndex ) {
		
		List<Integer> connected = new ArrayList<Integer>(8);
		IBond bond = graph.getBond(bondIndex);
		
		for( IBond b : graph.bonds() ) {
			if( bond.isConnectedTo(b) ) {
				connected.add( graph.getBondNumber(b) );
			}
		}
		
		return connected;
	}
	
	
	public static List<IBond> getConnectedBonds( IAtomContainer graph, int bondIndex ) {
		
		List<IBond> connected = new ArrayList<IBond>(8);
		
		for( int b : getConnectedBondsIndices(graph, bondIndex) ) {
			connected.add( graph.getBond(b) );
		}
		
		return connected;
	}
	
	
	
	/**
    *
    * Returns bond map between source and target molecules based on the atoms
    * @param ac1 source molecule
    * @param ac2 target molecule
    * @param mapping mappings between source and target molecule atoms
    * @return bond map between source and target molecules based on the atoms
    * 
    * XXX NOTE - may need to be synchronized
    */
   public static Map<IBond, IBond> makeBondMapOfAtomMap(IAtomContainer ac1, IAtomContainer ac2, Map<IAtom, IAtom> mapping, boolean matchBondTypes) {
       Map<IBond, IBond> maps = new HashMap<IBond, IBond>();
       
       
       

       for (Map.Entry<IAtom, IAtom> mapS : mapping.entrySet()) {
           //IAtom indexI = getAtomViaID( ac1, mapS.getKey() );
           //IAtom indexJ = getAtomViaID( ac2, mapS.getValue() );
           IAtom indexI = mapS.getKey();
           IAtom indexJ = mapS.getValue();

           for (Map.Entry<IAtom, IAtom> mapD : mapping.entrySet()) {
               //IAtom indexIPlus = getAtomViaID( ac1, mapD.getKey() );
               //IAtom indexJPlus = getAtomViaID( ac2, mapD.getValue() );
               IAtom indexIPlus = mapD.getKey();
               IAtom indexJPlus = mapD.getValue();
               
               if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {
                   IBond ac1Bond = ac1.getBond(indexI, indexIPlus);
                   //if( ac1.contains(indexI) ) {
                	 //  boolean pie = true;
                   //}
                   
                   if (ac1Bond != null) {
                       IBond ac2Bond = ac2.getBond(indexJ, indexJPlus);
                       
                       if (ac2Bond != null) {
                    	   
	                    	// The bond matcher below doesn't take aromaticity into account, evidently
	                    	if( ConvenienceTools.isAromatic(ac1Bond) != ConvenienceTools.isAromatic(ac2Bond) )   
	                    		continue;
	                    	
	                    	//Bond Matcher
	               	    	DefaultBondMatcher bondMatcher = null;
	               	        
	               	        if( ac1Bond instanceof IQueryBond && matchBondTypes )
	               	        	bondMatcher = new DefaultBondMatcher( (IQueryBond) ac1Bond );
	               	        else 
	               	        	bondMatcher = new DefaultBondMatcher( ac1, ac1Bond, matchBondTypes );
	               	        
	               	        
	               	        if( bondMatcher.matches(ac2, ac2Bond) ) {
	               	        	maps.put(ac1Bond, ac2Bond);
	               	        }
	                    	   
	                    	   /*
	                    	   if( matchBondTypes ) {
	                    	        
		                    	   if( ac1Bond instanceof IQueryBond ) {
		                    		   if( ((IQueryBond) ac1Bond).matches(ac2Bond) ) {
		                        		   maps.put(ac1Bond, ac2Bond);
		                        	   }
		                    	   } else {
			                    	   if( ac1Bond.getOrder().equals( ac2Bond.getOrder() ) ) {
			                    		   maps.put(ac1Bond, ac2Bond);
			                    	   }
			                       }
	                    	   } else {
	                    		   maps.put(ac1Bond, ac2Bond);
	                    	   }*/
                    	   
                       }
                   }
               }
           }
       }

  

       return maps;

   }
   
   
   // XXX NOTE - may need to be synchronized
   public static Map<IBond, IBond> makeBondMapOfAtomMap(IAtomContainer ac1, IAtomContainer ac2, Map<IAtom, IAtom> mapping ) {
	   return makeBondMapOfAtomMap( ac1, ac2, mapping, true );
   }
	
	
	
	public static IAtomContainer createCommonSubgraph( IAtomContainer g1, IAtomContainer g2, ArrayList<Integer> match ) {
		HashMap<IAtom, IAtom> atomMap = new HashMap<IAtom, IAtom>();
	   
	   //System.out.println("stuff: " + mapChr.size() + " " + mapper.getQueryMol().getAtomCount() + " " + mapper.getMainMol().getAtomCount() + " | " + mapChr );
	   
	   // translate mapping from array to hash
	   for( int m = 0; m < match.size(); m++ ) {
		   if( match.get(m) >= 0 ) {
			   if( g2.getAtom( m ) != null && g1.getAtom( match.get(m) ) != null )
				   atomMap.put( g2.getAtom( m ), g1.getAtom(  match.get(m) ) );
		   }
	   }
	   
	   Map<IBond, IBond> mapping = makeBondMapOfAtomMap(g2, g1, atomMap);
	   
	   return createCommonSubgraph( g2, g1, mapping );
	}
	
	
	public static IAtomContainer createCommonSubgraph( IAtomContainer g1, IAtomContainer g2, Map<IBond, IBond> mapping ) {
		/*
		if( g1 instanceof IQueryAtomContainer && g1.getNotification() ) {
			for( IAtom a : g1.atoms() ) {
				if( a.getListenerCount() > 1 ) {
				System.err.println("Listeners are evil");
				//Thread.currentThread().getStackTrace();
				Thread.dumpStack();
				//System.exit(1);
				}
			}
		}*/
		
		IAtomContainer common;
		
		if( g1 instanceof IQueryAtomContainer )
			common = new QueryAtomContainer2( DefaultChemObjectBuilder.getInstance() );
		else 
			common = new AtomContainer();
		
		common.setNotification(false);  // hoping this'll remove concurrency problems I've been getting with listener objects
		
		/*
		for( IAtom at : g1.atoms() ) {
			common.addAtom( at );
		}*/
		
		BitSet bondsPresent = new BitSet( g1.getBondCount() );
		
		//int counter = 0;
		for( Entry<IBond, IBond> e : mapping.entrySet() ) {
			IBond b1 = e.getKey();
			Map<Object, Object> propertyUnion = new HashMap<Object, Object>( b1.getProperties() );
			propertyUnion.putAll( e.getValue().getProperties() );
			b1.setProperties(propertyUnion);
			//System.out.println( b1.getProperties() );
			
			
			common.addBond(b1);
			
			if( ! common.contains( b1.getAtom(0) ) ) {
				common.addAtom( b1.getAtom(0) );
				b1.getAtom(0).removeListener( common );
				b1.getAtom(0).removeListener( g1 );
				b1.getAtom(0).removeListener( g2 );
			}
			
			if( ! common.contains( b1.getAtom(1) ) ) {
				common.addAtom( b1.getAtom(1) );
				b1.getAtom(1).removeListener( common );
				b1.getAtom(1).removeListener( g1 );
				b1.getAtom(1).removeListener( g2 );
			}
			
			int g1bn = g1.getBondNumber(b1);
			
			//if( g1bn < 0 )
			//	System.out.println( g1.getBondNumber(b1) );
			
			bondsPresent.set( g1bn );
			//++counter;
		}
		
		common.setProperty(origBondIndicesProperty, bondsPresent);
		
		for( IAtom at : common.atoms() ) {
			at.setNotification(false);
		}
		
		calculateImplicitHydrogens(common);
		
		//System.out.println( "bitset - " + bondsPresent );
		
		return common;
		
	}
	
	
	
	public static IAtomContainer createCommonSubgraphDeepCopy( IAtomContainer g1, IAtomContainer g2, Map<IBond, IBond> mapping ) {
		
		IAtomContainer common;
		
		if( g1 instanceof IQueryAtomContainer )
			common = new QueryAtomContainer( DefaultChemObjectBuilder.getInstance() );
		else 
			common = new AtomContainer();
		
		
		Map<Integer, IAtom> indexToAtom = new HashMap<Integer, IAtom>();
		
		
		
		BitSet bondsPresent = new BitSet( g1.getBondCount() );
		
		//int counter = 0;
		for( Entry<IBond, IBond> e : mapping.entrySet() ) {
			IBond b1 = e.getKey();
			Map<Object, Object> propertyUnion = b1.getProperties();
			propertyUnion.putAll( e.getValue().getProperties() );
			//System.out.println( b1.getProperties() );
			
			
			try {
				
			
				// clone atoms
				IAtom at1 = b1.getAtom(0);
				IAtom at2 = b1.getAtom(1);
				
				if( ! indexToAtom.containsKey( g1.getAtomNumber(at1) ) ) {
					IAtom newAt1;
					newAt1 = at1.clone();
					common.addAtom( newAt1 );
					indexToAtom.put( g1.getAtomNumber(at1), newAt1 );
				}
				
				if( ! indexToAtom.containsKey( g1.getAtomNumber(at2) ) ) {
					IAtom newAt2;
					newAt2 = at2.clone();
					common.addAtom( newAt2 );
					indexToAtom.put( g1.getAtomNumber(at2), newAt2 );
				}
				
				// clone bond
				IBond newBond = b1.clone();
				newBond.setAtom( indexToAtom.get( g1.getAtomNumber(at1) ), 0 );
				newBond.setAtom( indexToAtom.get( g1.getAtomNumber(at2) ), 1 );
				common.addBond(newBond);
			
			} catch (CloneNotSupportedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
			int g1bn = g1.getBondNumber(b1);
			
			 
			bondsPresent.set( g1bn );
		 
		}
		
		common.setProperty(origBondIndicesProperty, bondsPresent);
		//System.out.println( "bitset - " + bondsPresent );
		
		return common;
		
	}
	
	

	/**
	 * Calculates the common subgraphs from both molecules.  This function tests for commonality 
	 * (i.e the degree sequences of the two subgraphs are the same).  If they're not, then
	 * the subgraphs are not common and thus the clique is invalid
	 * 
	 * @param currentClique
	 * @return
	 */
	public static boolean deltaYExchangeOccured( GenerateCompatibilityGraphEdges modProd, IAtomContainer ac1, IAtomContainer ac2, Collection<Integer> clique ) {
		
		
		
		Map<IBond, IBond> commonSubgraphMap = new HashMap<IBond, IBond>( clique.size() );
		Map<IBond, IBond> reverseSubgraphMap = new HashMap<IBond, IBond>( clique.size() );
	
		for( Integer mpn : clique ) {
			int[] node = modProd.getNodes().get(mpn);
			
			IBond hsBond = ac1.getBond(node[0]);
			IBond qBond = ac2.getBond(node[1]);
			
			commonSubgraphMap.put(hsBond, qBond);
		}
		
		//System.out.println( c + " clique size - " + findMCS.bestCliques.get(c).size() );
		
		
		
		for( Entry<IBond, IBond> match : commonSubgraphMap.entrySet() ) {
			reverseSubgraphMap.put( match.getValue() , match.getKey() );
		}
		
		IAtomContainer sGraph = createCommonSubgraph( ac1, ac2, commonSubgraphMap );
		IAtomContainer sGraph2 = createCommonSubgraph( ac2, ac1, reverseSubgraphMap );
		
		
		/*ArrayList<Integer> degrees1 = new ArrayList<Integer>( clique.size() );
		ArrayList<Integer> degrees2 = new ArrayList<Integer>( clique.size() );
		
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
		
		if( degrees1.equals( degrees2 ) )
			return false;*/
		
		
		// Delta-Y exchange normally compares identity of degree sequences
		// However I don't see what's wrong with comparing degree frequencies instead
		int dfSize = Math.min( ac1.getAtomCount(), ac2.getAtomCount() );
		ArrayList<Integer> degreeFreqs1 = new ArrayList<Integer>( dfSize );
		ArrayList<Integer> degreeFreqs2 = new ArrayList<Integer>( dfSize );
		
		for( int n = 0; n <= dfSize; n++ ) {
			degreeFreqs1.add(0);
			degreeFreqs2.add(0);
		}
		
		for( IAtom b : sGraph.atoms() ) {
			int d = sGraph.getConnectedAtomsCount( b );
			degreeFreqs1.set( d, degreeFreqs1.get(d) + 1 );
		}
		
		for( IAtom b : sGraph2.atoms() ) {
			int d = sGraph2.getConnectedAtomsCount( b );
			degreeFreqs2.set( d, degreeFreqs2.get(d) + 1 );
		}
		
		
		if( degreeFreqs1.equals( degreeFreqs2 ) )
			return false;
		
		/*System.out.println( "clique - " + clique );
		for( Integer cInd : clique ) {
			System.out.print( "[" + modProd.getNodes().get(cInd)[0] + "," + modProd.getNodes().get(cInd)[1] + "] " );
		}
		System.out.println();
		
		System.out.println( " 1 degrees - " + degrees1 );
		System.out.println( " 2 degrees - " + degrees2 );*/
		
		/*try {
			System.out.println(" subgraph SMILES - " + sGenerator.create(sGraph) );
			System.out.println(" subgraph SMILES 2 - " + sGenerator.create(sGraph2) );
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
		
		
		return true;
	}
	
	
	/**
	 * As above but when we have the molecule indices, as opposed to the modular product indices
	 * 
	 * @param modProd
	 * @param ac1
	 * @param ac2
	 * @param clique
	 * @return
	 */
	public static boolean deltaYExchangeOccured( IAtomContainer ac1, IAtomContainer ac2, Collection<int[]> clique ) {
		
		
		
		Map<IBond, IBond> commonSubgraphMap = new HashMap<IBond, IBond>( clique.size() );
		Map<IBond, IBond> reverseSubgraphMap = new HashMap<IBond, IBond>( clique.size() );
	
		for( int[] node : clique ) {
			IBond hsBond = ac1.getBond(node[0]);
			IBond qBond = ac2.getBond(node[1]);
			
			commonSubgraphMap.put(hsBond, qBond);
			//System.out.println( node[0] + " " + node[1] );
		}
		
		
		
		for( Entry<IBond, IBond> match : commonSubgraphMap.entrySet() ) {
			reverseSubgraphMap.put( match.getValue() , match.getKey() );
		}
		
		IAtomContainer sGraph = createCommonSubgraph( ac1, ac2, commonSubgraphMap );
		IAtomContainer sGraph2 = createCommonSubgraph( ac2, ac1, reverseSubgraphMap );
		
		
		// Delta-Y exchange normally compares identity of degree sequences
		// However I don't see what's wrong with comparing degree frequencies instead
		int dfSize = Math.min( ac1.getAtomCount(), ac2.getAtomCount() );
		ArrayList<Integer> degreeFreqs1 = new ArrayList<Integer>( dfSize );
		ArrayList<Integer> degreeFreqs2 = new ArrayList<Integer>( dfSize );
		
		for( int n = 0; n <= dfSize; n++ ) {
			degreeFreqs1.add(0);
			degreeFreqs2.add(0);
		}
		
		for( IAtom b : sGraph.atoms() ) {
			int d = sGraph.getConnectedAtomsCount( b );
			degreeFreqs1.set( d, degreeFreqs1.get(d) + 1 );
		}
		
		for( IAtom b : sGraph2.atoms() ) {
			int d = sGraph2.getConnectedAtomsCount( b );
			degreeFreqs2.set( d, degreeFreqs2.get(d) + 1 );
		}
		
		
		if( degreeFreqs1.equals( degreeFreqs2 ) )
			return false;
	
		
		
		return true;
	}
	
	
	public static boolean bondsMatch( IAtomContainer queryMol, IAtomContainer targetMol, IBond qb, IBond tb, boolean compareBondAtoms ) {
		
		//boolean matches = true;
		
		DefaultBondMatcher dbf = null;
		
		if( qb instanceof SmartsBondExpression ) {
			SmartsBondExpression sbe = (SmartsBondExpression) qb;
		
			if( sbe.tokens.size() > 0 )
				dbf = new DefaultBondMatcher( (SmartsBondExpression) qb);
			else
				dbf = new DefaultBondMatcher(queryMol, qb, true);
		} else { 
			dbf = new DefaultBondMatcher(queryMol, qb, true);
		}
		
		return dbf.matches(targetMol, tb);
		/*
		if( bond1 instanceof IQueryBond ) {
			IQueryBond qBond1 = (IQueryBond) bond1; 
			
			matches = qBond1.matches( bond2 );
		} else {
			matches = bond1.getOrder() == bond2.getOrder();
		}
		
		if( compareBondAtoms ) {
			
			
			
		}
		
		return matches;
		*/
		
	}
	
	
	public static boolean atomsMatch( IAtom at1, IAtom at2 ) {
		
		//boolean matches = true;
		
		DefaultAtomMatcher daf = new DefaultAtomMatcher(at1, false);
		
		return daf.matches(at2);
		
		
	}
	


	/**
	 * Calculates ring information using the latest CDK ring perception algorithms (2015).
	 * 
	 * This was designed to work with the MCS algorithms I produced.  As of such, in addition to marking whether
	 * an atom/bond is in a ring, it determines which rings an atom/bond are in, the number of rings per molecule,
	 * and which rings contain which bonds.
	 * 
	 * Thanks to John May for advice on how to use the new Ring Perception classes
	 * 
	 * XXX NOTE - may need to be synchronized
	 * 
	 * @author Edmund Duesbury
	 * @param atomContainer
	 * @return
	 */
	public static int countRings(IAtomContainer atomContainer) {

		int minCycleBasisCount = 0;
		int sum = 0;
		// ring search finds the biconnected components (i.e. ring systems) and separates
		// out isolated/spiro rings from fused/bridged systems
		RingSearch ringSearch = new RingSearch(atomContainer);

		// for finding fused bonds
		CycleFinder cf = Cycles.mcb();
		//int maxPathLen = atomContainer.getBondCount();

		EdgeToBondMap bondMap = EdgeToBondMap.withSpaceFor(atomContainer);
		int[][] adjList = GraphUtil.toAdjList(atomContainer, bondMap);

		ArrayList<Integer> bondRingCounts = new ArrayList<Integer>( atomContainer.getBondCount() );
		for( int b = 0; b < atomContainer.getBondCount(); b++ ) {
			bondRingCounts.add( 0 );
		}

		// determine number of bond which would need to be broken to make fused system acyclic
		List<IAtomContainer> frFragments = ringSearch.fusedRingFragments() ;
		for (IAtomContainer m : frFragments ) {
			sum += m.getBondCount() - (m.getAtomCount() - 1);
		}


		for (IAtom atom : atomContainer.atoms()) {
			atom.setFlag(CDKConstants.ISINRING, false);
		}

		int[][] isolatedRings = ringSearch.isolated();
		int[][] fusedRings = ringSearch.fused();
		
		ArrayList<ArrayList<Integer>> ringBondsList = new ArrayList<ArrayList<Integer>>( sum );
		
		minCycleBasisCount = configureMolRingInfo(atomContainer, isolatedRings, cf, adjList, minCycleBasisCount, bondMap, bondRingCounts, ringBondsList);
		minCycleBasisCount = configureMolRingInfo(atomContainer, fusedRings, cf, adjList, minCycleBasisCount, bondMap, bondRingCounts, ringBondsList);
				
		sum += ringSearch.isolated().length;

		


/*
		for( int[] rSystem : isolatedRings ) {

			HashSet<Integer> ringBonds = new HashSet<Integer>( rSystem.length );

			// configure ring settings for all atoms that exist in a ring
			for( int rAtom = 0; rAtom < rSystem.length; rAtom++ ) {
				IAtom atom = atomContainer.getAtom( rSystem[rAtom] );

				List<Integer> ringsizes = null;
				if( atom.getProperty( CDKConstants.RING_SIZES ) != null )
					ringsizes = (List<Integer>) atom.getProperty( CDKConstants.RING_SIZES );
				else
					ringsizes = new ArrayList<Integer>();

				ringsizes.add( rSystem.length );

				atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
				atom.setFlag(CDKConstants.ISINRING, true);

				// do same for ring bonds
				for( int otherAt = 1; otherAt < rSystem.length; otherAt++ ) {

					IBond bond = atomContainer.getBond( 
							atomContainer.getAtom(rSystem[rAtom]) ,
							atomContainer.getAtom(rSystem[otherAt])
							);

					if( bond != null ) {
						bond.setFlag(CDKConstants.ISINRING, true);
						ringBonds.add( atomContainer.getBondNumber(bond) );

						if( bond.getProperty( CDKConstants.RING_CONNECTIONS ) == null )
							bond.setProperty( CDKConstants.RING_CONNECTIONS, new ArrayList<Integer>() );

						ArrayList<Integer> ringsOfBond = ((ArrayList<Integer>) bond.getProperty( CDKConstants.RING_CONNECTIONS ));

						if( ! ringsOfBond.contains( minCycleBasisCount ) )
							ringsOfBond.add( minCycleBasisCount );

						//System.err.println( "ring bond set - " + bond );
					}
				}


			}

			minCycleBasisCount++;
			ringBondsList.add( new ArrayList<Integer>( ringBonds ) );
		}


		

		for( int[] fSystem : fusedRings ) {

			HashSet<Integer> ringBonds = new HashSet<Integer>( fSystem.length );

			// take the subgraph and find the MCB
			int[][] subAdjList = GraphUtil.subgraph(adjList, fSystem);
			int[][] paths = null;
			try {
				paths = cf.find(atomContainer, subAdjList, maxPathLen).paths();
			} catch (Intractable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			// for each ring in this system...
			for (int[] path : paths) {
				IBond[] bondPath = new IBond[path.length - 1];

				// map each bond back to the parent
				for (int i = 0, last = path.length - 1; i < last; i++) {
					bondPath[i] = bondMap.get(fSystem[path[i]], fSystem[path[i+1]]);
					Integer bondIndex = atomContainer.getBondNumber(bondPath[i]);
					bondRingCounts.set( bondIndex, bondRingCounts.get(bondIndex) + 1 );

					if( bondPath[i].getProperty( CDKConstants.RING_CONNECTIONS ) == null )
						bondPath[i].setProperty( CDKConstants.RING_CONNECTIONS, new ArrayList<Integer>() );

					((ArrayList<Integer>) bondPath[i].getProperty( CDKConstants.RING_CONNECTIONS )).add( minCycleBasisCount );

					bondPath[i].setFlag(CDKConstants.ISINRING, true);
					ringBonds.add( bondIndex );
					//System.out.print( bondIndex + ", " +  bondRingCounts.get(bondIndex) + "| " );
				}
				//System.out.println(bondPath.length + " ");

				minCycleBasisCount++;
			}

			ringBondsList.add( new ArrayList<Integer>( ringBonds ) );
			//System.out.println();

			// configure ring settings for all atoms that exist in a ring
			for( int rAtom = 0; rAtom < fSystem.length; rAtom++ ) {
				IAtom atom = atomContainer.getAtom( fSystem[rAtom] );

				List<Integer> ringsizes = null;
				if( atom.getProperty( CDKConstants.RING_SIZES ) != null )
					ringsizes = (List<Integer>) atom.getProperty( CDKConstants.RING_SIZES );
				else
					ringsizes = new ArrayList<Integer>();

				ringsizes.add( fSystem.length );

				atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
				atom.setFlag(CDKConstants.ISINRING, true);

				 
			}


		}
*/

		atomContainer.setProperty( CDKConstants.RELEVANT_RINGS, ringBondsList );

	


		
		// number of rings per atom
		for (IAtom atom : atomContainer.atoms()) {
			List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);

			int counter = 0;
			IAtom any;
			for (IAtom connectedAtom : connectedAtoms) {
				any = connectedAtom;
				if (any.getFlag(CDKConstants.ISINRING)) {
					counter++;
				}
			}
			atom.setProperty(CDKConstants.RING_CONNECTIONS, counter);
		}


		ArrayList<Integer> fusedBonds = new ArrayList<Integer>( atomContainer.getBondCount() ); // find fused ring bonds

		// Assign number of rings each bond is in
		for( int b = 0; b < atomContainer.getBondCount(); b++ ) {

			if( bondRingCounts.get(b) > 1 )
				fusedBonds.add(b);
		}

		atomContainer.setProperty(CDKConstants.RING_CONNECTIONS , fusedBonds );
		
		
		// XXX remove references to fragments to destroy listeners
		for( IAtom at : atomContainer.atoms() ) {
			for( IAtomContainer rf : frFragments )
				at.removeListener(rf);
		}

		return sum;
	}
	
	
	
	
	private static int configureMolRingInfo( IAtomContainer atomContainer, int[][] ringPaths, CycleFinder cf, int[][] adjList, int startPoint, EdgeToBondMap bondMap, ArrayList<Integer> bondRingCounts, ArrayList<ArrayList<Integer>> ringBondsList ) {
		
		int minCycleBasisCount = startPoint;
		
		for( int[] fSystem : ringPaths ) {

			HashSet<Integer> ringBonds = new HashSet<Integer>( fSystem.length );
			
			// take the subgraph and find the MCB
			int[][] subAdjList = GraphUtil.subgraph(adjList, fSystem);
			//System.out.println( new DenseDoubleMatrix2D( int2DtoDouble2D(subAdjList) ) );
			int[][] paths = null;
			try {
				paths = cf.find(atomContainer, subAdjList, atomContainer.getBondCount() ).paths();
			} catch (Intractable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			// for each ring in this system...
			for (int[] path : paths) {
				IBond[] bondPath = new IBond[path.length - 1];

				// map each bond back to the parent
				for (int i = 0, last = path.length - 1; i < last; i++) {
					bondPath[i] = bondMap.get(fSystem[path[i]], fSystem[path[i+1]]);
					Integer bondIndex = atomContainer.getBondNumber(bondPath[i]);
					bondRingCounts.set( bondIndex, bondRingCounts.get(bondIndex) + 1 );

					if( bondPath[i].getProperty( CDKConstants.RING_CONNECTIONS ) == null )
						bondPath[i].setProperty( CDKConstants.RING_CONNECTIONS, new ArrayList<Integer>() );

					((ArrayList<Integer>) bondPath[i].getProperty( CDKConstants.RING_CONNECTIONS )).add( minCycleBasisCount );

					bondPath[i].setFlag(CDKConstants.ISINRING, true);
					ringBonds.add( bondIndex );
					//System.out.print( bondIndex + ", " +  bondRingCounts.get(bondIndex) + "| " );
				}
				//System.out.println(bondPath.length + " ");

				minCycleBasisCount++;
			}

			ringBondsList.add( new ArrayList<Integer>( ringBonds ) );
			//System.out.println();

			// configure ring settings for all atoms that exist in a ring
			for( int rAtom = 0; rAtom < fSystem.length; rAtom++ ) {
				IAtom atom = atomContainer.getAtom( fSystem[rAtom] );

				List<Integer> ringsizes = null;
				if( atom.getProperty( CDKConstants.RING_SIZES ) != null )
					ringsizes = (List<Integer>) atom.getProperty( CDKConstants.RING_SIZES );
				else
					ringsizes = new ArrayList<Integer>();

				ringsizes.add( fSystem.length );

				atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
				atom.setFlag(CDKConstants.ISINRING, true);

				 
			}


		}
		
		
		return minCycleBasisCount;
	}
	
	
	
	public static boolean isRingBond( IBond bond ) {
	
		
		
		
		/*if( bond instanceof SmartsBondExpression ) {
			SmartsBondExpression sBond = (SmartsBondExpression) bond;
			
			for (int i = 0; i< sBond.tokens.size(); i++)
			{
				int tok = sBond.tokens.get(i).intValue();
				
				if( i > 0 && sBond.tokens.get(i-1).intValue() == (SmartsConst.LO + SmartsConst.LO_NOT) && tok <= SmartsConst.BT_RING )
					return false;
				else if (tok == SmartsConst.BT_RING )
					return true;
			}
		}*/
	
		
		if( bond.getFlag( CDKConstants.ISINRING ) == true  )
			return true;
		
		if( bond.getFlag( CDKConstants.ISINRING ) == false  ) 
			return false;
		
		if( bond instanceof RingBond || bond instanceof RingQueryBond || bond instanceof AromaticQueryBond || bond instanceof AromaticOrSingleQueryBond )
			return true;
		
		
		
		return false;
	}
	
	
	
	
	/**
	 * Find minimum cost pairs of a matrix (Hungarian Algorithm/Linear assignment)
	 * 
	 * @param simMat
	 * @return
	 */
	public static List<int[]> linearAssignment( double[][] mat ) {
		
		List<int[]> matchList = new ArrayList<int[]>( mat.length );
		
		HungarianAlgorithm ha = new HungarianAlgorithm( mat );
		int[] matches = ha.execute();
		
		// obtain the best pairs - convert array to ArrayList of pairs
		for( int n = 0; n < matches.length; n++ ) {
			
			//if( matches[n] <= 0 )
			//	continue;
			
			//int[] matchPair = new int[]{ n, matches[n] };
			int[] matchPair = new int[]{ n, matches[n] };
			matchList.add( matchPair );
		}
		
		//System.out.println( bigraphMatrix );
		
		
		return matchList;
	}

	
	
	/**
	 * Morgan algorithm-like radial gathering of atoms (including any bonds in the path of expansion), node-induced
	 * 
	 * @param mol
	 * @param atIndex centre of desired subgraph to expand out from
	 * @param order		radius (order of 1 means no expansion), so specify 2 to get a radius of 2 atoms (including central atom)
	 * @return radial subgraph
	 */
	public static IAtomContainer getNeighbourhoodGraph( IAtomContainer mol, int atIndex, int order ) {
		
		IAtomContainer subgraph = new AtomContainer();

		
		LinkedList<IAtom> neighbourhood = new LinkedList<IAtom>();
		LinkedList<IAtom> tNeighbourhood = new LinkedList<IAtom>();
		
		// add initial atom
		neighbourhood.add( mol.getAtom(atIndex) );
		tNeighbourhood.add( mol.getAtom(atIndex) );
		
		//List<IAtom> atomNeighbours = mol.getConnectedAtomsList( mol.getAtom(atIndex) );
		
		// add initial neighbours
		neighbourhood.addAll( mol.getConnectedAtomsList( mol.getAtom(atIndex) ) );
		tNeighbourhood.addAll( mol.getConnectedAtomsList( mol.getAtom(atIndex) ) );
		
		
		// neighbours of neighbours
		for( int o = 0; o < order - 1; o++ ) {
			
			LinkedList<IAtom> tNeighbourhood2 = new LinkedList<IAtom>();
			
			for( int a = 0; a < tNeighbourhood.size(); a++ ) {
				tNeighbourhood2.addAll( mol.getConnectedAtomsList( tNeighbourhood.get(a) ) );
			} 
			
			tNeighbourhood.clear();
			
			tNeighbourhood.addAll( tNeighbourhood2 );
			neighbourhood.addAll( tNeighbourhood );
			
		}
		
		IAtom[] uniqueNeighbours = new HashSet<IAtom>( neighbourhood ).toArray(new IAtom[0]);
		
		for( int n = 0; n < uniqueNeighbours.length; n++ ) {
			subgraph.addAtom( uniqueNeighbours[n] );
		}
		
		
		for( int i = 0; i < uniqueNeighbours.length; i++ ) {
			for( int j = 0; j < uniqueNeighbours.length; j++ ) {
				
				if( i != j ) {
					IBond b = mol.getBond( uniqueNeighbours[i], uniqueNeighbours[j] );
					
					if( b != null && subgraph.getBond( b.getAtom(0), b.getAtom(1) ) == null ){
						subgraph.addBond( subgraph.getAtomNumber( b.getAtom(0) ), subgraph.getAtomNumber( b.getAtom(1) ), b.getOrder() );
					}
				}
				
			}
		}

		return subgraph;
	}
	
	
	
	/**
	 * Morgan algorithm-like radial gathering of atoms (including any bonds in the path of expansion), edge-induced
	 * 
	 * - starts from central atom
	 * - gets connected bonds, adds non-used ones to subgraph
	 * - same with atoms
	 * - the most recently added atoms are then expanded from in a similar manner.
	 * 
	 * @param mol
	 * @param atIndex centre of desired subgraph to expand out from
	 * @param order		radius (order of 1 means no expansion), so specify 2 to get a radius of 2 atoms (including central atom)
	 * @return radial subgraph
	 */
	public static IAtomContainer[] getNeighbourhoodGraphEdgeInduced( IAtomContainer mol, int atIndex, int order ) {
		
		IAtomContainer[] subgraphs = new IAtomContainer[ order ];
		
		IAtomContainer subgraph;
		
		if( mol instanceof IQueryAtomContainer )
			subgraph = new QueryAtomContainer(null);
		else
			subgraph = new AtomContainer();
		
		//LinkedList<IAtom> neighbourhood = new LinkedList<IAtom>();

		// store indices of atoms/bonds instead of the actual objects.  
		// Doing this with bitsets because we want the things uniquified 
		ArrayList<BitSet> orderAtoms = new ArrayList<BitSet>( order );
		ArrayList<BitSet> orderBonds = new ArrayList<BitSet>( order );
		
		BitSet initialContainer = new BitSet();
		initialContainer.set( atIndex );
		orderAtoms.add(0, initialContainer );
		orderBonds.add(0, new BitSet() );
		
		
		
		

		BitSet atomsPresent = new BitSet( mol.getBondCount() );  // record of all atoms in the subgraph
		BitSet bondsPresent = new BitSet( mol.getBondCount() );  // record of all bonds in the subgraph
		
		
		// add initial atom
		subgraph.addAtom( mol.getAtom(atIndex) );
		atomsPresent.set( atIndex );
		
		try {
			subgraphs[0] = subgraph.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// neighbours of neighbours
		// take atoms of latest radius
		// find neighbouring atoms and bonds
		// remove the atoms and bonds from these neighbours to obtain a radius
		for( int o = 1; o < order; o++ ) {
			
			BitSet radialAtomsSet = new BitSet();
			BitSet radialBondsSet = new BitSet();
			
			for (int a = orderAtoms.get(o-1).nextSetBit(0); a >= 0; a = orderAtoms.get(o-1).nextSetBit(a + 1) ) {
				IAtom at = mol.getAtom(a);
				
				for( IAtom radAt : mol.getConnectedAtomsList( at ) ) {
					radialAtomsSet.set( mol.getAtomNumber(radAt) );
					atomsPresent.set( mol.getAtomNumber(radAt) );
				}
				
				for( IBond radBond : mol.getConnectedBondsList( at ) ) {
					radialBondsSet.set( mol.getBondNumber(radBond) );
					bondsPresent.set( mol.getBondNumber(radBond) );
				}
				
				
			}
			
			
			radialAtomsSet.andNot( orderAtoms.get( o-1 ) );  // remove all previously used atoms
			radialBondsSet.andNot( orderBonds.get( o-1 ) );  // remove all previously used bonds
			/*for( int n = 0; n < o-1; n++ ) {
				radialBondsSet.andNot( orderBonds.get( n ) );  // remove all previously used bonds
			}*/
			
			orderAtoms.add(o, radialAtomsSet);
			orderBonds.add(o, radialBondsSet);
			
			
			if( mol instanceof IQueryAtomContainer )
				subgraph = new QueryAtomContainer(null);
			else
				subgraph = new AtomContainer();
			
			
			for (int a = atomsPresent.nextSetBit(0); a >= 0; a = atomsPresent.nextSetBit(a + 1) ) {
				subgraph.addAtom( mol.getAtom(a) );
			}
			
			for (int b = bondsPresent.nextSetBit(0); b >= 0; b = bondsPresent.nextSetBit(b + 1) ) {
				IBond nBond = mol.getBond( b );
				subgraph.addBond( nBond );
				
				//System.out.print( b + "  " + mol.getBond( b ).getProperty(CDKSMARTSHyperstructureFitness.bondMolOriginType) + "  ");
			}
			
			//System.out.println( "subgraph bonds - " + subgraph.getBondCount() );
			subgraph.setProperty(origBondIndicesProperty, bondsPresent.hashCode() );  // uniquification purposes
			String s = CDKSMARTSHyperstructureFitness.bondMolOriginType;
			subgraphs[o] = subgraph;
		}
		

		return subgraphs;
	}
	
	
	public static boolean isGhostSubstructure( IAtomContainer graph ) {
		
		List<Integer> bonds = new ArrayList<Integer>( 
				(List<Integer>) graph.getBond(0).getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType ) 
		);
		
		for( int b = 1; b < graph.getBondCount(); b++ ) {
			bonds.retainAll( (List<Integer>) graph.getBond(b).getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType ) );
			
			if( bonds.size() == 0 )
				return true;
		}
		
		return false;
	}
	
	
	/*
	 * is the query a subgraph of the target molecule?
	 */
	public static boolean isSubgraph( IAtomContainer query, IAtomContainer target ) {
		
			//Pattern pattern = VentoFoggia.findSubstructure(query);
			org.cisrg.mapping.VentoFoggia.Pattern pattern = VentoFoggia2.findSubstructure(query, false);
			return pattern.matches(target);
		
	}
	
	
	/**
	 * Takes atom map, and bonds from g1, and finds the corresponding bonds (from the atom map) in g2.
	 * 
	 * loop through mapped bonds in g1
	 *   loop through all bonds in g2
	 *     if g2 bond has both atoms mapped to g1 bond, note the correspondance
	 *     
	 * FIXME  This doesn't always get all the desired bonds for a given atom map, though so far it seems to work for the majority of cases (as detailed in my thesis)    
	 * 
	 * @param atomMap
	 * @param g2
	 * @param bonds
	 * @author edmund duesbury
	 */
	public static Map<IBond, IBond> bondMapFromOtherGraph( Map<IAtom, IAtom> atomMap, IAtomContainer g1, IAtomContainer g2, List<Integer> mappedBonds, boolean swap ) {
		
		Map<IBond, IBond> commonBonds = new HashMap<IBond, IBond>();
		

		
		if( swap ) {  // that is, the mapped bonds are from graph 2, not graph 1
			
			// need a reversed atom map from g2 to g1 (instead of the current g1 to g2 atoms)
			Map<IAtom, IAtom> atomMap2 = new HashMap<IAtom, IAtom>();
			for( Entry<IAtom, IAtom> e : atomMap.entrySet() ) {
				atomMap2.put( e.getValue() , e.getKey() );
			}
			
			for( Integer mbIndex : mappedBonds ) {
				IBond mappedBond = g2.getBond(mbIndex);
				
				IAtom g1At1 = atomMap2.get( mappedBond.getAtom(0) );
				IAtom g1At2 = atomMap2.get( mappedBond.getAtom(1) );
				
				for( IBond g1Bond : g1.bonds() ) {

					if( g1Bond.contains(g1At1) && g1Bond.contains(g1At2) ) {
						commonBonds.put( g1Bond, mappedBond );
					}
				}
			}
		} else {
			for( Integer mbIndex : mappedBonds ) {
				IBond mappedBond = g1.getBond(mbIndex);
				IAtom mappedAtom1 = mappedBond.getAtom(0);
				IAtom mappedAtom2 = mappedBond.getAtom(1);
				
				IAtom g2At1 = atomMap.get( mappedAtom1 );
				IAtom g2At2 = atomMap.get( mappedAtom2 );
				
				if( g2At1 == null || g2At2 == null )
					continue;
				
				IBond g2BondInferred = g2.getBond(g2At1, g2At2);
				
				if( g2BondInferred != null )
					commonBonds.put( mappedBond, g2BondInferred );
				
				/*for( IBond g2Bond : g2.bonds() ) {
					
					if( g2Bond.contains(g2At1) && g2Bond.contains(g2At2) ) {
						commonBonds.put( mappedBond, g2Bond );
					}
				}*/
			}
		}
		
		return commonBonds;
	}
	
	
	//map must be a bijection in order for this to work properly
	public static <K,V> HashMap<V,K> reverse(Map<K,V> map) {
	    HashMap<V,K> rev = new HashMap<V, K>();
	    for(Map.Entry<K,V> entry : map.entrySet())
	        rev.put(entry.getValue(), entry.getKey());
	    return rev;
	}
	

	
	 
	
	public static double calculateTversky(BitSet bitset1, BitSet bitset2, double alpha, double beta) throws CDKException
    {
		double _bitset1_cardinality = bitset1.cardinality();
        double _bitset2_cardinality = bitset2.cardinality();
        
       
        
        if (bitset1.size() != bitset2.size()) {
            throw new CDKException("Bitsets must have the same bit length");
        }
        BitSet one_and_two = (BitSet)bitset1.clone();
        one_and_two.and(bitset2);
        double _common_bit_count = one_and_two.cardinality();
       
        //similarity = (matchCount) / ( (1 - alpha) * (rBondCount - matchCount) + (beta * ( pBondCount - matchCount )) + matchCount  );
        
        
        double similarity = _common_bit_count/
        		( 
        				alpha * ( _bitset1_cardinality - _common_bit_count ) + 
        				beta  * ( _bitset2_cardinality - _common_bit_count ) +
        				_common_bit_count
        		);
        
        System.out.println( "cardins - " + _bitset1_cardinality + " " + _bitset2_cardinality + " , common - " + _common_bit_count + " | similarity = " + similarity );
        
        return similarity;
        //return _common_bit_count/(_bitset1_cardinality + _bitset2_cardinality - _common_bit_count);
    }
	
	
	public static double calculateTversky(double matchCount, double rBondCount, double pBondCount, double alpha, double beta  ) throws IOException {
	       int decimalPlaces = 4;
	       //double rBondCount = 0;
	       //double pBondCount = 0;
	       
	       //double alpha = 0.1;  
	       //double beta = 1.0 - alpha;
	       double similarity = 0.0;

	           //System.out.println( matchCount + " / " + alpha * (rBondCount - matchCount) + " + " + (beta * ( pBondCount - matchCount )) + " + " + matchCount );
	       similarity = (matchCount) / ( alpha * (rBondCount - matchCount) + (beta * ( pBondCount - matchCount )) + matchCount  );
	       BigDecimal tan = new BigDecimal(similarity);
	       tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
	       similarity = tan.doubleValue();
	       
	       return similarity;
	   }
	
	/*
	 * Try this equation (not quite Tversky but still asymmetric):
	 * 
	 * c / (c-a) + (c-b) + a + b.  a = number of database molecule bonds.  b = weighted sum of hyperstructure bonds
	 */
	public static double calculateWeightedTversky(double commonWeighted, double rBondCount, double pBondCount, double alpha, double beta  ) throws IOException {
	       int decimalPlaces = 4;
	       
	       double similarity = 0.0;

	       if( commonWeighted < 1 )  // pointless calculation if there're no bonds in common
	    	   return similarity;
	       
	           //System.out.println( matchCount + " / " + alpha * (rBondCount - matchCount) + " + " + (beta * ( pBondCount - matchCount )) + " + " + matchCount );
	       similarity = (commonWeighted) / ( alpha * (commonWeighted - rBondCount) + (beta * (commonWeighted - pBondCount)) + commonWeighted );
	      
	       try {
			BigDecimal tan = new BigDecimal(similarity);
			   tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
			   similarity = tan.doubleValue();
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println( "error for this similarity, returning 0.0.  Equation: " + commonWeighted + " / " + alpha * (rBondCount - commonWeighted) + " + " + (beta * ( pBondCount - commonWeighted )) + " + " + commonWeighted );
				return 0.0;
			}
	       
	       return similarity;
	   }
	
	//private static HashMap<String, Integer> symbolToNumber = new HashMap<String, Integer>();
	public static double[][] int2DtoDouble2D( int[][] intArray ) {
		
		double[][] dArray = new double[ intArray.length ][ intArray[0].length ];
		
	
	    for(int i = 0; i < intArray.length; i++)
	    {
	        for(int j = 0; j < intArray[0].length; j++)
	        	dArray[i][j] = (double) intArray[i][j];
	    }
		
	    return dArray;
	}
	
	
	public static IAtomContainer createSubgraph( IAtomContainer sourceMol, Collection<Integer> bonds ) {
			
			IAtomContainer nMol = new AtomContainer();
			for( Integer bi : bonds ) {
				IBond b = sourceMol.getBond(bi);
				nMol.addBond( b );
				nMol.addAtom( b.getAtom(0) );
				nMol.addAtom( b.getAtom(1) );
			}
			
			return nMol;
	}


	public static Collection<Integer> createCollection( Collection<Integer> c ) {
		return new BitSetExtended<Integer>( c );
		//return new BitSetCollectionAL<Integer>( c );
		//return new HashSet<Integer>( c );
	}

	public static Collection<Integer> createCollection( int s ) {
		return new BitSetExtended<Integer>( s );
		//return new HashSet<Integer>( s );
	}
	
	
	// for MaxCliqueSeq C++ program
	public static void writeDIMACSGraph( GenerateCompatibilityGraphEdges mp ) {
			
			List<Collection<Integer>> edges = mp.getAdjacencyList();
			
			Writer writer = null;
	
			try {
			    writer = new BufferedWriter(new OutputStreamWriter(
			          new FileOutputStream("/opt/source_MaxCliquePara_v2.2/DIMACS_subset/tempGraph4.clq"), "utf-8"));
			    
			    writer.write("p edge " + edges.size() + " " + mp.numberOfEdges + "\n" );
			    
			    for( int i = 0; i < edges.size(); i++ ) {
			    	for( Integer j : edges.get(i) ) {
			    		writer.write("e " + (i+1) + " " + (j+1) + "\n" );
			    	}
			    }
			} catch (IOException ex) {
			  // report
			} finally {
			   try {writer.close();} catch (Exception ex) {}
			} 
			
		}

	
	public static <T> void printList( List<T> list ) {
		
		/*for( int[] m : list ) {
			System.out.print("[");
			for( int e : m ) {
				System.out.print( e + ", " );
			}
			System.out.print("], ");
		}*/
		if( list.isEmpty() ) {
			System.out.println("no list");
			return;
		}
		
		for( T e1 : list ) {
			if( e1 instanceof List ) {
				printList( (List<T>) e1 );
			} else if ( e1 instanceof int[] ) {
				System.out.print("[");
				for( int e2 : (int[]) e1 ) {
					System.out.print( e2 + ", " );
				}
				System.out.print("], ");
			} else {
				System.out.println("no list");
			}
		}
		
		System.out.println();
	}
	
	
	
	
	public static void printHyperstructureStats( IAtomContainer hs, IAtomContainer[] queries ) {
		
		int queryAtoms = 0;
		int queryBonds = 0;
		
		if( queries != null ) {
			for( int q = 0; q < queries.length; q++ ) {
				queryAtoms += queries[q].getAtomCount();
				queryBonds += queries[q].getBondCount();
			}
		}
		
		System.out.println( "Number of hyperstructure atoms: " + hs.getAtomCount() );
		System.out.println( "Number of hyperstructure bonds: " + hs.getBondCount() );
		System.out.println( "Atom compression ratio: " + (float) queryAtoms / hs.getAtomCount() );
		System.out.println( "Bond compression ratio: " + (float) queryBonds / hs.getBondCount() );
	}
	
	
	public static String origBondIndicesProperty = "_pBonds";
	public static String atomTypeProperty = "_at";


	


	




	
	
}
