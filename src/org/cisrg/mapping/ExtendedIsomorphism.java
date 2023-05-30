/* 
 * Copyright (C) 2009-2014  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.cisrg.mapping;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;

import org.cisrg.mapping.ChemAxonMCS.ChemAxonMCSOptions;
//import org.openscience.cdk.annotations.TestClass;
//import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.smsd.Isomorphism;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
//import org.openscience.smsd.BaseMapping;
//import org.openscience.smsd.Substructure;
//import org.openscience.smsd.algorithm.mcsplus.MCSPlusHandler;
//import org.openscience.smsd.algorithm.rgraph.CDKMCSHandler;
import org.openscience.smsd.algorithm.single.SingleMappingHandler;
//import org.openscience.smsd.algorithm.vflib.VF2MCS;
//import org.openscience.smsd.helper.MoleculeInitializer;
import org.openscience.smsd.interfaces.Algorithm;
//import static org.openscience.smsd.interfaces.Algorithm.CDKMCS;
//import static org.openscience.smsd.interfaces.Algorithm.DEFAULT;
//import static org.openscience.smsd.interfaces.Algorithm.MCSPlus;
//import static org.openscience.smsd.interfaces.Algorithm.VFLibMCS;

/**
 * <p>
 * This class implements the Isomorphism- a multipurpose structure comparison
 * tool. It allows users to, i) find the maximal common substructure(s) (MCS);
 * ii) perform the mapping of a substructure in another structure, and; iii) map
 * two isomorphic structures.</p>
 *
 * <p>
 * It also comes with various published algorithms. The user is free to choose
 * his favorite algorithm to perform MCS or substructure search. For
 * example:</p> <OL> <lI>0: Default, <lI>1: MCSPlus, <lI>2: VFLibMCS, <lI>3:
 * CDKMCS </OL>
 * <p>
 * It also has a set of robust chemical filters (i.e. bond energy, fragment
 * count, stereo & bond match) to sort the reported MCS solutions in a
 * chemically relevant manner. Each comparison can be made with or without using
 * the bond sensitive mode and with implicit or explicit hydrogens.</p>
 *
 * <p>
 * If you are using <font color="#FF0000">Isomorphism, please cite Rahman
 * <i>et.al. 2009</i></font> {
 *
 * @cdk.cite SMSD2009}. The Isomorphism algorithm is described in this paper.
 * </p>
 *
 * <p>
 * An example for <b>MCS search</b>:</p> <font color="#003366">  <pre>
 *
 *
 * SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 * // Benzene
 * IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
 * // Napthalene
 * IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
 * //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
 * //Algorithm is VF2MCS
 * //Bond Sensitive is set True
 * //Ring Match is set True
 * Isomorphism comparison = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);
 * // set chemical filter true
 * comparison.setChemFilters(true, true, true);
 * //Get similarity score
 * System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
 * Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
 * Assert.assertEquals(12, comparison.getAllAtomMapping().size());
 * // Print the mapping between molecules
 * System.out.println(" Mappings: ");
 * for (AtomAtomMapping atomatomMapping : comparison.getAllAtomMapping()) {
 *      for (Map.Entry<IAtom, IAtom> mapping : atomatomMapping.getMappingsByAtoms().entrySet()) {
 *          IAtom sourceAtom = mapping.getKey();
 *          IAtom targetAtom = mapping.getValue();
 *          System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
 *          System.out.println(atomatomMapping.getQueryIndex(sourceAtom) + " " + atomatomMapping.getTargetIndex(targetAtom));
 *      }
 *      System.out.println("");
 *  }
 *
 *
 *  </pre> </font>
 *
 * @cdk.require java1.5+
 *
 * @cdk.module smsd
 * @cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
//@TestClass("org.openscience.cdk.smsd.factory.SubStructureSearchAlgorithmsTest")
public class ExtendedIsomorphism implements Serializable {
 
    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(ExtendedIsomorphism.class);
    static final long serialVersionUID = 0x24845e5c5ae878L;
    private ExtendedAlgorithm algorithmType;
    private double bondSensitiveMcGregorOut = -1;//mins
    private double bondInSensitiveMcGregor = -1;//mins
    
    private List<Map<Integer, Integer>> allAtomIndexMCS = null;
    private List<List<int[]>> allBondIndexMCS = null;
    private List<AtomAtomMapping> allAtomMCS = null;
    private List<Map<IBond, IBond>> allBondMCS = null;
    //private List<Map<IAtom, IAtom>> allAtomMCS = null;
    
    private boolean matchBonds = true;
    private boolean matchRings = false;
    private boolean cliqueDetectionRaymondHeuristics = false;
    private boolean cliqueDetectionRingHeuristics = false;
    private int topologicalDistanceLimit = -1;
    private int expansionTimeOut = 10000;
    
    
    private int mcsSize = 0;
    private int[] fragmentSizes = null;
    private String mcsSMARTS = null;
    
    private int modProdNodeCount = 0;
    private long modProdConstructionTime = 0;
    private double modProdEdgeDensity = 0.0;
    private long mcsTime = 0;
    
    private IAtomContainer mol1, mol2;

    /**
     * Initialize query and target molecules.
     *
     * Note: Here its assumed that hydrogens are implicit and user has called
     * these two methods percieveAtomTypesAndConfigureAtoms and
     * CDKAromicityDetector before initializing calling this method.
     *
     * @param query query molecule
     * @param target target molecule This is the algorithm factory and entry
     * port for all the MCS algorithm in the Isomorphism supported algorithm
     * {@link org.openscience.cdk.smsd.interfaces.Algorithm} types: <OL> <lI>0:
     * Default,
     * <lI>1: MCSPlus, <lI>2: VFLibMCS, <lI>3: CDKMCS </OL>
     * @param algorithmType
     * {@link org.openscience.cdk.smsd.interfaces.Algorithm}
     */
    //@TestMethod("testIsomorphismTest")
    public ExtendedIsomorphism(
            IQueryAtomContainer query,
            IAtomContainer target,
            ExtendedAlgorithm algorithmType) {
        
    	this.mol1 = query;
        this.mol2 = target;
        
        this.algorithmType = algorithmType;
        mcsBuilder(query, target);
        //setSubgraph(isSubgraph());
    }

    /**
     * Initialize query and target molecules.
     *
     * Note: Here its assumed that hydrogens are implicit and user has called
     * these two methods percieveAtomTypesAndConfigureAtoms and
     * CDKAromicityDetector before initializing calling this method.
     *
     * @param query query mol
     * @param target target mol This is the algorithm factory and entry port for
     * all the MCS algorithm in the Isomorphism supported algorithm
     * {@link org.openscience.cdk.smsd.interfaces.Algorithm} types: <OL> <lI>0:
     * Default,
     * <lI>1: MCSPlus, <lI>2: VFLibMCS, <lI>3: CDKMCS </OL>
     * @param algorithmType
     * {@link org.openscience.cdk.smsd.interfaces.Algorithm}
     * @param bondTypeFlag Match bond types (i.e. double to double etc)
     * @param matchRings Match ring atoms and ring size
     * @param matchAtomType
     */
    //@TestMethod("testIsomorphismTest")
    public ExtendedIsomorphism(
            IAtomContainer query,
            IAtomContainer target,
            ExtendedAlgorithm algorithmType,
            boolean bondTypeFlag,
            boolean matchRings,
            boolean matchAtomType,
            boolean rHeuristics,
            boolean ringHeuristics,
            int topoDistance,
            int timeLim
    		) {
        //super(query, target, bondTypeFlag, matchRings, matchAtomType);
        this.algorithmType = algorithmType;
        
        setUseRaymondHeuristics(rHeuristics);
        setUseRingHeuristics(ringHeuristics);
        setTopologicalDistanceLimit(topoDistance);
        setTimeLimit(timeLim);
        
        this.mol1 = query;
        this.mol2 = target;
        
       /* if (isMatchRings()) {
            try {
                MoleculeInitializer.initializeMolecule(getQuery());
                MoleculeInitializer.initializeMolecule(getTarget());
            } catch (CDKException ex) {
            }
        }*/
        mcsBuilder(query, target);
        //setSubgraph(isSubgraph());
    }

    
    private synchronized void mcsBuilder(IAtomContainer mol1, IAtomContainer mol2) {
    	
    	allAtomIndexMCS = Collections.synchronizedList(new ArrayList<Map<Integer, Integer>>());
        allBondIndexMCS =  Collections.synchronizedList(new ArrayList<List<int[]>>());
        allAtomMCS = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        allBondMCS = Collections.synchronizedList(new ArrayList<Map<IBond, IBond>>());
    	
        int rBondCount = mol1.getBondCount();
        int pBondCount = mol2.getBondCount();

        int rAtomCount = mol1.getAtomCount();
        int pAtomCount = mol2.getAtomCount();

        long mcsTimeTemp = System.currentTimeMillis();
        if ((rBondCount == 0 && rAtomCount > 0) || (pBondCount == 0 && pAtomCount > 0)) {
            singleMapping();
        } else {
            chooseAlgorithm();
        }
        
        if( mcsTime == 0 ) {
        	mcsTime = System.currentTimeMillis() - mcsTimeTemp;
        }
        
        /*if( allAtomMCS.isEmpty() ) {
        	allAtomMCS = getAllAtomMapping();
        }*/
        
        // if no bond maps were made
        if( ! allAtomMCS.isEmpty() && allBondMCS.isEmpty() ) {
        	
        	
        	List<List<Integer>> atomIndexChromosomes = new ArrayList<>();
        	for( AtomAtomMapping aam : allAtomMCS ) {
        		atomIndexChromosomes.add( MCSMethods.atomMapToChromosome(mol1, mol2, aam.getMappingsByAtoms() ) );
        	}
        	
        	allAtomMCS.clear();
        	mapChromosomesToMainMaps(atomIndexChromosomes, mol1, mol2);
        	
        	allBondMCS.addAll( makeBondMapsOfAtomMaps(mol1, mol2 , allAtomMCS) );
            allBondIndexMCS.addAll( bondMapsToIndices(allBondMCS, false) );
            
            mcsSize = allBondMCS.get(0).size();
            fragmentSizes = new int[]{ allAtomMCS.get(0).getCount() };
            mcsSMARTS = "not implemented";  // FIXME  temporarily null
        }
    }




    private synchronized void chooseAlgorithm() {

        switch (algorithmType) {
            case CDKMCS:
            case DEFAULT:
            case MCSPlus:
            case VFLibMCS:
                SMSDMCSAlgorithm(algorithmType);
                break;
            case ChemAxon_cMCES:
         	   chemAxon_MCES(true);
         	   break;
            case ChemAxon_dMCES:
         	   chemAxon_MCES(false);
         	   break;
            case BK_dMCES:
            case CP_dMCES:
            case RASCAL_dMCES:
            case Depolli_dMCES:
            case fMCS:
            case kCombu_dMCES:
            case kCombu_cMCES:
            case consR_dMCES:
            default:
         	   alternativeMCES();
         	   //System.out.println("Algorithm selected");
         	   break;
        }
    }
    
    private synchronized void SMSDMCSAlgorithm( ExtendedAlgorithm exAlg ) {
    	Algorithm alg = Algorithm.valueOf( exAlg.name() );
    	Isomorphism isomorphism = new org.openscience.smsd.Isomorphism(mol1, mol2, alg, matchBonds, matchRings, true);
    	clearMaps();
    	allAtomMCS = new ArrayList<>();
    	allAtomMCS.addAll( isomorphism.getAllAtomMapping() );
    	
    }



    private synchronized void singleMapping() {
        SingleMappingHandler mcs;
        mcs = new SingleMappingHandler( mol1, mol2, matchRings );
        clearMaps();
        allAtomMCS.addAll(mcs.getAllAtomMapping());
    }

    
    private synchronized void alternativeMCES() {
 	   
 	   MCSMethods mcs = null;
 	   CliqueDetection.ModularProductOptions mpOpts = new CliqueDetection.ModularProductOptions(cliqueDetectionRaymondHeuristics, cliqueDetectionRingHeuristics, topologicalDistanceLimit);
 	   
 	   switch (algorithmType) {
 	   	case BK_dMCES:
 	   		mcs = new BronKerboschCazals( mpOpts );
 	   		break;
 	    case BK_cMCES:
	   		mcs = new BronKerboschHariharan(mpOpts);
	   		break;
 	   	case CP_dMCES:
 	   		mcs = new CarraghanPardalosMCS( mpOpts );
 	   		break;
 	   	case RASCAL_dMCES:
 	   		mcs = new RASCALCliqueMCS( mpOpts, true );
 	   		break;
 	   	case Depolli_dMCES:
 	   		mcs = new DepolliCliqueDetection( mpOpts );
 	   		break;
 	   	case fMCS:
 	   		mcs = new FMCS();
 	   		break;
 	   	case kCombu_dMCES:
 	   		mcs = new KawabataBuildupMCS( topologicalDistanceLimit, false, cliqueDetectionRaymondHeuristics );  // dMCES
 	   		break;
 	   	case kCombu_cMCES:
 	   		mcs = new KawabataBuildupMCS( topologicalDistanceLimit, true, cliqueDetectionRaymondHeuristics );  // cMCES
 	   		break;
 	    case consR_dMCES:
	   		mcs = new ZhuSpectralMCES();   
	   		break;
 	   }
 	   
 	 
 	   mcs.setExpansionTimeLimit(expansionTimeOut);
 	   mcs.setMatchBonds( matchBonds );
 		
        
 	   mcs.setMainMol( mol1 );
 	   mcs.setQueryMol( mol2 );
        
 	   clearMaps();
 	   
        //mcs.setSmartsHandling(true);
        mcsTime = System.currentTimeMillis();
        mcs.execute();
        mcsTime = (System.currentTimeMillis() - mcsTime);
        //System.out.println( "ei mcs time - " + (System.currentTimeMillis() - mcsTime) );
        //mcs.searchMCS(isMatchBonds());
        
        
        
        modProdConstructionTime = mcs.getModularProductConstructionTime();
        modProdEdgeDensity = mcs.getModularProductEdgeDensity();
        modProdNodeCount = mcs.getModularProductNodeCount();
        //mcsTime = mcs.getMCSSearchTime();
        
         
        mapChromosomesToMainMaps( mcs.getBestAtomIndexMatches(), mol1, mol2 );
        
        // translate and add all maps
        //mapChromosomesToMainMaps(mcs, mcs.getMainMol(), mcs.getQueryMol() );

       // firstBondMCS = mcs.getBestBondMatches().get(0);
        allBondMCS = mcs.getBestBondMatches();
        allBondIndexMCS = mcs.getBestBondIndexMatches(); 
     		   
        //allBondIndexMCS = bondMapsToIndices(allBondMCS, false);
        calculateSubgraphInformation( mcs );
        
        
        
        //mcsSize = mcs.mcsSize;
        //fragmentSizes = mcs.fragmentSizes;
        //mcsSMARTS = mcs.mcsSMARTS;
    }
    
    
    private void chemAxon_MCES( boolean connectedMode ) {
 	   
 	   ChemAxonMCS mcs = null;
 	   ChemAxonMCSOptions mcsOpts = new ChemAxonMCSOptions();
 		mcsOpts.connectedMode = connectedMode;
 		mcsOpts.matchBonds = matchBonds;
 		mcsOpts.SMARTSHandling = true;
 		mcsOpts.verbose = false;
 		mcsOpts.ringEnforcement = true;
 		
 		mcs = new ChemAxonMCS( mol1, mol2, mcsOpts );
        //mcs.setSmartsHandling(true);
        
 		
        
        //mcsTime = System.currentTimeMillis();
        mcs.execute();
       // mcsTime = (System.currentTimeMillis() - mcsTime);
        
        clearMaps();
        
        
        mcsTime = mcs.getMCSSearchTime();
        //mcs.searchMCS(isMatchBonds());

        
        
        // translate and add all maps
        //mapChromosomesToMainMaps(mcs.getBestAtomIndexMatches(), mcs.getQueryMol(), mcs.getMainMol() );
        
        mapChromosomesToMainMaps( mcs.getBestAtomIndexMatches(), mol1, mol2 );
        
       // firstBondMCS = mcs.getBestBondMatches().get(0);
        allBondMCS = mcs.getBestBondMatches();
        //System.out.println( allBondMCS );
        
        allBondIndexMCS = bondMapsToIndices(allBondMCS, false);
        //calculateSubgraphInformation( allBondIndexMCS.get(0) );
        
        
        //firstBondMCS = makeBondMapOfAtomMap( mcs.getQueryMol(), mcs.getMainMol(), firstAtomMCS);
        //allBondMCS = makeBondMapsOfAtomMaps(mcs.getQueryMol(), mcs.getMainMol(), allAtomMCS);
        calculateSubgraphInformation( mcs );
        
    }
    
    
    
    private void mapChromosomesToMainMaps( List<List<Integer>> atomIndexChromosomeMaps, IAtomContainer mainMol, IAtomContainer queryMol ) {
 	   int sCount = 0;
 	   
 	   for( List<Integer> mapChr : atomIndexChromosomeMaps ) {
 	    	// instantiate map variables
 	    	   Map<Integer, Integer> map = new HashMap<>();
 	    	   //HashMap<IAtom, IAtom> atomMap = new HashMap<IAtom, IAtom>();
 	    	  AtomAtomMapping atomMap = new AtomAtomMapping( mainMol, queryMol );
 	    	   
 	    	   //System.out.println("stuff: " + mapChr.size() + " " + queryMol.getAtomCount() + " " + mainMol.getAtomCount() + " | " + mapChr );
 	    	   
 	    	   // translate mapping from array to hash
 	    	   for( int m = 0; m < mapChr.size(); m++ ) {
 	    		   if( mapChr.get(m) >= 0 ) {
 	    			   map.put( m,  mapChr.get(m)  );
 	    			   
 	    			   if( mainMol.getAtom( m ) != null && queryMol.getAtom(  mapChr.get(m) ) != null )
 	    				   atomMap.put( mainMol.getAtom( m ), queryMol.getAtom(  mapChr.get(m) ) );
 	    		   }
 	    	   }
 	    	   
 	    	   
 	    	   // refuse to deal with over a certain amount of matches
 	    	   if( ++sCount > 1000 ) {
 	    		   break;
 	    	   }
 	    	   
 	    	   allAtomIndexMCS.add(map);
 	    	   allAtomMCS.add(atomMap);
 	       }
    }
    
    
    
    
    
    public int getCAMCSSize() {
 	   return mcsSize;
    }
    
    
    
    /**
     * Several methods (i.e. the pre-existing SMSD MCS methods) only find atom mappings, not bond mappings
     * 
     * this converts a list of bond maps to bond index (integers) maps
     * 
     * 
     * @param bondMaps
     * @param reverse
     * @return
     */
    private List<List<int[]>> bondMapsToIndices( List<Map<IBond, IBond>> bondMaps, boolean reverse ) {
 	   ArrayList<List<int[]>> indicesList = new ArrayList<List<int[]>>( );
 	   
 	   
 	   IAtomContainer ac1, ac2;
 	   
 	   ac1 = mol1;
	   ac2 = mol2;
	   
 	   for( Map<IBond, IBond> bondMap : bondMaps ) {
 		   List<int[]> indices = new ArrayList<int[]>( bondMap.size() );
 		   
 		   for( Entry<IBond, IBond> entry : bondMap.entrySet() ) {
 			   if( reverse )
 				   indices.add( new int[]{ 
 						   ac1.getBondNumber( entry.getValue() ) , 
 						   ac2.getBondNumber( entry.getKey() ) 
 				   } );
 			   else
 				   indices.add( new int[]{
 						   ac1.getBondNumber( entry.getKey() ) , 
 						   ac2.getBondNumber( entry.getValue() )  
 				   });
 		   }
 		   
 		   indicesList.add(indices);
 	   }
 	   
 	   return indicesList;
    }

    
    /**
     * Determine miscellaneous information regarding an MCS
     * 
     * @param bondMapIndices
     */
    private void calculateSubgraphInformation( MCSMethods mapper ) {
 	   mcsSize = mapper.mcsSize;
 	   fragmentSizes = mapper.fragmentSizes;
 	   mcsSMARTS = mapper.mcsSMARTS;
    }
    
    
    /**
     *
     * @return true if query is a subgraph of the target
     */
    //@TestMethod("testIsSubgraph")
    public synchronized boolean isSubgraph() {

        float mappingSize;
        if (getMappingCount() > 0) {
            mappingSize = getAllAtomMapping().iterator().next().getCount();
        } else {
            return false;
        }
        int sourceAtomCount = mol1.getAtomCount();
        int targetAtomCount = mol2.getAtomCount();

        if (mappingSize == sourceAtomCount && mappingSize <= targetAtomCount) {
            if (mappingSize == 1) {
                return true;
            } else if (!getAllBondMaps().isEmpty()
                    && getAllBondMaps().iterator().next().size() == mol1.getBondCount()) {
                return true;
            }
        }
        return false;
    }
    
    
    /**
     * Returns bond maps between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     *
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mappings mappings between sourceAtomCount and targetAtomCount
     * molecule atoms
     * @return bond maps between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     */
    public synchronized List<Map<IBond, IBond>> makeBondMapsOfAtomMaps(IAtomContainer ac1,
            IAtomContainer ac2, List<AtomAtomMapping> mappings) {
        List<Map<IBond, IBond>> bondMaps = Collections.synchronizedList(new ArrayList<Map<IBond, IBond>>());
        for (AtomAtomMapping mapping : mappings) {
            bondMaps.add(makeBondMapOfAtomMap(ac1, ac2, mapping));
        }
        return bondMaps;
    }

    /**
     *
     * Returns bond map between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     *
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mapping mappings between sourceAtomCount and targetAtomCount
     * molecule atoms
     * @return bond map between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     */
    private synchronized Map<IBond, IBond> makeBondMapOfAtomMap(IAtomContainer ac1, IAtomContainer ac2,
            AtomAtomMapping mapping) {

        Map<IBond, IBond> bondbondMappingMap = Collections.synchronizedMap(new HashMap<IBond, IBond>());

        for (Map.Entry<IAtom, IAtom> map1 : mapping.getMappingsByAtoms().entrySet()) {
            for (Map.Entry<IAtom, IAtom> map2 : mapping.getMappingsByAtoms().entrySet()) {
                if (map1.getKey() != map2.getKey()) {
                    IBond bond1 = ac1.getBond(map1.getKey(), map2.getKey());
                    IBond bond2 = ac2.getBond(map1.getValue(), map2.getValue());
                    if (bond1 != null && bond2 != null && !bondbondMappingMap.containsKey(bond1)) {
                        bondbondMappingMap.put(bond1, bond2);
                    }
                }
            }
        }
//        System.out.println("Mol Map size:" + bondbondMappingMap.size());
        return bondbondMappingMap;
    }
    
    
    
    /**
     * Utilise node and edge-deletion heuristics for modular product (clique detection algorithms only)
     * 
     * @param use
     */
    public void setUseRaymondHeuristics( boolean use ) {
 	   cliqueDetectionRaymondHeuristics = use;
    }
    
    /**
     * Weak ring and Strong ring edge-deletion heuristics for modular product (clique detection algorithms only)
     * 
     * @param use
     */
    public void setUseRingHeuristics( boolean use ) {
 	   cliqueDetectionRingHeuristics = use;
    }
    
    
    /**
     * Set topological distance limit for tdMCS (clique detection and build-up algorithms only)
     * 
     * @param tdl
     */
    public void setTopologicalDistanceLimit( int tdl ) {
 	   topologicalDistanceLimit = tdl;
    }
    
    
    
    public ExtendedAlgorithm getAlgorithmType() {
		return algorithmType;
	}

	public boolean isMatchBonds() {
		return matchBonds;
	}

	public boolean isMatchRings() {
		return matchRings;
	}

	public boolean isCliqueDetectionRaymondHeuristics() {
		return cliqueDetectionRaymondHeuristics;
	}

	public boolean isCliqueDetectionRingHeuristics() {
		return cliqueDetectionRingHeuristics;
	}

	public int getTopologicalDistanceLimit() {
		return topologicalDistanceLimit;
	}

	public int getExpansionTimeOut() {
		return expansionTimeOut;
	}

	public void setTimeLimit( int tl ) {
 	   expansionTimeOut = tl;
    }

    /**
     * @return the bondSensitiveMcGregorOut
     */
    public double getBondSensitiveMcGregorOut() {
        return bondSensitiveMcGregorOut;
    }

    /**
     * @param bondSensitiveMcGregorOut the bondSensitiveMcGregorOut to set
     */
    public void setBondSenSitiveMcGregorOut(double bondSensitiveMcGregorOut) {
        this.bondSensitiveMcGregorOut = bondSensitiveMcGregorOut;
    }

    /**
     * @return the bondInSensitiveMcGregor
     */
    public double getBondInSensitiveMcGregor() {
        return bondInSensitiveMcGregor;
    }

    /**
     * @param bondInSensitiveMcGregor the bondInSensitiveMcGregor to set
     */
    public void setBondInSenSitiveMcGregor(double bondInSensitiveMcGregor) {
        this.bondInSensitiveMcGregor = bondInSensitiveMcGregor;
    }
    
    
    
    
    public synchronized void clearMaps() {
    	this.allAtomMCS.clear();
    	this.allAtomIndexMCS.clear();
    	this.allBondMCS.clear();
    	this.allBondIndexMCS.clear();
    	
    	fragmentSizes = null;
    	mcsSMARTS = null;
    	
    	modProdConstructionTime = 0;
    	modProdEdgeDensity = 0.0;
    	modProdNodeCount = 0;
    	mcsTime = 0;
    	mcsSize = 0;
    }

    public synchronized int getMappingCount() {
        return allBondMCS.isEmpty() ? 0 : allBondMCS.size();
    }
    
    
    /**
     * @return the allBondMCS
     */

    public synchronized List<Map<IBond, IBond>> getAllBondMaps() {
    	
    	if( ! allBondMCS.isEmpty() )
    		return allBondMCS;
    	
        if (!allAtomMCS.isEmpty()) {
        	allBondMCS.addAll( makeBondMapsOfAtomMaps(mol1, mol2, allAtomMCS) );
        	return allBondMCS;
        }
        return new ArrayList<>();
    }
    
    /**
     * @return the allBondIndexMCS
     */
    public List<List<int[]>> getAllBondIndexMaps() {
        return allBondIndexMCS;
    }
    
    
    
    /**
     * {@inheritDoc}
     *
     * @return
     */
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        //return Collections.unmodifiableList(new ArrayList<>(getMCSList()));
    	return allAtomMCS;
    }
    
    /** {@inheritDoc}
     */
    public synchronized List<Map<Integer, Integer>> getAllAtomIndexMapping() {
        return allAtomIndexMCS.isEmpty() ? null : allAtomIndexMCS;
    }
    

    public synchronized int[] getFragmentSizes() {
        return (fragmentSizes != null && fragmentSizes.length != 0 )
                ? fragmentSizes : null;
    }
    
    public String getMCSSMARTS() {
 	   if( mcsSMARTS == null )
 		   return "Error - MCS SMARTS is non-existant";
 	   
 	   return mcsSMARTS;
    }

    
    public long getModularProductConstructionTime() {
 		return modProdConstructionTime;
 	}
 	
 	public int getModularProductNodeCount() {
 		return modProdNodeCount;
 	}
 	
 	public double getModularProductEdgeDensity() {
 		return modProdEdgeDensity;
 	}
 	
 	public long getElapsedTime() {
 	   return mcsTime;
    }
}
